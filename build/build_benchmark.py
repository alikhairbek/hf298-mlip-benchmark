#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_benchmark.py — verified experimental gas-phase ΔHf°(298 K) test set from the Bains et al. (2022) compilation
(Zenodo 10.5281/zenodo.4661783, Measured_Enthalpy_V2.7.xlsx; values in kJ/mol from up to 24 source columns).

  python build_benchmark.py Measured_Enthalpy_V2.7.xlsx audit_index/ stage1.csv      # then: python apply_atct_and_classes.py stage1.csv ../benchmark_verified.csv

Rules
  * HCNOF only; RDKit-parsable; net-neutral; closed-shell (no radical electrons); duplicates merged by full InChIKey (sources pooled)
  * consensus value = MEDIAN of all source values (robust to a single wrong column — e.g. the compilation's "ATcT" column is
    off by 4–5 kcal/mol for formic and acetic acid); each source deviating > 2 kcal/mol from the median is counted as an outlier
  * uncertainty proxy = half the inter-source spread (NaN for single-source entries); n_sources and outlier sources recorded
  * robust spread: with ≥ 3 sources, sources deviating > 2 kcal/mol from the median are set aside (listed in outlier_sources)
  * tiers:  A = ≥ 2 agreeing sources and robust spread ≤ 1 kcal/mol   B = single source, or robust spread ≤ 2 kcal/mol
            C = robust spread > 2 kcal/mol (kept, flagged)              X = spread > 12 kcal/mol with no majority → excluded
  * stereochemistry: if the SMILES leaves stereo elements unspecified, all stereoisomers are enumerated (≤ 16), embedded and
    MMFF94-optimised; when ≥ 2 distinct, embeddable DIASTEREOMERS exist, the lowest-MMFF-energy one is assigned
    (column smiles = assigned isomer, stereo_assumed = True) so that a sensitivity analysis can exclude these molecules
  * strata: n_heavy ≤ 9 (QM9-size) vs ≥ 10 (size extrapolation)
  * leakage-audit columns in_<set>: InChIKey-skeleton match against audit_index/<set>.csv
"""
import pandas as pd, numpy as np, sys
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi, AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
RDLogger.DisableLog("rdApp.*")
KJ = 1/4.184
SRC = ['Add Lit 1','Add Lit 2','Add Lit 3','Add Lit 4','Add Lit 5','Pedley','Cox-1','Cox-2','Cox-3','Cox-4','Cox-5','Winget','Stewart',
       'Jorgensen','ATcT','Sandler','Cioslowski','Yaws','NIST-1','NIST-2','NIST-3','NIST-4','NIST-5','JANAF']

def mirror(m):
    """invert every tetrahedral centre → enantiomer (double-bond geometry unchanged)"""
    mm = Chem.Mol(m)
    for a in mm.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return Chem.MolToSmiles(mm)

def resolve_stereo(m):
    """returns (assigned_mol, n_diastereomers, assumed_flag)"""
    si = Chem.FindPotentialStereo(m)
    if not any(x.specified == Chem.StereoSpecified.Unspecified for x in si):
        return m, 1, False
    opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=16, unique=True)
    cands = []
    for iso in EnumerateStereoisomers(m, options=opts):
        mh = Chem.AddHs(iso); ps = AllChem.ETKDGv3(); ps.randomSeed = 7; ps.enforceChirality = True
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=4, params=ps))
        if not cids: continue
        res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=2000)
        e = min(en for conv, en in res)
        cands.append((e, Chem.MolToSmiles(iso), iso))
    if not cands: return m, 1, False
    # distinct diastereomers: merge enantiomeric pairs
    seen = set(); dia = []
    for e, smi, iso in sorted(cands, key=lambda t: t[0]):
        if smi in seen: continue
        seen.add(smi); seen.add(mirror(iso)); dia.append((e, smi, iso))
    if len(dia) >= 2:
        return dia[0][2], len(dia), True
    return dia[0][2], 1, False

def main(xlsx, audit_dir, out_csv):
    df = pd.read_excel(xlsx, header=2); df.columns = [str(c).strip() for c in df.columns]
    for c in SRC: df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ['C','H','As','B','Br','Cl','F','Ge','I','N','O','P','S','Se','Si']: df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df = df[df[['As','B','Br','Cl','Ge','I','P','S','Se','Si']].sum(axis=1) == 0].copy()
    log = {"HCNOF rows": len(df)}
    rows = {}
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(str(r["SMILES"]))
        if m is None: continue
        if sum(a.GetFormalCharge() for a in m.GetAtoms()) != 0: continue
        if any(a.GetNumRadicalElectrons() > 0 for a in m.GetAtoms()): continue
        if not {a.GetSymbol() for a in m.GetAtoms()} <= {"C","H","N","O","F"}: continue
        key = inchi.MolToInchiKey(m); vals = [(c, float(r[c])) for c in SRC if pd.notna(r[c])]
        if not vals: continue
        if key in rows: rows[key]["vals"] += vals; rows[key]["names"].add(str(r["Name"]).strip()); continue
        rows[key] = dict(mol=m, names={str(r["Name"]).strip()}, vals=list(vals), bains_filter=(r.get("FILTER") == 1))
    log["parsable, neutral, closed-shell, unique"] = len(rows)
    idx = {p.name.split(".")[0]: set(pd.read_csv(p)["skeleton"]) for p in sorted(Path(audit_dir).glob("*.csv*"))} if audit_dir else {}
    out = []
    for key, d in rows.items():
        m0 = d["mol"]; v = np.array([x[1] for x in d["vals"]], float); srcs = [x[0] for x in d["vals"]]
        med = float(np.median(v)); spread = float(v.max() - v.min()) if len(v) > 1 else np.nan; sp_kcal = spread*KJ if len(v) > 1 else np.nan
        outl = sorted({s for s, x in d["vals"] if abs(x - med)*KJ > 2})
        atct = [x for s, x in d["vals"] if s == "ATcT"]
        # robust spread: when >= 3 sources, drop sources deviating > 2 kcal/mol from the median (documented in outlier_sources)
        vin = v[np.abs(v - med)*KJ <= 2] if len(v) >= 3 else v
        sp_rob = float(vin.max() - vin.min())*KJ if len(vin) > 1 else np.nan
        if len(v) > 1 and sp_kcal > 12 and len(vin) < 3: tier = "X"      # discordant, no majority
        elif len(vin) > 1 and sp_rob <= 1: tier = "A"
        elif len(v) == 1 or (len(vin) > 1 and sp_rob <= 2): tier = "B"
        else: tier = "C"
        m, n_dia, assumed = resolve_stereo(m0)
        counts = {e: sum(1 for a in Chem.AddHs(m).GetAtoms() if a.GetSymbol() == e) for e in "HCNOF"}
        n_heavy = m.GetNumHeavyAtoms(); sk = key[:14]
        out.append(dict(name=sorted(d["names"])[0], smiles=Chem.MolToSmiles(m), smiles_original=Chem.MolToSmiles(m0), inchikey=key,
                        dHf_exp_kcal=med*KJ, dHf_exp_kJ=med, unc_kcal=(sp_rob/2 if len(vin) > 1 else np.nan), n_sources=len(v),
                        sources=";".join(sorted(set(srcs))), spread_kcal=sp_kcal, spread_robust_kcal=sp_rob, n_sources_agreeing=len(vin), outlier_sources=";".join(outl),
                        has_ATcT=bool(atct), ATcT_minus_median_kcal=((atct[0] - med)*KJ if atct else np.nan), tier=tier,
                        bains_filter_flag=bool(d["bains_filter"]), n_diastereomers=n_dia, stereo_assumed=assumed,
                        n_heavy=n_heavy, stratum=("small(<=9)" if n_heavy <= 9 else "large(>=10)"), n_rot=Descriptors.NumRotatableBonds(m),
                        n_arom_rings=rdMolDescriptors.CalcNumAromaticRings(m), n_rings=rdMolDescriptors.CalcNumRings(m),
                        **{f"n_{e}": counts[e] for e in "HCNOF"}, **{f"in_{k}": (sk in s) for k, s in idx.items()}))
    B = pd.DataFrame(out).sort_values(["n_heavy", "name"]).reset_index(drop=True)
    B.insert(0, "id", [f"B{i+1:04d}" for i in range(len(B))])
    elemental = B.smiles.isin(["[H][H]", "[HH]", "N#N", "O=O", "FF"])          # reference states: ΔHf° ≡ 0 by definition
    multifrag = B.smiles.str.contains(r"\.", regex=True)                      # disconnected SMILES = ion pair written as two fragments (the compilation's
    B["evaluate"] = (B.tier != "X") & (~elemental) & (~multifrag)             # "pyridinium dicyanomethylide" is C8H7N3, the real ylide is C8H5N3)
    B["exclude_reason"] = np.where(B.tier == "X", "discordant sources", np.where(elemental, "reference state",
                                   np.where(multifrag, "multi-fragment SMILES (wrong formula)", "")))
    B.to_csv(out_csv, index=False)
    log.update({"tier A": int((B.tier=="A").sum()), "tier B": int((B.tier=="B").sum()), "tier C": int((B.tier=="C").sum()), "tier X (excluded)": int((B.tier=="X").sum()),
                "stereo assumed (kept, flagged)": int(B.stereo_assumed.sum()), "ATcT column deviating >2 kcal from median": int((B.ATcT_minus_median_kcal.abs() > 2).sum()),
                "EVALUATION SET": int(B.evaluate.sum()), "  n_heavy<=9": int((B.evaluate & (B.n_heavy<=9)).sum()), "  n_heavy>=10": int((B.evaluate & (B.n_heavy>=10)).sum()),
                "  >=2 sources": int((B.evaluate & (B.n_sources>=2)).sum()), "  tier A": int((B.evaluate & (B.tier=="A")).sum())})
    for k in idx: log[f"  in {k}"] = int((B.evaluate & B[f"in_{k}"]).sum())
    for k, val in log.items(): print(f"{k:45s} {val}")
    return B

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None, sys.argv[3] if len(sys.argv) > 3 else "benchmark_verified.csv")
