#!/usr/bin/env python3
"""Stage 2 of the benchmark build: ATcT 1.222 override, charge classes, formula, exclusions.

    python apply_atct_and_classes.py stage1.csv atct_v1222_gas_HCNOF_subset.psv atct_matches.csv ../benchmark_verified.csv

stage1.csv is the output of build_benchmark.py.  The ATcT subset (name|formula|dHf298_kJ|unc_kJ) holds every closed-shell gas-phase
HCNOF species of ATcT version 1.222 whose Hill formula occurs in the set (read from atct.anl.gov, 1 Oct 2026).

Matching (atct_matches.csv, 177 rows): a formula that is unique in both the ATcT subset and the benchmark ("unique-formula"); equality of
the normalised names ("name"); hand-vetted synonyms ("name-manual").  Isomer ambiguities were left unmatched.
For every match the ATcT value and uncertainty replace the compilation consensus and the molecule is graded A+; the compilation value
is kept in dHf_compilation_kcal.  Charge classes are assigned from the SMILES (formal charges in standard dipolar groups vs. other
charge-separated species).  Exclusions: multi-fragment SMILES and the nine structure–value mismatches of the reference verification.
"""
import sys, re, numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
KJ = 1/4.184

# The hand-vetted ATcT matching table (id, benchmark name, formula, ATcT name, how matched) ships with the repository as
# build/atct_matches.csv; it was produced by unique-formula matching, normalised-name matching and manual synonym vetting,
# and is applied verbatim so that the released benchmark is reproduced exactly.
def norm(s):
    s = str(s).lower().strip()
    s = re.sub(r"[\s\-–,'′()\[\]]", "", s)
    s = s.replace("ortho", "o").replace("para", "p").replace("meta", "m")
    return s

def hill(smiles):
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))

STD = ["[N+][O-]", "[n+][O-]", "[N-]=[N+]=N", "N=[N+]=[N-]", "[C-]#[N+]", "[C-]#[O+]"]   # nitro / N-oxide / nitronate, azide, isocyanide, CO
def charge_class(smiles):
    m = Chem.MolFromSmiles(smiles)
    charged = [a for a in m.GetAtoms() if a.GetFormalCharge() != 0]
    if not charged: return "neutral"
    covered = set()
    for p in STD:
        q = Chem.MolFromSmarts(p)
        for match in m.GetSubstructMatches(q): covered.update(match)
    if all(a.GetIdx() in covered for a in charged): return "standard dipolar group (nitro/azide/diazo/isocyanide/N-oxide/CO)"
    return "other charge-separated (ylide/zwitterion)"

MISMATCH = {"B0055": "structure–value mismatch (imidic-acid SMILES, formamide value)", "B0582": "structure–value mismatch (non-aromatic tautomer SMILES)",
            "B0925": "structure–value mismatch (value is for benzofurazan, SMILES is the 1,2,3-isomer)", "B0849": "structure–value mismatch (SMILES is an O-methyl oxime, not dimethylglyoxime; duplicate of B0710)",
            "B1345": "structure–value mismatch (7H-adenine SMILES; gas-phase species is 9H)", "B0646": "structure–value mismatch (SMILES is propyl acetate, value is isopropyl acetate)",
            "B1480": "structure–value mismatch (hydroperoxide placed on the ring)", "B0900": "ambiguous structure (which tetramethyldiazetidine isomer)", "B1109": "ambiguous stereo (value probably for the cis isomer)"}

def main(stage1, atct_psv, matches, out):
    B = pd.read_csv(stage1); A = pd.read_csv(atct_psv, sep="|")
    B["formula"] = [hill(s) for s in B.smiles]; B["key"] = B.inchikey.str[:14]
    A["formula"] = A.formula.str.replace(" ", "")
    B["atct_name"] = pd.Series([None]*len(B), dtype="object"); B["match"] = pd.Series([None]*len(B), dtype="object")
    B["atct_dHf_kcal"] = np.nan; B["atct_unc_kcal"] = np.nan
    T = pd.read_csv(matches).set_index("id"); A = A.set_index("name")
    for i, r in T.iterrows():
        j = B.index[B.id == i]
        if len(j) != 1 or r.atct_name not in A.index: raise SystemExit(f"cannot apply match {i} -> {r.atct_name}")
        a = A.loc[r.atct_name]; j = j[0]
        B.loc[j, "atct_name"] = r.atct_name; B.loc[j, "match"] = r["match"]; B.loc[j, "atct_dHf_kcal"] = a.dHf298_kJ*KJ; B.loc[j, "atct_unc_kcal"] = a.unc_kJ*KJ
    B["dHf_compilation_kcal"] = B.dHf_exp_kcal; B["unc_compilation_kcal"] = B.unc_kcal
    m = B.match.notna()
    B.loc[m, "dHf_exp_kcal"] = B.loc[m, "atct_dHf_kcal"]; B.loc[m, "dHf_exp_kJ"] = B.loc[m, "atct_dHf_kcal"]/KJ; B.loc[m, "unc_kcal"] = B.loc[m, "atct_unc_kcal"]
    B.loc[m, "tier"] = "A+"; B["dev_vs_ATcT_kcal"] = np.where(m, B.dHf_compilation_kcal - B.atct_dHf_kcal, np.nan)
    x = m & (B.exclude_reason.fillna("") == "discordant sources")          # an ATcT value supersedes a discordant compilation (tier X)
    B.loc[x, "evaluate"] = True; B.loc[x, "exclude_reason"] = ""
    B["value_source"] = np.where(m, "ATcT v1.222", "compilation median")
    B["charge_class"] = [charge_class(s) for s in B.smiles]
    for i, why in MISMATCH.items():
        B.loc[B.id == i, "evaluate"] = False; B.loc[B.id == i, "exclude_reason"] = why
    B.to_csv(out, index=False)
    print(f"matched ATcT: {int(m.sum())} ({B.match.value_counts().to_dict()}); evaluation set: {int(B.evaluate.sum())}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
