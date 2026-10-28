#!/usr/bin/env python3
"""Second-level analyses on the evaluated runs: correlation by tier, family consensus, ensembles, leakage (matched) tests, outliers."""
import numpy as np, pandas as pd
from pathlib import Path
OUT = Path(__file__).parent / "output" / "eval"
MODELS = ["uma-m", "uma-s", "mace-off23-medium", "mace-off23-small", "aimnet2", "egret-1", "egret-1e", "ani-2x", "ani-1ccx", "xtb-gfn2"]
B = pd.read_csv(str(Path(__file__).resolve().parent.parent / "benchmark_verified.csv")).set_index("id")
E = pd.read_csv(OUT / "errors_cv_wide.csv", index_col=0)
ev = {m: pd.read_csv(OUT / f"{m}_evaluated.csv").set_index("id") for m in MODELS}
NB = 10000; rng = np.random.default_rng(1)
def boot_mean(x):
    x = np.asarray(x, float); idx = rng.integers(0, len(x), (NB, len(x))); b = x[idx].mean(1); return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)
pd.set_option("display.width", 250); pd.set_option("display.float_format", "{:.2f}".format)

# 1. error correlations restricted to the ATcT-verified / verified tiers (shared experimental error is smallest there)
lines = []
for label, mask in [("all", B.evaluate), ("tier A+", B.tier == "A+"), ("tiers A+ and A", B.tier.isin(["A+", "A"])), ("tier B", B.tier == "B")]:
    ids = E.index[E.index.isin(B.index[mask])]
    sub = E.loc[ids].drop(columns=["uma-s", "mace-off23-small"])          # one representative per model size pair (Egret-1 and -1e are different models, both kept)
    C = sub.dropna().corr()
    lines.append(f"\n=== Pearson r of CV errors — {label} (n = {len(sub.dropna())}) ===\n" + C.round(2).to_string())
    C.to_csv(OUT / f"error_correlation_{label.replace(' ', '_').replace('+', 'plus')}.csv")
print("".join(lines))

# 2. family consensus: the ωB97M family (UMA-m, AIMNet2, MACE-OFF23 medium), ANI-1ccx (CCSD(T)*/CBS, HCNO) and ANI-2x (ωB97X/6-31G*)
fam = E[["uma-m", "aimnet2", "mace-off23-medium"]].median(1).rename("wB97M_family_median")
F = pd.concat([fam, E["ani-1ccx"], E["ani-2x"]], axis=1)
thr = 3.0
pos = (F > thr); neg = (F < -thr)
avail = F.notna().sum(1)
agree = ((pos.sum(1) == avail) | (neg.sum(1) == avail)) & (avail >= 2)
sus = B.loc[F.index[agree], ["name", "formula", "n_heavy", "dHf_exp_kcal", "unc_kcal", "tier", "value_source", "n_sources", "sources", "charge_class", "stereo_assumed"]].copy()
sus = pd.concat([sus, F.loc[sus.index]], axis=1)
sus["families_available"] = avail[sus.index]
sus = sus.reindex(sus.wB97M_family_median.abs().sort_values(ascending=False).index)
sus.to_csv(OUT / "family_consensus_suspect_references.csv")
print(f"\n=== reference values contradicted by every available model family (> {thr} kcal/mol, same sign): {len(sus)} molecules "
      f"(A+: {(sus.tier=='A+').sum()}, A: {(sus.tier=='A').sum()}, B: {(sus.tier=='B').sum()}, C: {(sus.tier=='C').sum()}; "
      f"single-source: {(sus.n_sources==1).sum()}) ===")
print(sus[["name", "formula", "dHf_exp_kcal", "tier", "n_sources", "wB97M_family_median", "ani-1ccx", "ani-2x"]].head(40).to_string())

# 3. ensembles of weakly correlated families
P = {m: ev[m].pred_cv for m in MODELS}; y = B.dHf_exp_kcal
def mae_ci(pred, ids):
    e = (pred.loc[ids] - y.loc[ids]).values; m, lo, hi = boot_mean(np.abs(e)); return m, lo, hi, len(ids)
ids_hcno = E[["uma-m", "ani-1ccx"]].dropna().index
res = []
for name, pred in [("UMA-m", P["uma-m"]), ("ANI-1ccx", P["ani-1ccx"]), ("mean(UMA-m, ANI-1ccx)", (P["uma-m"] + P["ani-1ccx"]) / 2),
                   ("mean(UMA-m, AIMNet2)", (P["uma-m"] + P["aimnet2"]) / 2), ("mean(UMA-m, MACE-OFF23 medium)", (P["uma-m"] + P["mace-off23-medium"]) / 2),
                   ("mean(UMA-m, ANI-1ccx, AIMNet2)", (P["uma-m"] + P["ani-1ccx"] + P["aimnet2"]) / 3)]:
    m, lo, hi, n = mae_ci(pred, ids_hcno); res.append(dict(predictor=name, n=n, MAE=m, CI_lo=lo, CI_hi=hi))
ENS = pd.DataFrame(res); ENS.to_csv(OUT / "ensembles_HCNO.csv", index=False); print("\n=== ensembles on the HCNO set ===\n" + ENS.to_string(index=False))

# 4. leakage: matched paired test.  For MACE-OFF23 (trained on SPICE) the per-molecule |error| difference to UMA-m (not trained on SPICE)
#    is compared between molecules inside and outside the SPICE v1 inputs; a training-set advantage would show as a more negative difference inside.
rows = []
def paired_by(mask_name, mask, a, b):
    d = (E[a].abs() - E[b].abs()).dropna(); m = mask.reindex(d.index).fillna(False).astype(bool)
    for lab, sel in [(f"in {mask_name}", m), (f"not in {mask_name}", ~m)]:
        mu, lo, hi = boot_mean(d[sel]); rows.append(dict(test=f"|err {a}| − |err {b}|", subset=lab, n=int(sel.sum()), mean_diff=mu, CI_lo=lo, CI_hi=hi))
paired_by("SPICE v1", B.in_spice_v1_inputs, "mace-off23-medium", "uma-m")
paired_by("SPICE v1", B.in_spice_v1_inputs, "mace-off23-small", "uma-m")
paired_by("SPICE v1", B.in_spice_v1_inputs, "mace-off23-medium", "aimnet2")
paired_by("SPICE v1", B.in_spice_v1_inputs, "egret-1", "uma-m")
paired_by("QM9", B.in_qm9, "mace-off23-medium", "uma-m")
# size-matched comparison inside/outside SPICE for each model on its own: restrict to <= 9 heavy atoms, tiers A+/A
m_small = (B.n_heavy <= 9) & B.tier.isin(["A+", "A"])
for mod in ["mace-off23-medium", "mace-off23-small", "egret-1", "uma-m", "aimnet2", "ani-2x", "ani-1ccx"]:
    e = E[mod].dropna(); sel = m_small.reindex(e.index).fillna(False)
    for lab, msk in [("in SPICE v1", B.in_spice_v1_inputs), ("not in SPICE v1", ~B.in_spice_v1_inputs)]:
        ids = e.index[sel & msk.reindex(e.index).fillna(False)]
        mu, lo, hi = boot_mean(np.abs(e.loc[ids])); rows.append(dict(test=f"MAE {mod} (≤9 heavy, tiers A+/A)", subset=lab, n=len(ids), mean_diff=mu, CI_lo=lo, CI_hi=hi))
LK = pd.DataFrame(rows); LK.to_csv(OUT / "leakage_tests.csv", index=False); print("\n=== leakage tests ===\n" + LK.to_string(index=False))

# 5. largest |errors| per model (CV), with the other models' errors alongside
for mod in ["uma-m", "mace-off23-medium", "aimnet2", "egret-1", "egret-1e", "ani-1ccx", "ani-2x"]:
    w = E[mod].abs().sort_values(ascending=False).head(15).index
    T = B.loc[w, ["name", "formula", "dHf_exp_kcal", "tier", "n_sources"]].copy()
    for m2 in ["uma-m", "mace-off23-medium", "aimnet2", "egret-1", "egret-1e", "ani-1ccx", "ani-2x"]: T[m2] = E.loc[w, m2]
    T.to_csv(OUT / f"top15_errors_{mod}.csv"); print(f"\n=== 15 largest |CV errors| of {mod} ===\n" + T.to_string())

# 6. vibrational consistency: ZPVE relative to the median of UMA-m / MACE-OFF23 medium / AIMNet2, and the "energy-only" check for Egret-1e
import json
load = lambda p: {json.loads(l)["id"]: json.loads(l) for l in open(p)}
J = {m: load(Path(__file__).resolve().parent.parent / "raw_results" / f"{m}.jsonl") for m in MODELS}
ok_ids = [i for i in B.index[B.evaluate] if all(i in J[m] and "zpve" in J[m][i] for m in ["uma-m", "mace-off23-medium", "aimnet2"])]
ref = pd.Series({i: np.median([J[m][i]["zpve"] for m in ["uma-m", "mace-off23-medium", "aimnet2"]]) for i in ok_ids})
zrows = []
for m in MODELS:
    r = pd.Series({i: J[m][i]["zpve"] / ref[i] for i in ok_ids if i in J[m] and "zpve" in J[m][i]})
    zrows.append(dict(model=m, n=len(r), median_ratio=round(r.median(), 4), pct_gt2=round(100*np.mean(np.abs(r-1) > 0.02), 1),
                      pct_gt5=round(100*np.mean(np.abs(r-1) > 0.05), 1), pct_gt10=round(100*np.mean(np.abs(r-1) > 0.10), 1)))
Z = pd.DataFrame(zrows); Z.to_csv(OUT / "zpve_consistency.csv", index=False); print("\n=== ZPVE consistency vs 3-model median ===\n" + Z.to_string(index=False))
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import run_thermo_protocol as P
def hybrid(model_E, model_vib):
    rows = []
    for i in B.index[B.evaluate]:
        a = J[model_E].get(i, {}); v = J[model_vib].get(i, {})
        if "E_model_kcal" not in a or "zpve" not in v: continue
        rec = dict(counts=a["counts"], E_model_kcal=a["E_model_kcal"], zpve=v["zpve"], H298_minus_E=v["H298_minus_E"])
        rows.append(dict(id=i, raw=P.dHf_from_record(rec, a["atom_E_kcal"])[0]))
    H = pd.DataFrame(rows).set_index("id").join(B[["dHf_exp_kcal", "n_H", "n_C", "n_N", "n_O", "n_F"]])
    N = H[["n_H", "n_C", "n_N", "n_O", "n_F"]].values.astype(float); y = (H.dHf_exp_kcal - H.raw).values
    rng = np.random.default_rng(2026); perm = rng.permutation(len(B[B.evaluate])); FOLD = pd.Series(np.empty(len(perm), int), index=B.index[B.evaluate])
    for k, te in enumerate(np.array_split(perm, 10)): FOLD.iloc[te] = k
    f = FOLD.loc[H.index].values; p = np.zeros(len(H))
    for k in range(10):
        te = f == k; p[te] = N[te] @ np.linalg.lstsq(N[~te], y[~te], rcond=None)[0]
    e = H.raw.values + p - H.dHf_exp_kcal.values
    return len(e), np.abs(e).mean(), np.median(np.abs(e))
hrows = []
for mE, mV in [("egret-1e", "egret-1e"), ("egret-1e", "uma-m"), ("egret-1e", "egret-1"), ("egret-1", "egret-1"), ("egret-1", "uma-m"), ("uma-m", "egret-1e"), ("mace-off23-medium", "uma-m")]:
    n, mae, med = hybrid(mE, mV); hrows.append(dict(electronic_energy=mE, vibrational_terms=mV, n=n, MAE=round(mae, 2), MedAE=round(med, 2)))
HY = pd.DataFrame(hrows); HY.to_csv(OUT / "hybrid_energy_vibration_check.csv", index=False); print("\n=== energy / vibration decomposition (CV atom equivalents) ===\n" + HY.to_string(index=False))
