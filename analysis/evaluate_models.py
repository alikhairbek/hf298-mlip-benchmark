#!/usr/bin/env python3
"""Common evaluation layer for every model run (merged CSVs from postprocess_runs.py).

Error definitions (kcal/mol, model - experiment):
  err_raw      atomisation route with the model's own atomic energies (no fitted parameter)
  err_cv       + element-wise atom-equivalent corrections (5 parameters, H C N O F), 10-fold cross-validated with the SAME folds for
               every model (folds assigned once from the benchmark ids), so paired comparisons are like-for-like
  err_loo      same corrections, leave-one-out (analytic hat-matrix formula) — robustness check on the fold choice
  err_sizefit  corrections fitted on the small stratum (<= 9 heavy atoms) only and applied to the large stratum (>= 10):
               the size-extrapolation test
Bootstrap: 10 000 resamples of molecules (seed 0) for 95 % percentile CIs of the MAE; paired bootstrap for differences.
"""
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).parent / "output"; OUT = ROOT / "eval"; OUT.mkdir(exist_ok=True)
MODELS = ["uma-m", "uma-s", "mace-off23-medium", "mace-off23-small", "aimnet2", "egret-1", "egret-1e", "ani-2x", "ani-1ccx", "xtb-gfn2"]
LABEL = {"uma-m": "UMA-m 1.1", "uma-s": "UMA-s 1.1", "mace-off23-medium": "MACE-OFF23 medium", "mace-off23-small": "MACE-OFF23 small",
         "aimnet2": "AIMNet2", "ani-2x": "ANI-2x", "ani-1ccx": "ANI-1ccx", "xtb-gfn2": "GFN2-xTB", "egret-1": "Egret-1", "egret-1e": "Egret-1e"}
ELS = ["n_H", "n_C", "n_N", "n_O", "n_F"]
NBOOT = 10000

def metrics(e):
    e = np.asarray(e, float); a = np.abs(e)
    return dict(n=len(e), MAE=a.mean(), RMSE=np.sqrt((e**2).mean()), MedAE=np.median(a), MaxAE=a.max(), bias=e.mean(),
                within1=100*np.mean(a <= 1), within2=100*np.mean(a <= 2), within5=100*np.mean(a <= 5))
def boot_mae(e, seed=0):
    rng = np.random.default_rng(seed); e = np.abs(np.asarray(e, float)); idx = rng.integers(0, len(e), (NBOOT, len(e)))
    return np.percentile(e[idx].mean(1), [2.5, 97.5])
def lstsq(N, y): return np.linalg.lstsq(N, y, rcond=None)[0]
def loo_residuals(N, y):
    """leave-one-out residuals of the no-intercept linear fit y ~ N c :  r_i / (1 - h_ii)"""
    c = lstsq(N, y); r = y - N @ c
    H = N @ np.linalg.pinv(N.T @ N) @ N.T; h = np.clip(np.diag(H), 0, 0.999999)
    return r / (1 - h)

B = pd.read_csv(str(Path(__file__).resolve().parent.parent / "benchmark_verified.csv"))
B = B[B.evaluate].copy()
fold_rng = np.random.default_rng(2026); perm = fold_rng.permutation(len(B)); FOLD = pd.Series(np.empty(len(B), int), index=B.id.values)
for k, te in enumerate(np.array_split(perm, 10)): FOLD.iloc[te] = k

evaluated, rows, atomeq = {}, [], []
def add(model, label, e, extra=None):
    e = np.asarray(e, float); e = e[np.isfinite(e)]
    if len(e) < 5: return
    lo, hi = boot_mae(e); d = dict(model=model, subset=label, **metrics(e), MAE_CI_lo=lo, MAE_CI_hi=hi); d.update(extra or {}); rows.append(d)

for m in MODELS:
    M = pd.read_csv(ROOT / f"{m}_merged.csv")
    M = M[M.evaluate].copy()
    n_failed = int(M.failed.sum()); failed_ids = ";".join(M[M.failed].id)
    M = M[~M.failed].copy()
    N = M[ELS].values.astype(float); y = (M.dHf_exp_kcal - M.dHf298_raw).values      # correction target
    fold = FOLD.loc[M.id].values
    pred_cv = np.zeros(len(M))
    for k in range(10):
        te = fold == k; pred_cv[te] = N[te] @ lstsq(N[~te], y[~te])
    M["err_raw"] = M.dHf298_raw - M.dHf_exp_kcal
    M["err_cv"] = M.dHf298_raw + pred_cv - M.dHf_exp_kcal
    M["pred_cv"] = M.dHf298_raw + pred_cv
    M["err_loo"] = -loo_residuals(N, y)                       # exp - (raw + corr_loo) = r  ->  err = -r
    small = (M.n_heavy <= 9).values
    c_small = lstsq(N[small], y[small]); M["err_sizefit"] = np.where(small, np.nan, M.dHf298_raw + N @ c_small - M.dHf_exp_kcal)
    c_all = lstsq(N, y)
    # per-fold stability of the coefficients
    cf = np.array([lstsq(N[fold != k], y[fold != k]) for k in range(10)])
    atomeq.append(dict(model=m, **{e[2:]: c for e, c in zip(ELS, c_all)}, **{f"{e[2:]}_cv_sd": s for e, s in zip(ELS, cf.std(0))},
                       **{f"{e[2:]}_smallfit": c for e, c in zip(ELS, c_small)}))
    evaluated[m] = M
    M.to_csv(OUT / f"{m}_evaluated.csv", index=False)
    t = dict(hours=M.seconds.sum()/3600, median_s=M.seconds.median(), n_failed=n_failed, failed_ids=failed_ids,
             n_imag=int((M.n_imag > 0).sum()), n_floor=int((M.n_floor_modes > 0).sum()))
    add(m, "raw atomisation (no fit)", M.err_raw, t)
    add(m, "ALL — atom equivalents, 10-fold CV", M.err_cv, t)
    add(m, "ALL — atom equivalents, LOO", M.err_loo)
    add(m, "small (≤9 heavy) — CV", M.err_cv[small]); add(m, "large (≥10 heavy) — CV", M.err_cv[~small])
    add(m, "large (≥10) — fitted on small only (size extrapolation)", M.err_sizefit)
    for tt in ["A+", "A", "B", "C"]: add(m, f"tier {tt} — CV", M.err_cv[M.tier == tt])
    add(m, "tiers A+ and A — CV", M.err_cv[M.tier.isin(["A+", "A"])])
    add(m, "tier A+ (ATcT) — raw atomisation", M.err_raw[M.tier == "A+"])
    add(m, "stereo assumed — CV", M.err_cv[M.stereo_assumed]); add(m, "stereo unambiguous — CV", M.err_cv[~M.stereo_assumed])
    add(m, "in SPICE v1 — CV", M.err_cv[M.in_spice_v1_inputs]); add(m, "not in SPICE v1 — CV", M.err_cv[~M.in_spice_v1_inputs])
    add(m, "in QM9 — CV", M.err_cv[M.in_qm9]); add(m, "not in QM9 — CV", M.err_cv[~M.in_qm9])
    hc = (M.n_N == 0) & (M.n_O == 0) & (M.n_F == 0)
    add(m, "hydrocarbons — CV", M.err_cv[hc]); add(m, "O only (no N, no F) — CV", M.err_cv[(M.n_O > 0) & (M.n_N == 0) & (M.n_F == 0)])
    add(m, "contains N — CV", M.err_cv[M.n_N > 0]); add(m, "contains F — CV", M.err_cv[M.n_F > 0])
    add(m, "aromatic (≥1 ring) — CV", M.err_cv[M.n_arom_rings > 0]); add(m, "≥2 aromatic rings — CV", M.err_cv[M.n_arom_rings >= 2])
    add(m, "non-aromatic rings only — CV", M.err_cv[(M.n_rings > 0) & (M.n_arom_rings == 0)]); add(m, "acyclic — CV", M.err_cv[M.n_rings == 0])
    add(m, "charge class: neutral — CV", M.err_cv[M.charge_class == "neutral"])
    add(m, "charge class: standard dipolar group — CV", M.err_cv[M.charge_class.str.startswith("standard")])
    add(m, "charge class: other charge-separated — CV", M.err_cv[M.charge_class.str.startswith("other")])
    add(m, "residual imaginary mode — CV", M.err_cv[M.n_imag > 0]); add(m, "no imaginary mode — CV", M.err_cv[M.n_imag == 0])
    for lo, hi in [(1, 3), (4, 6), (7, 9), (10, 12), (13, 15), (16, 40)]:
        add(m, f"heavy atoms {lo}–{hi} — CV", M.err_cv[(M.n_heavy >= lo) & (M.n_heavy <= hi)])
    sl = np.polyfit(M.n_heavy, M.err_cv, 1)[0]; rows[-1]["slope_err_vs_nheavy"] = sl
S = pd.DataFrame(rows); S.to_csv(OUT / "summary_long.csv", index=False)
pd.DataFrame(atomeq).to_csv(OUT / "atom_equivalents.csv", index=False)

# ------------------------------------------------------------------ headline table
H = S[S.subset == "ALL — atom equivalents, 10-fold CV"].set_index("model")
R = S[S.subset == "raw atomisation (no fit)"].set_index("model")
A = S[S.subset == "tiers A+ and A — CV"].set_index("model")
X = S[S.subset == "large (≥10) — fitted on small only (size extrapolation)"].set_index("model")
L = S[S.subset == "ALL — atom equivalents, LOO"].set_index("model")
head = pd.DataFrame({"model": [LABEL[m] for m in MODELS], "n": [H.loc[m, "n"] for m in MODELS], "failed": [H.loc[m, "n_failed"] for m in MODELS],
                     "raw MAE": [R.loc[m, "MAE"] for m in MODELS], "raw bias": [R.loc[m, "bias"] for m in MODELS],
                     "CV MAE": [H.loc[m, "MAE"] for m in MODELS], "CI lo": [H.loc[m, "MAE_CI_lo"] for m in MODELS], "CI hi": [H.loc[m, "MAE_CI_hi"] for m in MODELS],
                     "LOO MAE": [L.loc[m, "MAE"] for m in MODELS], "RMSE": [H.loc[m, "RMSE"] for m in MODELS], "MedAE": [H.loc[m, "MedAE"] for m in MODELS],
                     "MaxAE": [H.loc[m, "MaxAE"] for m in MODELS], "bias": [H.loc[m, "bias"] for m in MODELS],
                     "%≤1": [H.loc[m, "within1"] for m in MODELS], "%≤2": [H.loc[m, "within2"] for m in MODELS], "%≤5": [H.loc[m, "within5"] for m in MODELS],
                     "A+/A MAE": [A.loc[m, "MAE"] for m in MODELS], "A+/A MedAE": [A.loc[m, "MedAE"] for m in MODELS],
                     "size-extrap. MAE (large)": [X.loc[m, "MAE"] for m in MODELS], "size-extrap. MedAE": [X.loc[m, "MedAE"] for m in MODELS],
                     "n imag": [H.loc[m, "n_imag"] for m in MODELS], "GPU-h": [H.loc[m, "hours"] for m in MODELS], "median s/mol": [H.loc[m, "median_s"] for m in MODELS]})
head.to_csv(OUT / "headline.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.float_format", "{:.2f}".format)
print(head.to_string(index=False))

# ------------------------------------------------------------------ paired comparisons (common molecules of each pair) and error correlations
E = pd.concat([evaluated[m].set_index("id").err_cv.rename(m) for m in MODELS], axis=1)
common_all = E.dropna().index; common_HCNOF = E.drop(columns="ani-1ccx").dropna().index
print(f"\ncommon molecules: all {len(MODELS)} models = {len(common_all)}; {len(MODELS)-1} models without ANI-1ccx = {len(common_HCNOF)}")
rng = np.random.default_rng(0)
pairs = []
for i, a in enumerate(MODELS):
    for b in MODELS[i+1:]:
        d = E[[a, b]].dropna(); ea, eb = np.abs(d[a].values), np.abs(d[b].values); diff = ea - eb
        idx = rng.integers(0, len(d), (NBOOT, len(d))); bs = diff[idx].mean(1)
        pairs.append(dict(model_a=a, model_b=b, n=len(d), MAE_a=ea.mean(), MAE_b=eb.mean(), diff_MAE=diff.mean(), CI_lo=np.percentile(bs, 2.5), CI_hi=np.percentile(bs, 97.5),
                          p_two_sided=2*min((bs <= 0).mean(), (bs >= 0).mean()), pearson_r=np.corrcoef(d[a], d[b])[0, 1], spearman_r=d[a].corr(d[b], method="spearman"),
                          frac_same_sign=np.mean(np.sign(d[a]) == np.sign(d[b]))))
PAIRS = pd.DataFrame(pairs); PAIRS.to_csv(OUT / "paired_comparisons.csv", index=False)
C = E.loc[common_HCNOF].drop(columns="ani-1ccx").corr(); C.to_csv(OUT / "error_correlation_pearson.csv")
E.loc[common_all].corr().to_csv(OUT / "error_correlation_pearson_all8_HCNO.csv")
E.loc[common_HCNOF].drop(columns="ani-1ccx").corr(method="spearman").to_csv(OUT / "error_correlation_spearman.csv")
print("\nPearson correlation of CV errors (models without ANI-1ccx, common HCNOF set):"); print(C.round(2).to_string())
E.to_csv(OUT / "errors_cv_wide.csv")
print("\nMAE on the common HCNO set (all models, n = %d):" % len(common_all)); print(E.loc[common_all].abs().mean().round(2).to_string())

# ------------------------------------------------------------------ model consensus -> suspect experimental values
good = ["uma-m", "uma-s", "mace-off23-medium", "aimnet2", "ani-2x"]
G = E[good]; cons = G.median(1); mn = G.min(1); mx = G.max(1)
same_dir = ((mn > 3) | (mx < -3)) & G.notna().all(1)              # every one of the five models misses by > 3 kcal/mol in the same direction
sus = B.set_index("id").loc[same_dir[same_dir].index, ["name", "formula", "n_heavy", "dHf_exp_kcal", "unc_kcal", "tier", "value_source", "n_sources", "sources", "stereo_assumed", "charge_class"]].copy()
sus["consensus_err"] = cons[sus.index]; sus["min_err"] = mn[sus.index]; sus["max_err"] = mx[sus.index]
for m in good: sus[m] = G.loc[sus.index, m]
sus = sus.sort_values("consensus_err", key=np.abs, ascending=False)
sus.to_csv(OUT / "consensus_suspect_references.csv")
print(f"\nmolecules where all five best models miss by > 3 kcal/mol in the same direction: {len(sus)} "
      f"(tier A+: {(sus.tier=='A+').sum()}, A: {(sus.tier=='A').sum()}, B: {(sus.tier=='B').sum()}, C: {(sus.tier=='C').sum()})")
print(sus[["name", "formula", "dHf_exp_kcal", "tier", "n_sources", "consensus_err", "min_err", "max_err"]].head(30).to_string())
