#!/usr/bin/env python3
"""Paper tables from the evaluation outputs (run after evaluate_models.py).

  Table 1  output/eval/table1_headline.csv                      headline accuracy, CV/LOO, size extrapolation, cost
  Table 2  output/eval/table2_subsets_MAE_only.csv              MAE per subset (one column per model)
           output/eval/table2_subsets_MAE_MedAE.csv             MAE [95% CI] and MedAE per subset
  Table 3  output/eval/table3_error_by_reference_provenance.csv  error vs. provenance of the reference value (ATcT / multi-source / single source)

Provenance groups are disjoint: a molecule whose value was replaced by ATcT (tier A+) is counted only in the ATcT row, so the
"single source" rows contain only molecules whose reference value really rests on that one compilation column.  A model's entry is
left blank when it does not cover the whole group and covers fewer than five of its molecules (ANI-1ccx, HCNO only).
Bootstrap 95 % CIs of the provenance MAEs: 10 000 resamples, numpy default_rng(0).
"""
import numpy as np, pandas as pd
from pathlib import Path

OUT = Path(__file__).parent / "output" / "eval"
MODELS = ["uma-m", "uma-s", "aimnet2", "mace-off23-medium", "mace-off23-small", "egret-1", "egret-1e", "ani-1ccx", "ani-2x", "xtb-gfn2"]
NBOOT = 10_000
CV = "ALL — atom equivalents, 10-fold CV"

S = pd.read_csv(OUT / "summary_long.csv")
head = pd.read_csv(OUT / "headline.csv")
head.to_csv(OUT / "table1_headline.csv", index=False)

# ---------------------------------------------------------------- Table 2: subsets
subsets = list(dict.fromkeys(S.subset))                       # order of first appearance
mae = S.pivot(index="subset", columns="model", values="MAE").loc[subsets, MODELS].round(2)
mae.to_csv(OUT / "table2_subsets_MAE_only.csv")
n = S[S.model == "uma-m"].set_index("subset").loc[subsets, "n"].astype(int)
cols = {"n (UMA-m)": n}
for m in MODELS:
    T = S[S.model == m].set_index("subset").reindex(subsets)
    cols[f"{m} MAE [95% CI]"] = [("" if np.isnan(r.MAE) else f"{r.MAE:.2f} [{r.MAE_CI_lo:.2f}, {r.MAE_CI_hi:.2f}]") for r in T.itertuples()]
    cols[f"{m} MedAE"] = T.MedAE.round(2).values
pd.DataFrame(cols, index=subsets).to_csv(OUT / "table2_subsets_MAE_MedAE.csv")

# ---------------------------------------------------------------- Table 3: provenance of the reference value
ev = {m: pd.read_csv(OUT / f"{m}_evaluated.csv").set_index("id") for m in MODELS}
B = ev["uma-m"]
groups = [("ATcT v1.222 (tier A+)", B.tier == "A+"),
          ("multi-source, tier A", (B.n_sources >= 2) & (B.tier == "A")),
          ("multi-source, tier B/C", (B.n_sources >= 2) & B.tier.isin(["B", "C"]))]
single = B[(B.n_sources == 1) & (B.tier != "A+")]
for src in ["Pedley", "Cox-1", "NIST-1", "Cioslowski", "Yaws", "Add Lit 1", "Sandler", "Winget"]:
    groups.append((f"single source: {src}", (B.n_sources == 1) & (B.tier != "A+") & (B.sources == src)))
rng = np.random.default_rng(0)
rows = []
for label, mask in groups:
    ids = B.index[mask]
    row = {"provenance": label, "n": int(len(ids))}
    for m in MODELS:
        e = ev[m].reindex(ids).err_cv.dropna().abs().values
        if len(e) == 0 or (len(e) < len(ids) and len(e) < 5):
            row.update({f"{m} MAE": np.nan, f"{m} CI_lo": np.nan, f"{m} CI_hi": np.nan, f"{m} MedAE": np.nan}); continue
        bs = e[rng.integers(0, len(e), (NBOOT, len(e)))].mean(1)
        row.update({f"{m} MAE": round(e.mean(), 2), f"{m} CI_lo": round(np.percentile(bs, 2.5), 2), f"{m} CI_hi": round(np.percentile(bs, 97.5), 2), f"{m} MedAE": round(np.median(e), 2)})
    rows.append(row)
pd.DataFrame(rows).to_csv(OUT / "table3_error_by_reference_provenance.csv", index=False)
print("tables written:", ", ".join(f"table{i}" for i in (1, 2, 3)))
