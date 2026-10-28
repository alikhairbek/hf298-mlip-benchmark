#!/usr/bin/env python3
"""Cross-model figures (Wong colour-blind palette, no black plot elements, 300 dpi PNG + PDF)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
OUT = Path(__file__).parent / "output" / "eval"; FIG = Path(__file__).parent / "output" / "figures"; FIG.mkdir(exist_ok=True, parents=True)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 9, "axes.edgecolor": "#444444", "axes.labelcolor": "#222222", "xtick.color": "#333333",
                     "ytick.color": "#333333", "text.color": "#222222", "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "savefig.dpi": 300, "figure.dpi": 110})
INK = "#333333"; GRID = "#DDDDDD"
MODELS = ["uma-m", "uma-s", "aimnet2", "mace-off23-medium", "mace-off23-small", "egret-1", "egret-1e", "ani-1ccx", "ani-2x", "xtb-gfn2"]
LABEL = {"uma-m": "UMA-m", "uma-s": "UMA-s", "mace-off23-medium": "MACE-OFF23 medium", "mace-off23-small": "MACE-OFF23 small",
         "aimnet2": "AIMNet2", "ani-2x": "ANI-2x", "ani-1ccx": "ANI-1ccx", "xtb-gfn2": "GFN2-xTB", "egret-1": "Egret-1", "egret-1e": "Egret-1e"}
COL = {"uma-m": "#0072B2", "uma-s": "#0072B2", "mace-off23-medium": "#D55E00", "mace-off23-small": "#D55E00", "aimnet2": "#009E73",
       "ani-2x": "#E69F00", "ani-1ccx": "#CC79A7", "xtb-gfn2": "#7F7F7F", "egret-1": "#56B4E9", "egret-1e": "#56B4E9"}
OPEN = {"uma-s", "mace-off23-small", "egret-1e"}                       # open markers for the small sibling of a model pair
MK = {m: ("o" if m not in OPEN else "o") for m in MODELS}
def style(m, **kw):
    d = dict(color=COL[m], marker="o", ms=5.5, mec=COL[m], mfc=(COL[m] if m not in OPEN else "white"), mew=1.3); d.update(kw); return d
def save(fig, name):
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight"); fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight"); plt.close(fig); print("saved", name)

S = pd.read_csv(OUT / "summary_long.csv"); H = pd.read_csv(OUT / "headline.csv")
B = pd.read_csv(str(Path(__file__).resolve().parent.parent / "benchmark_verified.csv")).set_index("id")
E = pd.read_csv(OUT / "errors_cv_wide.csv", index_col=0)
def sub(model, subset): r = S[(S.model == model) & (S.subset == subset)]; return r.iloc[0] if len(r) else None

# ------------------------------------------------------------------ Figure 1: headline + tiers
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.9), gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axes[0]; ys = np.arange(len(MODELS))[::-1]
for y, m in zip(ys, MODELS):
    r = sub(m, "ALL — atom equivalents, 10-fold CV")
    ax.plot([r.MAE_CI_lo, r.MAE_CI_hi], [y, y], color=COL[m], lw=1.6, solid_capstyle="round")
    ax.plot(r.MAE, y, **style(m, ls="none"))
    ax.plot(r.MedAE, y, marker="|", ms=9, mew=1.6, color=COL[m], ls="none")
    ax.text(r.MAE_CI_hi + 0.15, y, f"{r.MAE:.2f}", va="center", ha="left", fontsize=8, color=INK)
ax.set_yticks(ys); ax.set_yticklabels([LABEL[m] + (" (HCNO)" if m == "ani-1ccx" else "") for m in MODELS])
ax.set_xlabel("MAE, kcal mol$^{-1}$ (CV atom equivalents)"); ax.set_xlim(0, 9); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], marker="o", color=INK, ls="none", ms=5, label="MAE (95 % bootstrap CI)"), Line2D([], [], marker="|", color=INK, ls="none", ms=9, mew=1.6, label="median |error|")],
          loc="upper right", fontsize=7.5)
ax.set_title("a  All 1 579 molecules", loc="left", fontsize=9.5, fontweight="bold")
ax = axes[1]; tiers = ["tier A+ — CV", "tier A — CV", "tier B — CV", "tier C — CV"]; xt = np.arange(4)
SIBLING = {"uma-s", "mace-off23-small"}                   # size siblings left out of the multi-model panels
PLOT = [m for m in MODELS if m not in SIBLING]
for i, m in enumerate(PLOT):
    v = [sub(m, t).MAE for t in tiers]; lo = [sub(m, t).MAE_CI_lo for t in tiers]; hi = [sub(m, t).MAE_CI_hi for t in tiers]
    off = (i - (len(PLOT) - 1) / 2) * 0.08
    ax.errorbar(xt + off, v, yerr=[np.array(v) - lo, np.array(hi) - v], color=COL[m], lw=1.2, capsize=0, ls="-", **{k: s for k, s in style(m).items() if k not in ("color", "ls")}, label=LABEL[m])
ax.set_xticks(xt); ax.set_xticklabels(["A+ (ATcT)\nn = 177", "A\nn = 434", "B\nn = 936", "C\nn = 32"]); ax.set_xlabel("reference-value tier")
ax.set_ylabel("MAE, kcal mol$^{-1}$"); ax.set_ylim(0, 13.5); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.legend(fontsize=6.5, ncol=2, loc="upper left", handlelength=1.2, columnspacing=0.8)
ax.set_title("b  By reference quality", loc="left", fontsize=9.5, fontweight="bold")
save(fig, "fig1_headline_and_tiers")

# ------------------------------------------------------------------ Figure 2: error-correlation structure
REP = ["uma-m", "aimnet2", "mace-off23-medium", "egret-1", "egret-1e", "ani-1ccx", "ani-2x", "xtb-gfn2"]; SHORT = {**LABEL, "mace-off23-medium": "MACE-OFF23 med."}
fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), gridspec_kw=dict(width_ratios=[1, 1, 1.1], wspace=0.08), constrained_layout=True)
for k, (ax, (title, csv)) in enumerate(zip(axes[:2], [("a  All molecules (n = 1 421)", "error_correlation_all.csv"), ("b  Tier A+ only (n = 148)", "error_correlation_tier_Aplus.csv")])):
    C = pd.read_csv(OUT / csv, index_col=0).loc[REP, REP].values
    im = ax.imshow(C, cmap="RdBu", vmin=-1, vmax=1)
    for i in range(len(REP)):
        for j in range(len(REP)):
            ax.text(j, i, f"{C[i, j]:.2f}", ha="center", va="center", fontsize=5.8, color=("white" if abs(C[i, j]) > 0.65 else INK))
    ax.set_xticks(range(len(REP))); ax.set_yticks(range(len(REP))); ax.set_xticklabels([SHORT[m] for m in REP], rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels([SHORT[m] for m in REP] if k == 0 else [], fontsize=7)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold"); [sp.set_visible(False) for sp in ax.spines.values()]; ax.tick_params(length=0)
cb = fig.colorbar(im, ax=axes[:2], shrink=0.75, pad=0.02, aspect=28); cb.set_label("Pearson r of CV errors", fontsize=7.5); cb.ax.tick_params(labelsize=7); cb.outline.set_edgecolor("#888888")
ax = axes[2]
d = E[["uma-m", "mace-off23-medium"]].dropna(); t = B.loc[d.index, "tier"]
for lab, mask, c, z in [("tier B / C", t.isin(["B", "C"]), "#E69F00", 1), ("tier A+ / A", t.isin(["A+", "A"]), "#0072B2", 2)]:
    ax.plot(d["uma-m"][mask], d["mace-off23-medium"][mask], ls="none", marker="o", ms=2.6, mfc=c, mec="none", alpha=0.55, label=f"{lab} (n = {int(mask.sum())})", zorder=z)
ax.plot([-50, 50], [-50, 50], color="#999999", lw=0.8, ls="--", zorder=0); ax.set_xlim(-25, 25); ax.set_ylim(-25, 25); ax.set_aspect("equal")
ax.set_xlabel("UMA-m error, kcal mol$^{-1}$"); ax.set_ylabel("MACE-OFF23 medium error, kcal mol$^{-1}$"); ax.legend(fontsize=7, loc="upper left", handletextpad=0.2, markerscale=1.8)
ax.set_title("c  Shared errors follow the reference tier", loc="left", fontsize=9, fontweight="bold"); ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
save(fig, "fig2_error_correlation")

# ------------------------------------------------------------------ Figure 3: size and chemistry
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.3), gridspec_kw=dict(width_ratios=[1, 1.2]), constrained_layout=True)
ax = axes[0]; bins = ["heavy atoms 1–3 — CV", "heavy atoms 4–6 — CV", "heavy atoms 7–9 — CV", "heavy atoms 10–12 — CV", "heavy atoms 13–15 — CV", "heavy atoms 16–40 — CV"]
xb = np.arange(len(bins))
ax.axvspan(2.5, 5.5, color="#F0E442", alpha=0.18, lw=0)
for m in PLOT:
    v = [sub(m, b).MAE for b in bins]
    ax.plot(xb, v, **style(m, lw=1.3, ls="-"), label=LABEL[m])
n_b = [int(sub("uma-m", b).n) for b in bins]
ax.set_xticks(xb); ax.set_xticklabels([f"{l}\nn = {n}" for l, n in zip(["1–3", "4–6", "7–9", "10–12", "13–15", "16–33"], n_b)], fontsize=7.5)
ax.set_yscale("log"); ax.set_ylim(1, 60); ax.set_yticks([1, 2, 5, 10, 20, 50]); ax.set_yticklabels(["1", "2", "5", "10", "20", "50"])
ax.set_xlabel("heavy atoms"); ax.set_ylabel("MAE, kcal mol$^{-1}$ (log scale)"); ax.grid(axis="y", color=GRID, lw=0.6, which="major"); ax.set_axisbelow(True)
ax.text(4, 1.12, "size-extrapolation stratum\n(≥ 10 heavy atoms)", ha="center", va="bottom", fontsize=7, color=INK)
ax.text(5, 42, "perfluoro- and\npolynitro only", ha="center", va="top", fontsize=6.5, color=INK, style="italic")
ax.legend(fontsize=7, loc="upper left", ncol=2, handlelength=1.4, columnspacing=0.8); ax.set_title("a  Molecule size", loc="left", fontsize=9.5, fontweight="bold")
ax = axes[1]; classes = [("hydrocarbons — CV", "hydro-\ncarbons"), ("O only (no N, no F) — CV", "C,H,O\nonly"), ("contains N — CV", "contains\nN"), ("contains F — CV", "contains\nF"),
                         ("charge class: standard dipolar group — CV", "nitro, N-oxide,\nazide, isocyanide"), ("≥2 aromatic rings — CV", "≥ 2 aromatic\nrings")]
xc = np.arange(len(classes))
for i, m in enumerate(PLOT):
    v = [(sub(m, c).MAE if sub(m, c) is not None else np.nan) for c, _ in classes]
    off = (i - (len(PLOT) - 1) / 2) * 0.09
    ax.plot(xc + off, v, **style(m, ls="none"), label=LABEL[m])
n_c = [int(sub("uma-m", c).n) for c, _ in classes]
ax.set_xticks(xc); ax.set_xticklabels([f"{l}\nn = {n}" for (_, l), n in zip(classes, n_c)], fontsize=7); ax.set_ylabel("MAE, kcal mol$^{-1}$"); ax.set_ylim(0, 24)
ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_title("b  Chemistry", loc="left", fontsize=9.5, fontweight="bold")
ax.legend(fontsize=7, loc="upper left", ncol=3, handlelength=1.0, columnspacing=0.7, handletextpad=0.3)
save(fig, "fig3_size_and_chemistry")

# ------------------------------------------------------------------ Figure 4: UMA-m residuals and the contradicted reference values
U = pd.read_csv(OUT / "uma-m_evaluated.csv").set_index("id")
SUS = pd.read_csv(OUT / "family_consensus_suspect_references.csv", index_col=0)
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), gridspec_kw=dict(width_ratios=[1, 1.15]), constrained_layout=True)
ax = axes[0]
for lab, mask, c, z in [("tier B / C", U.tier.isin(["B", "C"]), "#E69F00", 1), ("tier A+ / A", U.tier.isin(["A+", "A"]), "#0072B2", 2)]:
    ax.plot(U.dHf_exp_kcal[mask], U.err_cv[mask], ls="none", marker="o", ms=2.6, mfc=c, mec="none", alpha=0.6, label=f"{lab} (n = {int(mask.sum())})", zorder=z)
sids = SUS.index.intersection(U.index)
ax.plot(U.dHf_exp_kcal[sids], U.err_cv[sids], ls="none", marker="o", ms=5, mfc="none", mec="#D55E00", mew=0.9, label=f"contradicted by all model families (n = {len(sids)})", zorder=3)
ax.axhline(0, color="#999999", lw=0.8, zorder=0); ax.set_xlim(-1200, 200); ax.set_ylim(-52, 52)
ax.set_xlabel("experimental ΔH$_f$°(298 K), kcal mol$^{-1}$"); ax.set_ylabel("UMA-m − experiment, kcal mol$^{-1}$"); ax.legend(fontsize=7, loc="lower left", handletextpad=0.2, markerscale=1.5)
ax.set_title("a  UMA-m residuals (1 579 molecules)", loc="left", fontsize=9.5, fontweight="bold"); ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
ax = axes[1]; top = SUS.head(20); yy = np.arange(len(top))[::-1]
ax.axvline(0, color="#999999", lw=0.8)
for y, (idx, r) in zip(yy, top.iterrows()):
    ax.plot(r.wB97M_family_median, y, marker="o", ms=5, color="#0072B2", ls="none", zorder=3)
    if np.isfinite(r["ani-1ccx"]): ax.plot(r["ani-1ccx"], y, marker="s", ms=4.2, color="#CC79A7", ls="none", zorder=3)
    if np.isfinite(r["ani-2x"]): ax.plot(r["ani-2x"], y, marker="^", ms=4.5, color="#E69F00", ls="none", zorder=3)
def short(n, k=28): return n if len(n) <= k else n[:k - 1] + "…"
ax.set_yticks(yy); ax.set_yticklabels([f"{short(r['name'])} ({r.tier}, {int(r.n_sources)} src)" for _, r in top.iterrows()], fontsize=6.5)
ax.set_xlabel("model − compilation value, kcal mol$^{-1}$"); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_xlim(-65, 45)
ax.legend(handles=[Line2D([], [], marker="o", color="#0072B2", ls="none", ms=5, label="ωB97M family median\n(UMA-m, AIMNet2, MACE-OFF23)"),
                   Line2D([], [], marker="s", color="#CC79A7", ls="none", ms=4.2, label="ANI-1ccx (CCSD(T)*/CBS)"), Line2D([], [], marker="^", color="#E69F00", ls="none", ms=4.5, label="ANI-2x (ωB97X/6-31G*)")],
          fontsize=6.5, loc="center left", bbox_to_anchor=(0.0, 0.55), labelspacing=0.6)
ax.set_title("b  20 largest contradicted reference values", loc="left", fontsize=9.5, fontweight="bold")
save(fig, "fig4_residuals_and_contradicted_references")

# ------------------------------------------------------------------ Figure 5: leakage tests
LK = pd.read_csv(OUT / "leakage_tests.csv")
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.3), gridspec_kw=dict(width_ratios=[1, 1.45]), constrained_layout=True)
ax = axes[0]; tests = ["|err mace-off23-medium| − |err uma-m|", "|err mace-off23-small| − |err uma-m|", "|err egret-1| − |err uma-m|", "|err mace-off23-medium| − |err aimnet2|"]
labs = ["MACE-OFF23 medium\nvs UMA-m", "MACE-OFF23 small\nvs UMA-m", "Egret-1\nvs UMA-m", "MACE-OFF23 medium\nvs AIMNet2"]
for i, (t, l) in enumerate(zip(tests, labs)):
    for j, (subn, c, mk) in enumerate([("in SPICE v1", "#D55E00", "o"), ("not in SPICE v1", "#0072B2", "s")]):
        r = LK[(LK.test == t) & (LK.subset == subn)].iloc[0]; y = i + (0.15 if j == 0 else -0.15)
        ax.plot([r.CI_lo, r.CI_hi], [y, y], color=c, lw=1.6); ax.plot(r.mean_diff, y, marker=mk, color=c, ms=5, ls="none", label=(f"{subn} (n = {int(r.n)})" if i == 0 else None))
ax.axvline(0, color="#999999", lw=0.8); ax.set_yticks(range(len(tests))); ax.set_yticklabels(labs, fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("paired difference of |error|, kcal mol$^{-1}$\n(> 0: first-named model worse)"); ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, columnspacing=1.0); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_title("a  Paired test on SPICE membership", loc="left", fontsize=9.5, fontweight="bold")
ax = axes[1]; mods = ["mace-off23-medium", "mace-off23-small", "egret-1", "uma-m", "aimnet2", "ani-1ccx", "ani-2x"]
for i, m in enumerate(mods):
    for j, (subn, mk, fc) in enumerate([("in SPICE v1", "o", COL[m]), ("not in SPICE v1", "s", "white")]):
        r = LK[(LK.test == f"MAE {m} (≤9 heavy, tiers A+/A)") & (LK.subset == subn)].iloc[0]; x = i + (-0.17 if j == 0 else 0.17)
        ax.plot([x, x], [r.CI_lo, r.CI_hi], color=COL[m], lw=1.5); ax.plot(x, r.mean_diff, marker=mk, ms=5.5, mfc=fc, mec=COL[m], mew=1.3, ls="none")
ax.set_xticks(range(len(mods))); ax.set_xticklabels([LABEL[m].replace(" ", "\n") for m in mods], fontsize=6.5); ax.set_xlim(-0.6, len(mods) - 0.4); ax.set_ylabel("MAE, kcal mol$^{-1}$"); ax.set_ylim(0, 10.5)
ax.legend(handles=[Line2D([], [], marker="o", color=INK, ls="none", ms=5.5, label="in SPICE v1 (n = 126)"), Line2D([], [], marker="s", mfc="white", mec=INK, color=INK, ls="none", ms=5.5, label="not in SPICE v1 (n = 376)")], fontsize=7, loc="upper left")
ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_title("b  Size- and tier-matched (≤ 9 heavy atoms, tiers A+/A)", loc="left", fontsize=9.5, fontweight="bold")
save(fig, "fig5_leakage_tests")
