#!/usr/bin/env python3
"""Unified post-processing of all model runs.

For every record the ZPVE and H(298.15 K) - E are recomputed from the stored harmonic frequencies with ONE rule
(complex magnitude of residual imaginary modes, 50 cm^-1 floor for all modes below it), so that runs made with
the first notebook version (no floor -> NaN enthalpies for purely imaginary modes) are treated exactly like the
later ones.  For runs that already used the rule the recomputation reproduces the stored values (checked below).

Rigid-rotor/harmonic-oscillator ideal gas, 298.15 K:
    H - E = ZPVE + sum_i h nu_i / (exp(h nu_i / kT) - 1) + 3/2 kT (trans) + (3/2 | 1) kT (rot, nonlinear | linear) + kT (pV)
(identical to ase.thermochemistry.IdealGasThermo.get_enthalpy with potentialenergy = 0).
"""
import json, sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import run_thermo_protocol as P

KT = 8.617333262e-5 * 298.15          # eV
CM = 1.0 / 8065.54                    # eV per cm^-1
FLOOR_CM = 50.0

def enthalpy_from_freqs(freqs_cm, geometry):
    f = np.asarray(freqs_cm, float)
    n_floor = int((f < FLOOR_CM).sum()); f = np.where(f < FLOOR_CM, FLOOR_CM, f)
    e = f * CM
    zpve = 0.5 * e.sum()
    vib = (e / (np.exp(e / KT) - 1.0)).sum()
    rot = {"monatomic": 0.0, "linear": 1.0, "nonlinear": 1.5}[geometry] * KT
    H = zpve + vib + 1.5 * KT + rot + KT
    return zpve * P.EV, H * P.EV, n_floor

MODELS = ["aimnet2", "ani-1ccx", "ani-2x", "mace-off23-medium", "mace-off23-small", "uma-m", "uma-s", "xtb-gfn2", "egret-1", "egret-1e"]
OUT = Path(__file__).parent / "output"; OUT.mkdir(exist_ok=True)
B = pd.read_csv(str(Path(__file__).resolve().parent.parent / "benchmark_verified.csv"))
KEEP = ["E_model_kcal", "zpve", "thermal", "H298_minus_E", "n_imag", "n_floor_modes", "lowest_cm", "geometry", "dHf298_raw",
        "seconds", "converged", "fmax_final", "reopt_tries", "conf_energy_spread", "mmff_rank_of_best", "n_conf_opt", "n_conf_embedded"]
report = []
for m in MODELS:
    recs = [json.loads(l) for l in open(Path(__file__).resolve().parent.parent / "raw_results" / f"{m}.jsonl")]
    rows, dz, dh, dd = [], [], [], []
    for r in recs:
        if "freqs_cm" not in r or r.get("E_model_kcal") is None:          # failed / unsupported molecule
            rows.append(dict(id=r["id"], supported=False)); continue
        zp, H, nfl = enthalpy_from_freqs(r["freqs_cm"], r["geometry"])
        old_H = r.get("H298_minus_E"); old_z = r.get("zpve"); old_d = r.get("dHf298_raw")
        r2 = dict(r); r2["zpve"] = zp; r2["H298_minus_E"] = H; r2["thermal"] = H - zp; r2["n_floor_modes"] = max(nfl, int(r.get("n_floor_modes") or 0))
        dHf, D0, dHf0 = P.dHf_from_record(r2, r["atom_E_kcal"])
        if old_H is not None and np.isfinite(old_H):                        # consistency check against the stored values
            dz.append(zp - old_z); dh.append(H - old_H); dd.append(dHf - old_d)
        row = {k: r2.get(k) for k in KEEP}; row.update(id=r["id"], supported=True, dHf298_raw=dHf, D0=D0, dHf0=dHf0,
                                                       recomputed=not (old_H is not None and np.isfinite(old_H)))
        rows.append(row)
    R = pd.DataFrame(rows)
    M = B.merge(R, on="id", how="left")
    M["supported"] = M.supported.fillna(False).astype(bool)
    M["evaluate"] = M.evaluate & M.supported & M.dHf298_raw.notna()
    # catastrophic model failure: optimiser did not converge and the final |force| is still > 0.5 eV/Å (structure blew up);
    # such molecules are counted and listed per model and excluded from the error statistics
    M["failed"] = M.supported & (~M.converged.fillna(True).astype(bool)) & (M.fmax_final.fillna(0) > 0.5)
    M.to_csv(OUT / f"{m}_merged.csv", index=False)
    rep = dict(model=m, n_records=len(recs), n_supported=int(R.supported.sum()), n_evaluated=int(M.evaluate.sum()), n_failed=int((M.failed & M.evaluate).sum()),
               failed_ids=";".join(M[M.failed & M.evaluate].id),
               n_recomputed_from_NaN=int(R.get("recomputed", pd.Series(dtype=bool)).fillna(False).sum()),
               max_abs_dZPVE_vs_stored=float(np.max(np.abs(dz))) if dz else np.nan,
               max_abs_dH298_vs_stored=float(np.max(np.abs(dh))) if dh else np.nan,
               max_abs_ddHf_vs_stored=float(np.max(np.abs(dd))) if dd else np.nan,
               n_changed_gt_0p01=int((np.abs(dd) > 0.01).sum()) if dd else 0,
               n_imag_gt0=int((R.n_imag > 0).sum()), n_floor_gt0=int((R.n_floor_modes > 0).sum()),
               hours=float(R.seconds.sum() / 3600), median_s=float(R.seconds.median()))
    report.append(rep); print(rep)
pd.DataFrame(report).to_csv(OUT / "postprocess_report.csv", index=False)
