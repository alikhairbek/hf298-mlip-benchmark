#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_thermo_protocol.py — physically complete ΔHf°(298 K) protocol for pretrained ML potentials (production version).

Usage (Kaggle / Colab GPU, or CLI):
    python run_thermo_protocol.py --model mace-off23-medium --input benchmark_verified.csv --out results/
    python run_thermo_protocol.py --model ani-2x          --input benchmark_verified.csv --out results/
    python run_thermo_protocol.py --model aimnet2         --input benchmark_verified.csv --out results/
    python run_thermo_protocol.py --model uma-s           --input benchmark_verified.csv --out results/   (needs `huggingface-cli login`)
    python run_thermo_protocol.py --model mace:/path/or/url/to/model.model --input ...             (any MACE checkpoint, e.g. MACE-OFF24)
    python run_thermo_protocol.py --model egret-1e        --input benchmark_verified.csv --out results/   (Egret-1 / -1e / -1t / -1s / -1m, auto-downloaded)
    python run_thermo_protocol.py --model xtb-gfn2        --input benchmark_verified.csv --out results/   (CPU only, no GPU needed)

Input CSV columns (required): id, smiles, dHf_exp_kcal   (optional: name, unc_kcal, source, stratum)
Output: <out>/<model>.jsonl (one JSON record per molecule, appended; restart-safe), <out>/<model>_summary.csv

Per molecule:
  1. RDKit ETKDGv3 conformers (n = 10 + 4·n_rot, ≤ 50), MMFF94 pre-optimisation, keep the K lowest (default 5)
  2. Optimise every kept conformer with the model (ASE LBFGS, fmax 0.01 eV/Å); keep the lowest-energy minimum
  3. Harmonic frequencies by central finite differences of model forces (δ = 0.01 Å); imaginary mode → perturb & re-optimise (≤2×)
  4. RRHO ideal-gas thermochemistry at 298.15 K / 1 bar (ZPVE, H298 − H0)
  5. Atomisation route: D0 = Σ n_i E_model(atom_i) − E_model(M) − ZPVE;
     ΔfH°(0 K) = Σ n_i ΔfH°(atom_i, 0 K) − D0;  ΔfH°(298) = ΔfH°(0 K) + [H298−H0](M) − Σ n_i [H298−H0](element_i, std. state)
  The record stores everything needed to re-derive ΔHf° with any atomic reference scheme (E_model, ZPVE, thermal, counts).
  Element-wise atom-equivalent corrections are fitted afterwards (cross-validated) by thermo_analysis.py — never here.
"""
import json, time, os, argparse, tempfile, shutil
import numpy as np, pandas as pd
from pathlib import Path

EV = 23.060548; KJ = 1/4.184
DHF0_ATOM = {"H": 216.034*KJ, "C": 711.38*KJ, "N": 470.57*KJ, "O": 246.844*KJ, "F": 77.27*KJ,     # ATcT, 0 K, kcal/mol
             "S": 274.73*KJ, "Cl": 119.621*KJ, "P": 315.6*KJ, "Br": 117.92*KJ}
H298_ELEM = {"H": 8.468/2*KJ, "C": 1.050*KJ, "N": 8.670/2*KJ, "O": 8.680/2*KJ, "F": 8.825/2*KJ,     # [H298−H0] per atom, std. state
             "S": 4.412*KJ, "Cl": 9.181/2*KJ, "P": 5.36*KJ, "Br": 24.52/2*KJ}
ELEMENTS = ["H", "C", "N", "O", "F", "S", "Cl", "P", "Br"]

# ------------------------------------------------------------------ models
def _torchscript_safe_mode():
    """TorchScript archives (the Egret-1 family, compiled MACE checkpoints) crash on CUDA with
    'RuntimeError: Global alloc not supported yet' (pytorch issues #69078, #144965): the profiling
    executor specialises the graph on the first two calls — here single atoms with *empty* edge
    tensors — and the fused kernel then fails on every later call.  Running the un-optimised graph
    (what torch.jit.optimized_execution(False) does) avoids the fuser altogether; the numbers are
    identical, only the kernel fusion speed-up is lost.  Harmless for pickled (non-scripted) models."""
    import torch
    for fn, arg in [(getattr(torch._C, "_set_graph_executor_optimize", None), False),
                    (getattr(torch._C, "_jit_set_texpr_fuser_enabled", None), False),
                    (getattr(torch._C, "_jit_override_can_fuse_on_gpu", None), False),
                    (getattr(torch._C, "_jit_set_nvfuser_enabled", None), False)]:
        try:
            if fn is not None: fn(arg)
        except Exception:
            pass

def make_calc(model, device):
    if model.startswith("mace-off23-"):
        from mace.calculators import mace_off
        return mace_off(model=model.split("-")[-1], device=device, default_dtype="float64")
    if model.startswith("mace:"):                                   # any MACE checkpoint (file path or URL), e.g. MACE-OFF24
        from mace.calculators import MACECalculator
        _torchscript_safe_mode()
        src = model[5:]
        if src.startswith("http"):
            import urllib.request
            path = Path.home()/".cache"/"mace_ckpt"/src.split("/")[-1]; path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists(): urllib.request.urlretrieve(src, path)
            src = str(path)
        return MACECalculator(model_paths=src, device=device, default_dtype="float64")
    if model == "ani-1ccx":
        import torchani, torch; return torchani.models.ANI1ccx(periodic_table_index=True).to(torch.device(device)).ase()
    if model == "ani-2x":
        import torchani, torch; return torchani.models.ANI2x(periodic_table_index=True).to(torch.device(device)).ase()
    if model == "aimnet2":                                          # pip install aimnet
        from aimnet.calculators import AIMNet2ASE
        return AIMNet2ASE("aimnet2")
    if model.startswith("uma-"):                                    # pip install fairchem-core ; huggingface-cli login (gated)
        from fairchem.core import pretrained_mlip, FAIRChemCalculator
        name = {"uma-s": "uma-s-1p1", "uma-m": "uma-m-1p1"}.get(model, model[4:])
        pred = pretrained_mlip.get_predict_unit(name, device=device)
        return FAIRChemCalculator(pred, task_name="omol")
    if model.startswith("egret-1"):                                 # Rowan Egret-1 family: MACE-architecture TorchScript models, MIT licence, no token needed
        import urllib.request
        from mace.calculators import MACECalculator
        _torchscript_safe_mode()                                    # required on GPU, see _torchscript_safe_mode()
        fname = {"egret-1": "EGRET_1.model", "egret-1e": "EGRET_1E.model", "egret-1t": "EGRET_1T.model",
                 "egret-1s": "EGRET_1S.model", "egret-1m": "EGRET_1M.model"}[model]
        path = Path.home()/".cache"/"egret"/fname; path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            urllib.request.urlretrieve(f"https://raw.githubusercontent.com/rowansci/egret-public/master/compiled_models/{fname}", path)
        return MACECalculator(model_paths=str(path), device=device, default_dtype="float64")
    if model == "xtb-gfn2":                                         # pip install tblite
        from tblite.ase import TBLite
        return TBLite(method="GFN2-xTB", verbosity=0)
    raise ValueError(f"unknown model {model}")

SUPPORTS = {"ani-1ccx": {"H","C","N","O"}, "ani-2x": {"H","C","N","O","F","S","Cl"}}

# ------------------------------------------------------------------ helpers
def conformers(smiles, k, seed=42):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    nrot = AllChem.CalcNumRotatableBonds(mol)
    n = int(min(50, 10 + 4*nrot))
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed; ps.pruneRmsThresh = 0.3; ps.enforceChirality = True
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n, params=ps))
    if not cids:                                                   # fallbacks: random coordinates, more attempts, relaxed chirality enforcement
        ps.useRandomCoords = True; ps.maxIterations = 2000; cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n, params=ps))
    if not cids:
        ps.enforceChirality = False; ps.ignoreSmoothingFailures = True; cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n, params=ps))
    if not cids:
        raise RuntimeError("RDKit could not embed any conformer")
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    ranked = sorted([(e, cid) for (conv, e), cid in zip(res, cids)])
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    return [(syms, np.array(mol.GetConformer(cid).GetPositions()), e) for e, cid in ranked[:k]], nrot, len(cids)

def new_atoms(symbols, positions):
    """All structures are neutral closed-shell species; charge/spin are set explicitly because the
    UMA (fairchem) calculator reads atoms.info and otherwise warns at every force call.  The same
    (charge 0, spin 1) setting is used for the atomic reference energies — identical to the fairchem
    default used in every completed run — so the raw atomisation route is model-defined and the
    element-wise constants cancel exactly in the atom-equivalent route."""
    from ase import Atoms
    a = Atoms(symbols, positions=positions); a.info.update(charge=0, spin=1)
    return a

def atom_energies(calc, elements):
    E = {}
    for el in elements:
        a = new_atoms(el, [[0, 0, 0]]); a.calc = calc
        try: E[el] = float(a.get_potential_energy())*EV
        except Exception: E[el] = None
    return E

def optimise(calc, symbols, pos, fmax=0.01, steps=800):
    from ase.optimize import LBFGS
    a = new_atoms(symbols, pos); a.calc = calc
    opt = LBFGS(a, logfile=None); converged = opt.run(fmax=fmax, steps=steps)
    return a, float(a.get_potential_energy())*EV, float(np.abs(a.get_forces()).max()), bool(converged)

def thermo(calc, atoms, symmetrynumber=1):
    from ase.vibrations import Vibrations
    from ase.thermochemistry import IdealGasThermo
    d = tempfile.mkdtemp()
    vib = Vibrations(atoms, delta=0.01, name=os.path.join(d, "vib")); vib.run()
    en = vib.get_energies(); vib.clean(); shutil.rmtree(d, ignore_errors=True)
    en = en[np.argsort(np.abs(en))]
    nat = len(atoms)
    if nat == 1: geometry, nskip = "monatomic", 0
    else:
        p = atoms.get_positions() - atoms.get_positions().mean(0); s = np.linalg.svd(p, compute_uv=False)
        geometry = "linear" if (s[1] < 1e-3*s[0]) else "nonlinear"; nskip = 5 if geometry == "linear" else 6
    tail = en[nskip:]
    imag = (np.abs(tail.imag) > 1e-6) | (tail.real < 0)
    vib_e = np.abs(tail.real)[~imag]
    expected = {"monatomic": 0, "linear": 3*nat - 5, "nonlinear": 3*nat - 6}[geometry]
    if len(vib_e) < expected:                      # residual imaginary / near-zero modes (free internal rotors): keep the mode count
        vib_e = np.abs(tail)[:expected]            # |iν| -> ν (complex magnitude)
    FLOOR = 50.0/8065.54                           # low-frequency floor (50 cm^-1) so that near-zero torsions do not produce NaN enthalpies
    n_floor = int((vib_e < FLOOR).sum()); vib_e = np.where(vib_e < FLOOR, FLOOR, vib_e)
    th = IdealGasThermo(vib_energies=vib_e, potentialenergy=0.0, atoms=atoms, geometry=geometry, symmetrynumber=symmetrynumber, spin=0)
    zpve = th.get_ZPE_correction()*EV
    H = th.get_enthalpy(temperature=298.15, verbose=False)*EV
    return dict(zpve=zpve, H298_minus_E=H, thermal=H - zpve, n_imag=int(imag.sum()), n_floor_modes=n_floor,
                lowest_cm=float(vib_e.min()*8065.54) if len(vib_e) else None, geometry=geometry,
                freqs_cm=[round(float(x)*8065.54, 1) for x in vib_e])

def dHf_from_record(rec, Eat):
    c = rec["counts"]; els = [e for e in ELEMENTS if c.get(e, 0)]
    Esum = sum(c[e]*Eat[e] for e in els)
    D0 = Esum - rec["E_model_kcal"] - rec["zpve"]
    dHf0 = sum(c[e]*DHF0_ATOM[e] for e in els) - D0
    return dHf0 + rec["H298_minus_E"] - sum(c[e]*H298_ELEM[e] for e in els), D0, dHf0

# ------------------------------------------------------------------ main loop
def run(model, inp, out_dir, k=5, fmax=0.01, device=None, limit=None, time_budget_h=None, checkpoint_every=50, checkpoint_zip=None):
    """time_budget_h: stop cleanly after this many hours (resume later from the JSONL/zip); checkpoint_zip: path of a zip rewritten every `checkpoint_every` molecules."""
    import torch
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    fj = out/f"{model.replace('/', '_').replace(':', '_')}.jsonl"
    done = set()
    if fj.exists():
        for line in open(fj):
            try: done.add(json.loads(line)["id"])
            except Exception: pass
    df = pd.read_csv(inp)
    if limit: df = df.head(limit)
    calc = make_calc(model, device)
    Eat = atom_energies(calc, ELEMENTS)
    print(f"[{model}] device={device} | atom energies (kcal/mol): { {k: (round(v,1) if v else None) for k,v in Eat.items()} }", flush=True)
    t_all = time.time(); n_ok = 0; n_since_ckpt = 0; stopped_early = False
    def checkpoint():
        if not checkpoint_zip: return
        import zipfile
        with zipfile.ZipFile(checkpoint_zip, "w", zipfile.ZIP_DEFLATED) as z:
            for pth in out.glob("*"):
                if pth.is_file(): z.write(pth, pth.name)
    with open(fj, "a") as fo:
        for _, r in df.iterrows():
            if r["id"] in done: continue
            if time_budget_h and (time.time() - t_all) > time_budget_h*3600:
                stopped_early = True; print(f"[{model}] time budget of {time_budget_h} h reached — stopping cleanly; re-run with the results zip attached to resume", flush=True); break
            if "evaluate" in r and str(r["evaluate"]).lower() == "false": continue
            t0 = time.time()
            rec = dict(id=r["id"], name=r.get("name", ""), smiles=r["smiles"], exp=float(r["dHf_exp_kcal"]), model=model)
            try:
                from rdkit import Chem
                els = {a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(r["smiles"])).GetAtoms()}
                if model in SUPPORTS and not els <= SUPPORTS[model]:
                    rec["skipped"] = f"elements {sorted(els - SUPPORTS[model])} not supported"; fo.write(json.dumps(rec)+"\n"); fo.flush(); continue
                confs, nrot, nemb = conformers(r["smiles"], k)
                best = None; conf_E = []
                for symbols, pos, emmff in confs:
                    a, E, fres, conv = optimise(calc, symbols, pos, fmax=fmax); conf_E.append(E)
                    if best is None or E < best[1]: best = (a, E, fres, conv, emmff)
                a, E, fres, conv, emmff = best
                th = thermo(calc, a); tries = 0
                while th["n_imag"] > 0 and tries < 2:
                    tries += 1; rng = np.random.default_rng(tries)
                    a2, E2, f2, c2 = optimise(calc, a.get_chemical_symbols(), a.get_positions() + rng.normal(0, 0.05, a.positions.shape), fmax=fmax/2, steps=1500)
                    th2 = thermo(calc, a2)
                    if th2["n_imag"] <= th["n_imag"]: a, E, fres, conv, th = a2, E2, f2, c2, th2
                counts = {el: a.get_chemical_symbols().count(el) for el in ELEMENTS if el in a.get_chemical_symbols()}
                rec.update(counts=counts, n_heavy=sum(v for e, v in counts.items() if e != "H"), n_atoms=len(a),
                           E_model_kcal=E, fmax_final=fres, converged=conv, n_conf_embedded=nemb, n_conf_opt=len(confs), n_rot=nrot,
                           conf_energy_spread=float(max(conf_E) - min(conf_E)), mmff_rank_of_best=int(np.argmin(conf_E)),
                           reopt_tries=tries, positions=np.round(a.get_positions(), 5).tolist(), symbols=a.get_chemical_symbols(), **th)
                rec["atom_E_kcal"] = {e: Eat[e] for e in counts}
                rec["dHf298_raw"], rec["D0"], rec["dHf0"] = dHf_from_record(rec, Eat)
                rec["error_raw"] = rec["dHf298_raw"] - rec["exp"]; rec["seconds"] = time.time() - t0; n_ok += 1
                print(f"[{model}] {str(rec['id']):>6s} {str(rec['name'])[:18]:18s} nheavy={rec['n_heavy']:2d} exp={rec['exp']:8.2f} raw={rec['dHf298_raw']:9.2f} imag={th['n_imag']} conv={conv} ({rec['seconds']:.0f}s)", flush=True)
            except Exception as e:
                rec["error"] = repr(e); print(f"[{model}] {r['id']} FAILED: {e!r}", flush=True)
            fo.write(json.dumps(rec) + "\n"); fo.flush(); n_since_ckpt += 1
            if n_since_ckpt >= checkpoint_every: checkpoint(); n_since_ckpt = 0
    checkpoint()
    print(f"[{model}] finished {n_ok} new molecules in {(time.time()-t_all)/60:.1f} min" + (" (STOPPED EARLY — resume needed)" if stopped_early else " (complete)"), flush=True)
    # light summary (raw route only; atom equivalents are fitted in thermo_analysis.py)
    recs = [json.loads(l) for l in open(fj)]
    ok = [x for x in recs if "dHf298_raw" in x]
    pd.DataFrame([{k: v for k, v in x.items() if k not in ("positions", "symbols", "freqs_cm")} for x in ok]).to_csv(out/f"{fj.stem}_summary.csv", index=False)
    return fj

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True); ap.add_argument("--input", required=True); ap.add_argument("--out", default="results")
    ap.add_argument("--k", type=int, default=5); ap.add_argument("--fmax", type=float, default=0.01); ap.add_argument("--device", default=None); ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(); run(a.model, a.input, a.out, k=a.k, fmax=a.fmax, device=a.device, limit=a.limit)
