#!/usr/bin/env python3
"""Build InChIKey-skeleton (first 14 chars) indices of training sets for the leakage audit.
Outputs: audit_index/<set>.csv with columns: skeleton, canonical_smiles, source_id"""
import h5py, glob, tarfile, csv, sys
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
RDLogger.DisableLog("rdApp.*")
OUT = Path("audit_index"); OUT.mkdir(exist_ok=True)

def skel(mol):
    try:
        return inchi.MolToInchiKey(mol)[:14]
    except Exception:
        return None

def clean(smi):
    """strip atom maps, return list of canonical component SMILES (largest first)"""
    m = Chem.MolFromSmiles(smi)
    if m is None: return []
    for a in m.GetAtoms(): a.SetAtomMapNum(0)
    frags = Chem.GetMolFrags(m, asMols=True)
    return sorted(frags, key=lambda f: -f.GetNumHeavyAtoms())

# ---------------- SPICE (openmm/spice-dataset inputs) ----------------
def build_spice():
    rows = {}
    files = sorted(glob.glob("spice-dataset/**/*.hdf5", recursive=True))
    for f in files:
        sub = f.split("/")[1]
        with h5py.File(f) as h:
            for k in h.keys():
                s = h[k]["smiles"][()]
                s = s[0].decode() if hasattr(s, "__len__") and not isinstance(s, str) else str(s)
                for frag in clean(s):
                    sk = skel(frag)
                    if sk and sk not in rows:
                        rows[sk] = (Chem.MolToSmiles(frag), f"{sub}:{k}")
    with open(OUT/"spice_v1_inputs.csv", "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["skeleton", "canonical_smiles", "source_id"])
        for sk, (smi, src) in rows.items(): w.writerow([sk, smi, src])
    print("SPICE skeletons:", len(rows), "from", len(files), "hdf5 files")

# ---------------- QM9 (dsgdb9nsd.xyz.tar.bz2) ----------------
def build_qm9(tar_path):
    rows = {}; n = 0; bad = 0
    with tarfile.open(tar_path, "r:bz2") as t:
        for mem in t:
            if not mem.name.endswith(".xyz"): continue
            n += 1
            lines = t.extractfile(mem).read().decode("utf8", "ignore").splitlines()
            try:
                smi_line = lines[-2].split()   # "SMILES_gdb  SMILES_relaxed"
                smi = smi_line[1] if len(smi_line) > 1 else smi_line[0]
                m = Chem.MolFromSmiles(smi)
                if m is None: bad += 1; continue
                sk = skel(m)
                if sk and sk not in rows: rows[sk] = (Chem.MolToSmiles(m), mem.name.replace(".xyz", ""))
            except Exception:
                bad += 1
    with open(OUT/"qm9.csv", "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["skeleton", "canonical_smiles", "source_id"])
        for sk, (smi, src) in rows.items(): w.writerow([sk, smi, src])
    print(f"QM9: {n} files, {bad} unparsable, {len(rows)} unique skeletons")

if __name__ == "__main__":
    which = sys.argv[1]
    if which == "spice": build_spice()
    elif which == "qm9": build_qm9(sys.argv[2])
