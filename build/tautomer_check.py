#!/usr/bin/env python3
"""Flags benchmark SMILES whose RDKit canonical tautomer differs from the input (candidates for a wrong tautomer; stereo is ignored).
    python build/tautomer_check.py benchmark_verified.csv tautomer_check.csv"""
import sys, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")
B = pd.read_csv(sys.argv[1]); B = B[B.evaluate]; te = rdMolStandardize.TautomerEnumerator(); rows = []
for _, r in B.iterrows():
    m = Chem.MolFromSmiles(r.smiles); can = Chem.MolToSmiles(te.Canonicalize(m)); s0 = Chem.MolToSmiles(m)
    if can != s0: rows.append(dict(id=r.id, name=r["name"], smiles=r.smiles, canonical_tautomer=can, tier=r.tier, sources=r.sources))
pd.DataFrame(rows).to_csv(sys.argv[2], index=False); print(len(rows), "candidates written")
