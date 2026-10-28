#!/usr/bin/env python3
"""Name -> structure audit of the benchmark with OPSIN (Java).  Flags entries whose compound name parses to a different connectivity
(InChIKey skeleton) than the SMILES, and entries that differ only in stereo/tautomer layers.

    python build/audit_names_opsin.py benchmark_verified.csv name_structure_audit.csv
Requires Java 11+ and the OPSIN CLI jar (https://github.com/dan2097/opsin/releases, opsin-cli-<version>-jar-with-dependencies.jar)
placed next to this script as opsin.jar."""
import sys, subprocess, pandas as pd
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog("rdApp.*")
B = pd.read_csv(sys.argv[1]); B = B[B.evaluate].reset_index(drop=True)
here = Path(__file__).parent; names = here / "names.txt"; out = here / "opsin_out.txt"
names.write_text("\n".join(B["name"].astype(str)) + "\n")
subprocess.run(["java", "-jar", str(here / "opsin.jar"), "-osmi", str(names), str(out)], check=True)
res = out.read_text().split("\n"); rows = []
for (_, r), smi in zip(B.iterrows(), res):
    smi = smi.strip(); m = Chem.MolFromSmiles(smi) if smi else None; m0 = Chem.MolFromSmiles(r.smiles)
    if m is None: rows.append(dict(id=r.id, name=r["name"], parsed=False)); continue
    k1, k0 = Chem.MolToInchiKey(m), Chem.MolToInchiKey(m0)
    rows.append(dict(id=r.id, name=r["name"], parsed=True, opsin_smiles=smi, bench_smiles=r.smiles, same_skeleton=k1[:14] == k0[:14], same_full=k1 == k0,
                     opsin_formula=rdMolDescriptors.CalcMolFormula(m), bench_formula=rdMolDescriptors.CalcMolFormula(m0)))
T = pd.DataFrame(rows); T.to_csv(sys.argv[2], index=False)
P = T[T.parsed]; print(f"parsed {len(P)} of {len(T)}; connectivity mismatches {(~P.same_skeleton.astype(bool)).sum()}; stereo/tautomer-only differences {(P.same_skeleton.astype(bool) & ~P.same_full.astype(bool)).sum()}")
