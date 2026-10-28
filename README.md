# hf298-mlip-benchmark

Data, protocol, raw results and analysis for

> A. A. Khairbek, R. Thomas, et al. *Pretrained machine-learning potentials as thermochemistry engines: a physically complete,
> leakage-audited benchmark of gas-phase enthalpies of formation*, submitted to *J. Phys. Chem. A* (2026).

Verified gas-phase ΔfH°(298.15 K) of 1 570 closed-shell H/C/N/O/F molecules, a physically complete thermochemistry protocol for
pretrained machine-learning potentials, the raw per-molecule records of ten model runs, and the scripts that regenerate every table and
figure of the paper. The tables and figures themselves are **not** stored here; `analysis/run_all.sh` rebuilds them in about a minute.

## Layout

```
benchmark_verified.csv          the benchmark (1 589 rows, 1 570 evaluated; one row per molecule, see benchmark_README.md)
benchmark_README.md             column definitions, tiers, exclusions, version history
run_thermo_protocol.py          the protocol: conformers -> model optimisation -> frequencies -> RRHO -> atomisation route
benchmark_pipeline_kaggle.ipynb one-notebook-per-model runner (Kaggle GPU, time budget, checkpoint/resume) that embeds the protocol
raw_results/<model>.jsonl       raw per-molecule records of the ten runs (energies, ZPVE, thermal terms, frequencies, geometry, timing)
analysis/                       post-processing, evaluation, figures (run_all.sh)
build/                          scripts and inputs that built benchmark_verified.csv (not needed to reproduce the paper)
docs/                           manual verification of the forty largest contradicted reference values
```

## Reproduce the paper

```bash
pip install -r requirements.txt        # numpy, pandas, matplotlib, rdkit, ase (CPU only)
bash analysis/run_all.sh
```

`run_all.sh` runs, in order,

| script | what it does | writes |
|---|---|---|
| `postprocess_runs.py` | uniform thermochemistry from the raw records: 50 cm⁻¹ floor, imaginary-mode handling, RRHO H(298)−H(0), atomisation enthalpy with ATcT atomic values | `analysis/output/<model>_merged.csv` |
| `evaluate_models.py` | five element-wise atom equivalents by 10-fold cross-validation (identical folds for every model, seed 2026) and leave-one-out; bootstrap 95 % CIs; paired model comparisons; error correlations; 37 subsets (tiers, size, chemistry, leakage, charge class) | `analysis/output/eval/summary_long.csv`, `headline.csv`, `paired_comparisons.csv`, `errors_cv_wide.csv`, `<model>_evaluated.csv`, … |
| `analysis_extra.py` | correlations by tier, family-consensus list of contradicted reference values, ensembles, matched leakage tests, ZPVE-consistency and energy/vibration hybrid checks | `leakage_tests.csv`, `family_consensus_suspect_references.csv`, `zpve_consistency.csv`, `hybrid_energy_vibration_check.csv`, … |
| `make_tables.py` | Tables 1–3 of the paper | `table1_headline.csv`, `table2_subsets_*.csv`, `table3_error_by_reference_provenance.csv` |
| `make_figures.py` | Figures 1–5 (PNG + PDF, 300 dpi) | `analysis/output/figures/` |

Every number in the paper comes from these outputs; the run is deterministic (fixed seeds), and differences between machines are at the
10⁻¹⁴ kcal/mol level (floating-point summation order).

## Run the protocol on a model

```bash
python run_thermo_protocol.py --model mace-off23-medium --input benchmark_verified.csv --out results/
```

Models: `uma-s`, `uma-m` (fairchem, OMol task; gated — `huggingface-cli login` first), `mace-off23-small`, `mace-off23-medium`,
`aimnet2`, `egret-1`, `egret-1e` (Rowan TorchScript models, auto-downloaded), `ani-2x`, `ani-1ccx` (HCNO only), `xtb-gfn2`.
Per-model installs are in the notebook's INSTALL cell. Each molecule is written to `results/<model>.jsonl` as soon as it finishes, so
an interrupted run resumes where it stopped. The ten runs here took 2–31 GPU-hours each (Table 1 of the paper, column GPU-h).

Protocol in one line: RDKit ETKDGv3 conformers (10 + 4·n_rot, ≤ 50) → MMFF94 → the five lowest re-optimised with the model
(ASE LBFGS, fmax 0.01 eV/Å) → finite-difference Hessian (δ = 0.01 Å) → up to two re-optimisations along imaginary modes → RRHO with a
50 cm⁻¹ floor → atomisation route with ATcT atomic ΔfH°(0 K) and elemental H(298)−H(0) → five atom equivalents fitted by cross-validation.

## Rebuild the benchmark (optional)

Needs `Measured_Enthalpy_V2.7.xlsx` from the Bains et al. compilation (Zenodo [10.5281/zenodo.4661783](https://doi.org/10.5281/zenodo.4661783)):

```bash
cd build
python build_benchmark.py Measured_Enthalpy_V2.7.xlsx audit_index/ stage1.csv          # HCNOF, closed-shell, neutral; median consensus; tiers; stereo; leakage flags
python apply_atct_and_classes.py stage1.csv atct_v1222_gas_HCNOF_subset.psv atct_matches.csv ../benchmark_verified.csv   # ATcT override, charge classes, exclusions
```

`audit_index/` holds the InChIKey-skeleton indices of QM9 and SPICE v1 (`build_audit_index.py` shows how they were made from the
original files). `audit_names_opsin.py` (needs Java and the OPSIN jar) and `tautomer_check.py` are the structure audits reported in
the paper; `atct_v1222_gas_HCNOF_subset.psv` is the ATcT 1.222 subset read on 1 Oct 2026.

## Citation and licence

Please cite the paper and the archived release of this repository (DOI in `CITATION.cff` once minted). Code is released under the MIT
licence; `benchmark_verified.csv` is derived from the Bains et al. (2022) compilation and is distributed under CC BY-SA 4.0 (see `LICENSE`).
ATcT values are from the Active Thermochemical Tables, version 1.222 (B. Ruscic and D. H. Bross, [atct.anl.gov](https://atct.anl.gov)).
