# benchmark_verified.csv — verified experimental gas-phase ΔHf°(298.15 K) test set (HCNOF)

Built by `build_benchmark.py` from the Bains, Petkowski, Zhan & Seager (2022) compilation
(*Data* **7**, 33; Zenodo 10.5281/zenodo.4661783, file `Measured_Enthalpy_V2.7.xlsx`, CC BY-SA 4.0),
which lists, per molecule, the values reported by up to 24 sources (Pedley, Cox & Pilcher, NIST WebBook, ATcT, JANAF,
Stewart, Winget, Jorgensen, Sandler, Cioslowski, Yaws, additional literature).

## Processing (all rules are in the script header)

1. HCNOF molecules only; RDKit-parsable; net-neutral; closed-shell; duplicates merged by full InChIKey (sources pooled).
2. **Consensus value = median of all source values.** The compilation's own "ATcT" column disagrees with every other source
   by 4–8 kcal/mol for ≈15 molecules (formic acid, acetic acid, carbonyl fluoride, …); the median is immune to a single wrong column.
3. **Robust spread:** with ≥ 3 sources, sources deviating > 2 kcal/mol from the median are set aside (`outlier_sources`);
   `unc_kcal` = half of the robust spread (NaN for single-source entries).
4. **Tiers:** A = ≥ 2 agreeing sources, robust spread ≤ 1 kcal/mol · B = single source or robust spread ≤ 2 · C = > 2 (kept, flagged) ·
   X = discordant (> 12 kcal/mol, no majority) → `evaluate = False`. Also excluded (`exclude_reason`): the elemental reference states and
   multi-fragment SMILES (the compilation's "pyridinium dicyanomethylide" is written as a disconnected ion pair, C8H7N3, whereas the
   real ylide is C8H5N3; every model fails on it, which is a data error, not a model error), and — v1.1 — nine structure–value mismatches
   found by an OPSIN name→structure audit and tautomer/isomer checks (formimidic acid carrying formamide's value, a non-aromatic 4-pyridone
   SMILES, "1,2,3-benzoxadiazole" carrying benzofurazan's value, a wrong dimethylglyoxime SMILES, 7H-adenine, propyl-acetate SMILES with
   isopropyl acetate's value, cumene hydroperoxide with –OOH on the ring, an ambiguous tetramethyldiazetidine, a probable cis/trans mix-up);
   see docs/reference_verification_top40.md.
5. **Stereochemistry:** where the SMILES leaves stereo elements unspecified, all stereoisomers are enumerated, embedded and MMFF94-optimised;
   if ≥ 2 distinct embeddable diastereomers exist, the lowest-energy one is assigned (`smiles`) and `stereo_assumed = True`
   (`smiles_original` keeps the source string). Sensitivity analyses should be run with and without these molecules.
6. **Strata:** `stratum` = small (≤ 9 heavy atoms, QM9-size) / large (≥ 10) for the size-extrapolation test.
7. **ATcT verification (v1.222, the current release, read from atct.anl.gov on 1 Oct 2026):** every gas-phase HCNOF ATcT species whose
   Hill formula occurs in the set (301 closed-shell entries; `atct_v1222_gas_HCNOF_subset.psv`) was matched to the benchmark by formula
   and name (unique-formula matches vetted by hand; isomer ambiguities left unmatched). For the 180 matched molecules the **ATcT value and
   uncertainty replace the compilation consensus** (`value_source = "ATcT v1.222"`, `tier = A+`); the compilation value is kept in
   `dHf_compilation_kcal` and the difference in `dev_vs_ATcT_kcal`. Finding: the compilation consensus deviates from ATcT by > 1 kcal/mol
   for 37 of 177 evaluated molecules (21 %) and by > 2 kcal/mol for 22 (12 %) — e.g. carbonyl fluoride (−7.7), ethyl hydroperoxide
   (−9.1), cyanic fluoride (+6.5), nitrosyl fluoride (+5.3), acrolein (−4.5), isocyanic acid (+4.1). ATcT covers 172 small and only
   5 large (≥ 10 heavy atoms) molecules, so the large stratum rests on Pedley / Cox–Pilcher / NIST / Roux et al. values.
8. **Training-set audit:** `in_qm9`, `in_spice_v1_inputs` = InChIKey-skeleton (connectivity) match against the QM9 (133 885 molecules) and
   SPICE v1 input (29 401 molecules) indices built by `build_audit_index.py`. Add further indices (ANI-1x/2x, OMol25, Egret-1) as CSVs
   with a `skeleton` column and re-run the builder.

## Columns

| column | meaning |
|---|---|
| id | B0001… (sorted by heavy-atom count, then name) |
| name, smiles, smiles_original, inchikey | identity; `smiles` carries the assigned stereochemistry |
| dHf_exp_kcal, dHf_exp_kJ | experimental ΔHf°(298 K, gas): ATcT v1.222 where matched, otherwise the compilation median |
| value_source | "ATcT v1.222" or "compilation median" |
| atct_name, atct_dHf_kcal, atct_unc_kcal, match, dev_vs_ATcT_kcal | ATcT match details (match = unique-formula / name / name-manual) |
| dHf_compilation_kcal, unc_compilation_kcal | the compilation consensus before the ATcT override |
| unc_kcal | half robust inter-source spread (uncertainty proxy) |
| n_sources, sources, n_sources_agreeing, outlier_sources | provenance |
| spread_kcal, spread_robust_kcal | full and robust inter-source spread |
| has_ATcT, ATcT_minus_median_kcal | whether the compilation's ATcT column exists and how far it sits from the consensus |
| tier | A+ (ATcT-verified) / A / B / C / X (see above) |
| bains_filter_flag | the compilation's own FILTER flag (spread > 50 kJ/mol) |
| n_diastereomers, stereo_assumed | stereo resolution |
| n_heavy, stratum, n_rot, n_arom_rings, n_rings, n_H … n_F | descriptors and atom counts |
| in_qm9, in_spice_v1_inputs | training-set overlap flags |
| evaluate, exclude_reason | True for the evaluation set; why a row is excluded (discordant sources / reference state / multi-fragment SMILES / structure–value mismatch / ambiguous structure) |

## Recommended reporting

* Headline metrics on the full evaluation set v1.1 (1 570 molecules: 1 096 small, 474 large); robustness on tiers A+/A only (610 molecules)
  and on `stereo_assumed == False` only.
* Size extrapolation: atom equivalents fitted on the small stratum, tested on the large stratum.
* Leakage: metrics stratified by `in_<training set>` for each model.
* Independent spot verification of the ≥ 10-heavy-atom stratum and of every outlier against the current ATcT release / Roux et al. 2008
  (PAHs) / Pedley 1994 before submission.
