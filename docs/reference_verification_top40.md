# Verification of the 40 most strongly contradicted reference values (18 Oct 2026)

Scope: the 40 molecules of `family_consensus_suspect_references.csv` with the largest |ωB97M-family median error| (UMA-m, AIMNet2,
MACE-OFF23 medium; ANI-1ccx and ANI-2x agree in sign for all of them). Independent evidence: NIST Chemistry WebBook (SRD 69, gas-phase
thermochemistry pages, read 18 Oct 2026), the ATcT 1.222 subset already in hand (none of the 40 is in ATcT), simple group-additivity /
isodesmic estimates, and a name→structure audit of the whole benchmark with OPSIN 2.8.0 (`name_structure_audit.csv`).
Values in kcal mol⁻¹ (kJ mol⁻¹ in parentheses where quoted from a source). "Models" = compilation value + ωB97M-family median error.

Verdict codes: **C** compilation value wrong (independent experiment or additivity agrees with the models) · **S** structure–value mismatch in
the benchmark (wrong tautomer / isomer / connectivity / misassigned value) · **E** independent experiment agrees with the compilation — model
error or unresolved · **N** no independent gas-phase value; single-source handbook/estimate (Yaws, Winget, "Add Lit") — cannot be adjudicated.

| id | molecule | compilation (source) | models | independent evidence | verdict |
|---|---|---|---|---|---|
| B0266 | 1,1,1,2-tetrafluoroethane (HFC-134a) | −166.6 (Sandler) | −214.5 | not in NIST WebBook; literature value ≈ −214 (−895.8 ± 1.5 kJ, Kolesov & Papina 1983; G3 −896, Haworth et al. 2000) [verify both] | C |
| B1577 | 1H,2H(?)-decafluorocyclohexane C6H2F10 | −462.7 (Cox-1) | −503.6 | not in NIST; perfluoro class, see note 1 | N (perfluoro) |
| B1570 | octafluorotoluene | −303.4 (Pedley) | −343.6 | NIST: −1269.4 ± 7.8 kJ (Krech et al. 1973, reanalysed Pedley 1986) = compilation; isodesmic estimate from C6F6/C2F6/CF4 ≈ −327 | E/N (perfluoro, note 1) |
| B0136 | propyleneimine (2-methylaziridine) | −16.8 (Winget) | +20.8 | NIST: no gas value; −16.8 = ΔfH of **n-propylamine** (−70.1 kJ); aziridine +30.2 − CH3 ≈ +22 | C (misassigned value) |
| B1194 | 2,4,6-trimethylpyridine N-oxide | +32.7 (Cioslowski) | 0.0 | not in NIST; pyridine N-oxide ≈ +20, three CH3 (−8 each) → ≈ −4 | C (likely; source is a computed set) |
| B0831 | cyclohexanone oxime | −46.4 (Yaws) | −16.9 | NIST: **−74.88 kJ = −17.9** (Kozyro, Kabo et al. 1992, Ccb) | C |
| B0925 | "1,2,3-benzoxadiazole" | +71.9 (Pedley) | +50.0 | NIST files 300.7 ± 2.1 kJ (Arshadi 1980) under 1,2,3-benzoxadiazole, but the paper measured **2,1,3-benzoxadiazole (benzofurazan)**; benchmark SMILES is the 1,2,3-isomer | S (isomer) |
| B1539 | isophthalamide | −91.4 (Pedley) | −69.7 | NIST: −382 ± 0.4 kJ (Hamilton & Witt 1971) = compilation; additivity from benzamide (−100.9 kJ) → ≈ −285 kJ = −68; PM6 −294 kJ | C (1971 value implausible by additivity) |
| B0692 | 1-hydroperoxyhexane | −80.5 (Cox-1) | −58.8 | NIST: no gas value; ROH→ROOH increment (+70 kJ) from 1-hexanol (−316 kJ) → −246 kJ = −58.8 | C |
| B1546 | phthalamide | −89.8 (Pedley) | −68.3 | NIST: no thermochemistry; same additivity as isophthalamide → ≈ −68 | C |
| B0849 | "dimethylglyoxime" | −19.3 (Cox-1) | +1.8 | benchmark SMILES CO/N=C(C)/C=N/O is an **O-methyl oxime isomer**, not dimethylglyoxime (which is B0710, same value, models +8.2) | S (connectivity); B0710 remains E/N |
| B0915 | 1,1,1,3,3,3-hexafluoropropane | −310.6 (Yaws) | −331.4 | NIST: no gas value | N |
| B0724 | tropone | −4.2 (Pedley) | +15.8 | NIST: ΔfH°liq ≈ −10 kJ (ΔcH −3602 ± 3, Jackson et al. 1971) + ΔvapH 54.0 ± 0.4 → **gas ≈ +44 kJ = +10.5**; compilation value is close to the liquid value | C (phase error; residual +5 model) |
| B0378 | bis(hydroxymethyl) peroxide | −137.1 (Cox-1; Pedley) | −117.8 | NIST: −571.7 ± 6.7 kJ (Jenkins & Style 1953, reanalysed Cox & Pilcher) = compilation; Benson groups → ≈ −494 kJ = −118; PM6 −521 kJ | C (1953 value implausible by additivity) |
| B1221 | 2-fluoro-2,2-dinitroethanol | −114.1 (Add Lit) | −95.7 | energetic-materials literature only | N |
| B1583 | tetranitrocyclobutane-1,3-dione | +9.0 (Add Lit) | −7.7 | energetic-materials literature only; ANI-1ccx −31 | N |
| B0148 | 1,2-difluoropropane | −132.2 (Yaws) | −115.7 | not in NIST | N |
| B1582 | perfluorocyclohexane | −566.4 (Cioslowski; Pedley) | −582.3 | NIST: **−2370.4 ± 8.4 kJ** (Price & Sapiano 1979, reanalysed Pedley 1986) = compilation | E (perfluoro, note 1) |
| B0834 | diacetone alcohol | −130.4 (Yaws) | −114.6 | NIST: −540.7 kJ = −129.2 (Vilcu et al. 1975, Ccb); Benson ≈ −121 to −123; PM6 −486 kJ | E/unresolved (experiment, additivity and models disagree) |
| B1053 | cyclohexyl azide | +36.9 (Cox-1; Pedley) | +52.7 | NIST: 154.4 kJ from ΔfH°liq 108.4 (Fagley & Myers 1954) + an **estimated** ΔvapH 46 kJ (citation missing) | E/unresolved (old liquid-phase combustion of an azide) |
| B0603 | butyl isocyanate | −30.4 (Yaws) | −45.8 | NIST: no gas value | N |
| B0582 | "4-pyridone" | −9.8 (Pedley) | +5.4 | benchmark SMILES O=C1C=CN=CC1 is a **non-aromatic tautomer**; NIST gas values −30.3 ± 5.3 / −40.8 ± 2.1 kJ refer to 4-hydroxypyridine/4-pyridone | S (tautomer) |
| B1113 | 1,1,1-trifluoroheptane | −187.6 (Yaws) | −202.4 | not checked (Yaws estimate) | N |
| B1571 | tetranitrocyclopropene | +102.4 (Add Lit) | +87.6 | energetic-materials literature only | N |
| B1310 | 4-methyl-2-isopropyl-1-pentanol | −107.0 (Yaws) | −92.6 | Yaws gives −107.0 and −81.6 for two C9H20O isomers that should differ by ≈ 1–2; models −92.6 / −92.4; 1-nonanol ≈ −90 | C (both Yaws values) |
| B1576 | DDCP (tetranitro-dimethylenecyclopropane) | +130.5 (Add Lit) | +116.4 | energetic-materials literature only | N |
| B0962 | "2,3-diethoxypropane" (1,2-diethoxypropane) | −121.2 (Cox-1) | −107.3 | NIST: not checked; 1,2-dimethoxypropane-type additivity ≈ −105 to −110 | C (likely) |
| B1083 | hexafluoropropene | −258.1 (Yaws) | −272.1 | NIST: **−1151.7 kJ = −275.3** (Papina & Kolesov 1987) | C |
| B0762 | squaric acid | −123.0 (Pedley) | −109.9 | NIST: −515 ± 17 kJ (Sellers 1971) = compilation, ±4 kcal uncertainty | E/unresolved |
| B0381 | butyl isocyanide | +9.6 (Yaws) | +22.2 | Yaws gives the **same value 9.56 for n-, sec- and tert-butyl isocyanide**; not in NIST | N (estimate) |
| B1345 | adenine | +48.9 (Pedley) | +61.4 | benchmark SMILES is the **7H-tautomer** (9H is ≈ 7–8 lower in the gas phase); NIST gas value absent on the species page; literature exp 205 kJ vs G3/G4 ≈ 220 kJ | S (tautomer) + residual |
| B1554 | N,N-difluoro-1,3-dinitroazetidin-3-amine | +7.0 (Add Lit) | +19.1 | energetic-materials literature only | N |
| B1054 | cyclohexyl isocyanate | −33.5 (Yaws) | −45.5 | not in NIST | N |
| B0465 | 1,3-difluoro-2-(fluoromethyl)propane | −172.0 (Yaws) | −160.3 | not checked (Yaws estimate) | N |
| B0274 | 1,2,3,6-tetrahydropyridine | +7.2 ("NIST-1") | +18.3 | NIST lists only the **PM6 value 30 kJ (Stewart 2004)** — the compilation's "NIST" value is semi-empirical; piperidine (−47.2 kJ) + C=C (+118) ≈ +71 kJ = +17 | C (PM6 value, not experiment) |
| B0369 | nitroethanol | −72.3 (Stewart; Winget) | −61.3 | NIST: no gas value; name ambiguous (1- vs 2-nitroethanol) | N |
| B1311 | 4-methyl-2-propyl-1-pentanol | −81.6 (Yaws) | −92.4 | see B1310 | C |
| B0824 | cis-4-methylcyclohexanol | −85.4 (Cox-1; Pedley) | −74.7 | NIST: no gas value; cyclohexanol (−286.2 kJ) + CH3 (−33) ≈ −76 | C (likely a liquid-phase value) |
| B0306 | 2,2-difluorobutane | −125.9 (Yaws) | −136.6 | not checked (Yaws estimate) | N |
| B0374 | bicyclo[1.1.0]butane-1-carbonitrile | +72.8 (Pedley) | +83.5 | NIST: 304.5 kJ (Hall & Baldt 1971, reaction calorimetry) = compilation; bicyclobutane (217 kJ) + CN increment (≈ +130) ≈ +347 kJ = +83 | C (likely; additivity agrees with models) |

Tally (40): C or S — 20 (compilation/benchmark errors confirmed or strongly indicated); N — 15 (unverifiable single-source estimates, mostly
Yaws / "Add Lit" energetic materials); E/unresolved — 5 (perfluorocyclohexane, octafluorotoluene, diacetone alcohol, cyclohexyl azide,
squaric acid). Among the first fifteen, eleven are C or S.

Note 1 — perfluoro compounds. The perfluoroalkane/-cycloalkane values rest on fluorine-/rotating-bomb combustion with HF(aq) as product; NIST
shows that the same measurement gives −3301 (original, Good et al. 1959) or −3385.4 ± 3.6 kJ (reanalysed, Kolesov & Kozina 1986) for
perfluoroheptane depending on the auxiliary data, and the compilation's values imply a CF2 increment of −126 kcal between C6F12 and C7F14
where ≈ −100 to −105 is expected. The ωB97M-family models deviate by −16 (C6F12) and +11 to +18 (C7F14, C7F16, C8F16, C10F22O) — opposite
signs for molecules that differ by one CF2 — so neither side can be adjudicated with the present data; the paper should report this class
separately and not count it in the headline.

Note 2 — structure–value mismatches found in the whole benchmark (OPSIN audit + manual checks; `name_structure_audit.csv`,
`tautomer_check.csv`): B0055 formimidic acid (SMILES N=CO carries **formamide's** value −186 kJ; models +9 to +11 = the imidic/amide gap),
B0582 4-pyridone (non-aromatic tautomer), B0925 benzoxadiazole (wrong isomer), B0849 dimethylglyoxime (wrong connectivity, duplicate of
B0710), B1345 adenine (7H instead of 9H), B0646 "isopropyl acetate" (SMILES is propyl acetate, value is isopropyl acetate's −481.6 kJ),
B1480 cumene hydroperoxide (SMILES puts –OOH on the ring instead of the α-carbon), B0900 tetramethyldiazetidine (ambiguous isomer),
B1109 trans-bicyclo[3.3.0]octan-2-one (value probably for the cis isomer). Recommended for benchmark v1.1: set `evaluate = False` with
`exclude_reason = "structure–value mismatch"` for these nine (B0710 stays), and re-derive the headline (the effect on the MAEs is < 0.05).

Note 3 — the compilation's "NIST-1" column contains Stewart 2004 PM6 values wherever the WebBook lists only a semi-empirical entry
(1,2,3,6-tetrahydropyridine is one case); single-source "NIST-1" rows should be re-checked for the "Semi" flag before being used as experiment.
