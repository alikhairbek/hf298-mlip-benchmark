#!/usr/bin/env bash
# Regenerates every table and figure of the paper from raw_results/*.jsonl (about one minute on a laptop).
set -e
cd "$(dirname "$0")"
python3 postprocess_runs.py     # uniform thermochemistry post-processing       -> output/<model>_merged.csv
python3 evaluate_models.py      # CV atom equivalents, bootstrap CIs, subsets    -> output/eval/
python3 analysis_extra.py       # correlations by tier, consensus list, leakage, ZPVE consistency
python3 make_tables.py          # Tables 1-3 of the paper                        -> output/eval/table*.csv
python3 make_figures.py         # Figures 1-5 (PNG + PDF, 300 dpi)               -> output/figures/
echo "done: tables in analysis/output/eval/, figures in analysis/output/figures/"
