# Y40 Pilot Experiments – Summary & Figure

A cleaned, publication-ready CLI that reproduces the key aggregation from your pilot experiments: it parses codon-count CSVs, computes TAG% (WT) and non‑TAT% (Y40X) at site 40, summarizes pre/post replicates, compares to input controls, runs simple t‑tests, and exports a compact PDF plus a CSV of metrics.

## Install
```bash
pip install pandas numpy scipy matplotlib
```

## Usage
```bash
python y40pilotexperiments_clean.py   --folder consensuseachUMI/filtered   --site 40   --out-pdf figures/Fig1B.pdf   --out-csv figures/Fig1B_summary.csv   --log-level INFO
```

**Expected inputs** in `--folder`:
- `*WT*codoncounts*.csv`
- `*Y40X*codoncounts*.csv`

Replicates are detected as `preWTR1..3`, `postWTR1..3`, and `preY40XR1..3`, `postY40XR1..3`. Input controls are detected by filenames containing `WTinput` and `Y40Xinput`.

### Outputs
- `figures/Fig1B.pdf` (two-panel bar chart with error bars)
- `figures/Fig1B_summary.csv` (summary metrics and p-values)

## Notes
- Non-interactive plotting backend (`Agg`) is used so this runs on servers.
- File patterns and defaults mirror your original workflow while being configurable.
- The CSV summation is robust to extra columns; it sums all codon columns at the selected site.

## License
MIT © 2025 Caroline Langley
