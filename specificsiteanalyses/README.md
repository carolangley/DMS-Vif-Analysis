# Specific Sites – AA Scatterplots

Clean CLI to reproduce your per-site, per-strain scatterplots from `enrich_df.csv`.
It mirrors your original logic (grouping LAI vs 1203 replicates, computing mean±SD
log enrichment per amino acid, excluding WT matches and optionally stop `*`) while
using a non-interactive backend and explicit arguments.

## Install
```bash
pip install pandas numpy matplotlib
```

## Usage
```bash
python specific_sites_clean.py   --enrich-csv enrich_df.csv   --sites 83 70   --samples-1203 post1203R1 post1203R2 post1203R3   --samples-LAI  postLAIR1  postLAIR2  postLAIR3   --outdir figures/scatterplots   --exclude-star   --labels   --lai-only
```

### Outputs
- `figures/scatterplots/Site_<N>_LAI_vs_1203.pdf`
- (if `--lai-only`) `figures/scatterplots/Site_<N>_LAI_only.pdf`

## Notes
- Y-limits default to `[-8, 8]` for LAI vs 1203 plots and `[-6, 6]` for LAI-only;
  override with `--ylim`.
- Colors and markers: LAI = light green squares; 1203 = light purple circles.
- Files are saved as PDF for print-friendly vector graphics.
- Requires only `enrich_df.csv` (produced by your site-level or enrichment pipeline).

## License
MIT © 2025 Caroline Langley
