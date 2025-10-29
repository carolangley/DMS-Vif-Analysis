# Site-Level Analyses (Clean CLI)

This repo converts your **site level analyses** notebook into a reproducible, modular CLI. It preserves your logic (background thresholding, designed-variant filtering, pooled-pre reference, WT centering, Wilcoxon+BH-FDR, replicate-consistency gating, conservation overlays, and cross-strain Mann–Whitney).

## Install
```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn
```

## Typical Run

1) **Prepare** (background → filtered codoncounts → pooled-pre tables)
```bash
python site_level_analyses_clean.py prepare
```

2) **Enrich** (build `counts_df.csv`, `enrich_df.csv`)
```bash
python site_level_analyses_clean.py enrich
```

3) **Replicate QC** (correlation grids + sign agreement as PDFs)
```bash
python site_level_analyses_clean.py replicate-qc
```

4) **Site stats** (site violins; Wilcoxon vs global median + BH-FDR; replicate-consistent lists)
```bash
python site_level_analyses_clean.py site-stats
```

5) **Conservation overlays** (median site enrichment vs %identity / Shannon entropy)
```bash
python site_level_analyses_clean.py conservation --alignment subtype_ref_protein.fasta --mode overlay
```

6) **Cross-strain** (LAI vs 1203 medians; export significant sites + plot)
```bash
python site_level_analyses_clean.py cross-strain
```

### Inputs expected
- `consensuseachUMI/filtered/*_codoncounts.csv`
- `LAIvariantproportions_with_wildtype.csv`, `1203variantproportions_with_wildtype.csv`
- (for conservation) `subtype_ref_protein.fasta`

### Outputs (high level)
- `qc_outputs/background_thresholds.csv`
- `consensuseachUMI/filteredcodoncounts/` filtered per-replicate + pooled-pre CSVs
- `counts_df.csv`, `enrich_df.csv`
- `figures/*.pdf` (QC grids, violins, site-median comparison)
- `site_significance_full.csv`, `replicate_consistent_BH_calls.csv`
- `significant_site_differences.xlsx`
- `supplementary/*svg` (conservation overlays)

## Notes
- WT centering is applied by default for stability: set `USE_WT_CENTERING` inside the script if you want raw log enrichment.
- The pipeline assumes replicate file names like `preLAIR1`, `postLAIR2`, etc., and creates pooled `mergedpreLAI/1203` tables.
- Figures use a non-interactive backend (`Agg`) so everything runs headless.

## License
MIT © 2025 Caroline Langley
