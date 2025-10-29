# Hypermutation Pipeline (Calculate → Analyze)

This repository packages your hypermutation workflow into a single, clean CLI that runs in the right order:
first **calculate** metrics from per-mutation FASTA files, then **analyze** those metrics to generate plots.

## Install
```bash
pip install pandas numpy scipy matplotlib biopython scikit-learn
```

## Usage

### 1) Calculate hypermutation metrics
```bash
python hypermutation_pipeline_clean.py calculate   --fasta-dir hypermutation/   --enrich-csv enrich_df.csv   --out-dir .   --lai-effects-csv LAI_GGtoAGmutation_effects.csv   --s1203-effects-csv 1203_GGtoAGmutation_effects.csv   --workers 22
```
This writes:
- `LAI_GGtoAGmutation_effects.csv` and `1203_GGtoAGmutation_effects.csv`
- `LAI_hypermutation_results.csv` and `1203_hypermutation_results.csv`

### 2) Analyze and plot
```bash
python hypermutation_pipeline_clean.py analyze   --results-dir .   --varprops-dir .   --out-dir figures
```
This writes:
- `figures/LAI_hypermutation.pdf`, `figures/1203_hypermutation.pdf`
- `figures/LAI_hypermutation_linear.pdf`, `figures/1203_hypermutation_linear.pdf`
- `figures/LAI_overall_fitness.pdf`, `figures/1203_overall_fitness.pdf`

## Notes
- The pipeline mirrors your original logic but is explicit, reproducible, and server-safe (non-interactive plotting).
- File naming for FASTA is expected as in your prior workflow: `{pre|post}{LAI|1203}R{rep}_filtered_{codon}_{wildtype}_{site}_{aa_mut}.fasta`.
- The **calculate** stage restricts to files containing `post` in the name (as before).
- Variant-proportion CSVs should be available as: `LAIvariantproportions_with_wildtype.csv` and `1203variantproportions_with_wildtype.csv`.

## License
MIT © 2025 Caroline Langley
