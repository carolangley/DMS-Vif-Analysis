# Hypermutation Read Extractor

Extracts reads from FASTQ files that carry specified codon mutations defined in an enrichment CSV. For each `(site, codon, wildtype, aa_mutation)` that matches a given sample, the script writes matching reads to a per‑mutation FASTA file.

This repository contains a cleaned, publication‑ready script with a consistent CLI, type hints, and structured logging.

---

## Requirements

- Python 3.8+
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)

Install via pip:
```bash
pip install biopython pandas
```

## Inputs

- **Enrichment CSV** (`--enrich-csv`): must include columns
  - `site` (1‑based codon index), `codon` (mutant triplet), `wildtype` (WT triplet), `sample_post` (sample identifier)
  - Optional: `aa_mutation` (string; used only in output file names)

- **FASTQ files**: one or more input reads (`.fastq` or `.fastq.gz`).
  - The sample id is derived from the FASTQ stem (suffix `_filtered` is stripped). Rows are selected where `sample_post` **contains** this sample id.

## Output

For each matching mutation, a FASTA file is written to `--out-dir`:
```
<sample>_<codon>_<wildtype>_<site>_<aa_mutation>.fasta
```

## Usage

```bash
python hypermutation_clean.py   --enrich-csv enrich_Y40Xonly40_df.csv   --out-dir hypermutationY40   --site-min 12   --site-max 115   --buffer-size 10000   --log-level INFO   consensuseachUMI/filtered/preY40XR1_filtered.fastq   consensuseachUMI/filtered/postY40XR1_filtered.fastq.gz
```

### CLI Options

- `--enrich-csv PATH` (required): enrichment CSV path
- `--out-dir PATH` (required): directory for FASTA outputs (created if needed)
- `--site-min INT` (default: 12): minimum codon site (inclusive)
- `--site-max INT` (default: 115): maximum codon site (inclusive)
- `--buffer-size INT` (default: 10000): per‑mutation write buffer
- `--log-level {DEBUG,INFO,WARNING,ERROR}` (default: INFO)
- `FASTQ ...`: one or more `.fastq` or `.fastq.gz` files

## Notes

- Only positions within `[site-min, site-max]` and where `codon != wildtype` are evaluated.
- Matching is an exact triplet comparison at the specified site.
- Gzipped FASTQ is supported transparently.
- If no enrichment rows match a sample, that FASTQ is skipped with a warning.

## Reproducibility

- Python and package versions should be recorded (e.g., `pip freeze > requirements.txt`).
- Keep the enrichment CSV committed or archived with the manuscript submission.

## Citation

If you use this code in a publication, please cite the associated manuscript.

## License

Add an appropriate license (e.g., MIT) before public release.
