# Filtering Consensus Reads

A cleaned, publication-ready CLI for processing post-UMI **consensus FASTQ** files:

- **`preprocess`** — normalize read length (default 366 nt) and handle specific 5' anomalies.
- **`filter`** — keep reads within a configurable amino-acid mismatch threshold to a reference (WT / Y40X / LAI / 1203); optional GG→AG codon exceptions.
- **`count-codons`** — write per-position codon counts to CSV for filtered reads.

## Install
```bash
pip install biopython pandas
```

## Quick usage

### 1) Preprocess
```bash
python filteringconsensusreads_clean.py preprocess \
  consensuseachUMI/consensus \
  consensuseachUMI/filtered/preprocessed \
  --target-len 366
```

### 2) Filter by reference
```bash
python filteringconsensusreads_clean.py filter \
  consensuseachUMI/filtered/preprocessed \
  consensuseachUMI/filtered \
  --pattern _merged_trimmed \
  --max-aa-mismatches 4 \
  --ggtoag-csv LAI_GGtoAGmutation_effects.csv
```
- Reference is inferred from filenames containing one of: `WT`, `Y40X`, `LAI`, `1203`.
- If `--ggtoag-csv` is omitted, filtering uses pure AA-mismatch thresholding.

### 3) Count codons
```bash
python filteringconsensusreads_clean.py count-codons \
  consensuseachUMI/filtered \
  --pattern "*_filtered*.fastq.gz" \
  --out consensuseachUMI/filtered
```

## Notes
- Defaults mirror the original script’s behavior (target length 366; max 4 AA mismatches).
- Output names follow the input naming with appropriate suffixes (`_filtered`, `_codoncounts.csv`).
- For GG→AG exception lists, provide CSV(s) with columns `site,variant`.

## License
MIT © 2025 Caroline Langley
