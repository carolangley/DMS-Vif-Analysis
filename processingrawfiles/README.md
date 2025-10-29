# Processing Raw FASTQ Files

This repository provides a cleaned, publication-ready script for merging and trimming paired-end FASTQ files and building UMI-based consensus reads.

## Contents
- `processingrawfiles_clean.py` — CLI tool with two subcommands:
  - `merge-trim` — merge paired reads with **BBMerge** and trim/length-filter.
  - `consensus` — build per-UMI consensus reads from merged/trimmed FASTQ files.
- `requirements.txt`
- `.gitignore`
- `LICENSE` (MIT)

## Installation
```bash
pip install biopython tqdm
# BBMerge must be installed and `bbmerge.sh` available on PATH for the merge step.
```

## Usage

### 1) Merge and trim
```bash
python processingrawfiles_clean.py merge-trim   fastqfiles/   consensuseachUMI/merged   --min-length 384   --max-length 384   --workers 8
```

- Expects files named like `SAMPLE_R1_...fastq.gz` and `SAMPLE_R2_...fastq.gz` in `fastqfiles/`.
- Outputs `SAMPLE_merged.fastq.gz` and `SAMPLE_merged_trimmed.fastq.gz` in `consensuseachUMI/merged/`.

### 2) UMI consensus from merged/trimmed reads
```bash
python processingrawfiles_clean.py consensus   consensuseachUMI/merged   consensuseachUMI/consensus   --umi-length 9   --workers 8   --min-gc 20 --max-gc 80   --max-n 10   --min-q 20
```

- Expects `*_merged_trimmed.fastq.gz` in the input folder.
- Writes `*_consensus.fastq.gz` to the output folder.

## Notes
- Logging level can be adjusted with `--log-level` on either subcommand (`DEBUG`, `INFO`, `WARNING`, `ERROR`).
- The merge step requires **BBMerge** (`bbmerge.sh`). Install from the BBMap suite.
- The consensus logic groups reads by 5' and 3' UMIs (default 9 nt each) and applies GC/N/quality filters prior to basewise voting.

## License
MIT © 2025 Caroline Langley
