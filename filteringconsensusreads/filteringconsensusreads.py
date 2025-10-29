#!/usr/bin/env python3
"""
filteringconsensusreads_clean.py

Three-stage toolkit for post-consensus FASTQ processing:

1) preprocess  : normalize read length and fix specific 5' anomalies (optional)
2) filter      : retain reads close to a given reference (AA mismatches threshold),
                 with optional GG→AG codon exception lists per site
3) count-codons: tally codon counts across positions for filtered FASTQ files


Requirements:
    - Python 3.8+
    - biopython
    - pandas
"""  # noqa: E501

from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# -------------------------
# Logging
# -------------------------

def init_logger(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


# -------------------------
# Preprocess
# -------------------------

def preprocess_reads(
    input_file: Path,
    output_file: Path,
    *,
    target_len: int = 366,
    add_A_if_prefix: str = "TGGAA",
    drop_first_if_prefix: str = "XATGGAA",
    pad_char: str = "X",
) -> Tuple[int, int, int]:
    """
    Normalize sequence length to `target_len` and fix specific 5' start patterns.

    Rules (matching original behavior):
      - If sequence starts with add_A_if_prefix (default "TGGAA"), prepend 'A'
      - Elif sequence starts with drop_first_if_prefix (default "XATGGAA"), drop first base
      - Then, if len > target_len, trim from 3' end to target_len
      - If len < target_len, right-pad with pad_char to target_len
    """
    added_A = trimmed = padded = 0

    with gzip.open(input_file, "rt") as in_fh, gzip.open(output_file, "wt") as out_fh:
        for record in SeqIO.parse(in_fh, "fastq"):
            seq = str(record.seq)

            if seq.startswith(add_A_if_prefix):
                seq = "A" + seq
                added_A += 1
            elif seq.startswith(drop_first_if_prefix):
                seq = seq[1:]

            if len(seq) > target_len:
                seq = seq[:target_len]
                trimmed += 1
            elif len(seq) < target_len:
                seq = seq + (pad_char * (target_len - len(seq)))
                padded += 1

            record.seq = Seq(seq)
            SeqIO.write(record, out_fh, "fastq")

    logging.info("Preprocess summary for %s: added_A=%d trimmed=%d padded=%d",
                 input_file.name, added_A, trimmed, padded)
    return added_A, trimmed, padded


def cmd_preprocess(input_folder: Path, output_folder: Path, *, target_len: int, workers: int) -> None:
    output_folder.mkdir(parents=True, exist_ok=True)
    files = [f for f in os.listdir(input_folder) if f.endswith("consensus.fastq.gz")]

    # simple parallelization
    from concurrent.futures import ThreadPoolExecutor, as_completed
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        futs = []
        for fn in files:
            inp = input_folder / fn
            out = output_folder / fn  # keep same name in preprocessed dir
            futs.append(ex.submit(preprocess_reads, inp, out, target_len=target_len))
        for fut in as_completed(futs):
            fut.result()


# -------------------------
# Filter by reference (AA mismatch threshold + optional GG→AG exceptions)
# -------------------------

REFS: Dict[str, str] = {
    "WT": "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATCTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGGAAT",  # noqa: E501
    "Y40X": "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATCTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGGAAT",  # noqa: E501
    "LAI": "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGGAAT",  # noqa: E501
    "1203": "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTGGACAGGATGAAGATTAAAACATGGAAAAGTTTAGTAAAGCATCATATGTATGTTTCAAAGAAGGCTAGGAGATGGTTTTATAGACATCACTATGAAAGCACTCATCCAAAAATAAGTTCAGAAGTACACATCCCACTAGAGAAGGGTGAATTGGTAGTAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAGGGAGATATAGCACACAAGTAGACCCTGACCTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGAATCTGCTATAAGGAAT",  # noqa: E501
}


def translate(nuc: str) -> str:
    return str(Seq(nuc).translate())


def read_ggtoag_csv(path: Path) -> Dict[int, List[str]]:
    """
    Load GG→AG variants CSV with columns: site, variant (codon). Returns {site: [codon, ...]}.
    """
    df = pd.read_csv(path)
    if "site" not in df.columns or "variant" not in df.columns:
        raise ValueError("GG→AG CSV must have columns: site, variant")
    d: Dict[int, List[str]] = df.groupby("site")["variant"].apply(list).to_dict()
    # ensure int keys
    return {int(k): v for k, v in d.items()}


def within_mismatch_threshold(
    merged_seq: str,
    ref_seq: str,
    *,
    max_aa_mismatches: int,
    ggtoag_allow: Dict[int, List[str]] | None = None,
) -> bool:
    """
    Compare translated AA sequences; allow up to max_aa_mismatches differences.
    If ggtoag_allow is provided, substitutions that produce a codon in ggtoag_allow[site]
    are not counted as mismatches (site is 1-based codon index).
    """
    merged_prot = translate(merged_seq)
    ref_prot = translate(ref_seq)
    mismatches = 0

    for i, (a, b) in enumerate(zip(merged_prot, ref_prot[: len(merged_prot)]), start=1):
        if a == b:
            continue
        if ggtoag_allow:
            codon_start = (i - 1) * 3
            codon = merged_seq[codon_start: codon_start + 3]
            wt_codon = ref_seq[codon_start: codon_start + 3]
            allowed = ggtoag_allow.get(i, [])
            if codon in allowed or codon == wt_codon:
                continue
        mismatches += 1
        if mismatches > max_aa_mismatches:
            return False
    return True


def filter_file(
    merged_file: Path,
    output_file: Path,
    *,
    ref_label: str,
    max_aa_mismatches: int,
    ggtoag_csv: Path | None,
) -> Tuple[int, int]:
    ref_seq = REFS[ref_label]
    ggtoag_allow = read_ggtoag_csv(ggtoag_csv) if ggtoag_csv else None

    kept = removed = 0
    with gzip.open(merged_file, "rt") as in_fh, gzip.open(output_file, "wt") as out_fh:
        for rec in SeqIO.parse(in_fh, "fastq"):
            if within_mismatch_threshold(str(rec.seq), ref_seq, max_aa_mismatches=max_aa_mismatches, ggtoag_allow=ggtoag_allow):
                SeqIO.write(rec, out_fh, "fastq")
                kept += 1
            else:
                removed += 1

    logging.info("%s: kept=%d removed=%d", output_file.name, kept, removed)
    return kept, removed


def infer_ref_from_name(name: str) -> str | None:
    for k in REFS:
        if k in name:
            return k
    return None


def cmd_filter(input_folder: Path, output_folder: Path, *, pattern: str, max_workers: int, max_aa_mismatches: int, ggtoag_csv: Path | None) -> None:
    output_folder.mkdir(parents=True, exist_ok=True)
    from concurrent.futures import ProcessPoolExecutor, as_completed

    files = [f for f in os.listdir(input_folder) if f.endswith(".fastq.gz") and pattern in f]
    futs = []
    with ProcessPoolExecutor(max_workers=max(1, max_workers)) as ex:
        for fn in files:
            ref_label = infer_ref_from_name(fn)
            if not ref_label:
                logging.warning("Skipping %s: could not infer reference from name", fn)
                continue
            inp = input_folder / fn
            out = output_folder / fn.replace("_merged_trimmed", "_filtered")
            futs.append(ex.submit(filter_file, inp, out, ref_label=ref_label, max_aa_mismatches=max_aa_mismatches, ggtoag_csv=ggtoag_csv))
        for fut in as_completed(futs):
            fut.result()


# -------------------------
# Codon counting
# -------------------------

EXHAUSTIVE_CODONS = [
    "AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT",
    "ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT",
    "CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT",
    "GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
    "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT",
    "TTA","TTC","TTG","TTT"
]


def get_reference_for_basename(base: str) -> str | None:
    return infer_ref_from_name(base)


def count_codons_in_fastq(fastq_path: Path, wildtype: str) -> pd.DataFrame:
    ncodons = len(wildtype) // 3
    counts: Dict[int, Dict[str, int]] = {i: {codon: 0 for codon in EXHAUSTIVE_CODONS} for i in range(ncodons)}

    with open(fastq_path, "rt") as fh:
        for record in SeqIO.parse(fh, "fastq"):
            seq = str(record.seq[: len(wildtype)])
            for i in range(0, len(wildtype) - 2, 3):
                codon = seq[i : i + 3]
                if codon in EXHAUSTIVE_CODONS:
                    counts[i // 3][codon] += 1

    df = pd.DataFrame(counts).T.fillna(0).astype(int)
    df.index += 1
    df.index.name = "site"
    wt_trios = [wildtype[i : i + 3] for i in range(0, len(wildtype), 3)]
    df["wildtype"] = wt_trios
    cols = ["wildtype"] + EXHAUSTIVE_CODONS
    return df[cols]


def decompress_gz(src: Path, dest: Path) -> None:
    with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def cmd_count_codons(input_folder: Path, pattern: str, output_folder: Path) -> None:
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    import glob
    files = []
    for p in [pattern]:
        files.extend(glob.glob(str(input_folder / p)))

    for gz in files:
        gz_path = Path(gz)
        base = gz_path.stem.replace(".fastq", "")  # handle .fastq.gz -> .fastq
        out_csv = Path(output_folder) / f"{base}_codoncounts.csv"
        if out_csv.exists():
            logging.info("Skipping (exists): %s", out_csv.name)
            continue

        # Decompress to temp .fastq in the same folder
        unzipped = gz_path.with_suffix("")  # drop .gz -> .fastq
        if not unzipped.exists():
            logging.info("Decompressing %s -> %s", gz_path.name, unzipped.name)
            decompress_gz(gz_path, unzipped)

        ref_label = get_reference_for_basename(base)
        if not ref_label:
            logging.warning("Unknown reference for %s; skipping.", base)
            continue
        df = count_codons_in_fastq(unzipped, REFS[ref_label])
        df.to_csv(out_csv)
        logging.info("Wrote codon counts: %s", out_csv.name)


# -------------------------
# CLI
# -------------------------

def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter and analyze consensus FASTQ reads.")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    sub = p.add_subparsers(dest="cmd", required=True)

    # preprocess
    pp = sub.add_parser("preprocess", help="Normalize read length and fix specific 5' anomalies.")
    pp.add_argument("input", type=Path, help="Folder containing *_consensus.fastq.gz files")
    pp.add_argument("out", type=Path, help="Output folder for preprocessed fastq.gz files")
    pp.add_argument("--target-len", type=int, default=366)
    pp.add_argument("--workers", type=int, default=16)

    # filter
    pf = sub.add_parser("filter", help="Filter reads by AA mismatch threshold vs. reference.")
    pf.add_argument("input", type=Path, help="Folder with input .fastq.gz files (e.g., preprocessed)")
    pf.add_argument("out", type=Path, help="Output folder for *_filtered.fastq.gz")
    pf.add_argument("--pattern", default="_merged_trimmed", help="Substring to select files (default: _merged_trimmed)")
    pf.add_argument("--max-aa-mismatches", type=int, default=4, help="Max allowed AA mismatches (default: 4)")
    pf.add_argument("--ggtoag-csv", type=Path, default=None, help="CSV of allowed GG→AG codons per site (columns: site, variant)")

    # count-codons
    pc = sub.add_parser("count-codons", help="Generate codon counts CSVs from filtered FASTQ files.")
    pc.add_argument("input", type=Path, help="Folder containing filtered .fastq.gz files")
    pc.add_argument("--pattern", default="*_filtered*.fastq.gz", help="Glob pattern to select files (default: *_filtered*.fastq.gz)")
    pc.add_argument("--out", type=Path, default=None, help="Output folder for CSVs (default: same as input)")

    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    init_logger(args.log_level)

    if args.cmd == "preprocess":
        cmd_preprocess(args.input, args.out, target_len=args.target_len, workers=args.workers)
    elif args.cmd == "filter":
        cmd_filter(args.input, args.out, pattern=args.pattern, max_workers=32, max_aa_mismatches=args.max_aa_mismatches, ggtoag_csv=args.ggtoag_csv)
    elif args.cmd == "count-codons":
        out = args.out if args.out is not None else args.input
        cmd_count_codons(args.input, pattern=args.pattern, output_folder=out)
    else:
        raise SystemExit("Unknown command")

    return 0


if __name__ == "__main__":
    sys.exit(main())
