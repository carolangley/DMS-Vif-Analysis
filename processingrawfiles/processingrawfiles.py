#!/usr/bin/env python3
"""
processingrawfiles.py

Pipeline utilities for FASTQ processing:
1) Merge paired-end reads with BBMerge (and trim/length-filter merged reads).
2) Build per-UMI consensus reads from merged-and-trimmed FASTQ files.

Requirements:
  - BBMerge (bbmerge.sh available on PATH)
  - Python 3.8+
  - Biopython
  - tqdm (optional; used for progress reporting)
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import logging
import os
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

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
# Merge + Trim
# -------------------------

def trim_and_filter_sequences(input_file: Path, output_file: Path, *, min_length: int = 384, max_length: int = 384) -> None:
    """Trim sequences longer than max_length and discard sequences shorter than min_length.

    Both input and output are gzipped FASTQ files.
    """
    with gzip.open(input_file, "rt") as in_fh, gzip.open(output_file, "wt") as out_fh:
        for record in SeqIO.parse(in_fh, "fastq"):
            if len(record.seq) >= min_length:
                trimmed_seq = record.seq[:max_length]
                trimmed_record = SeqRecord(trimmed_seq, id=record.id, description=record.description)
                trimmed_record.letter_annotations[
                    "phred_quality"
                ] = record.letter_annotations["phred_quality"][:max_length]
                SeqIO.write(trimmed_record, out_fh, "fastq")


def run_bbmerge(r1_file: Path, r2_file: Path, output_file: Path) -> None:
    """Run BBMerge to merge paired reads.

Requires 'bbmerge.sh' on PATH.
"""
    command = [
        "bbmerge.sh",
        f"in1={r1_file}",
        f"in2={r2_file}",
        f"out={output_file}",
        f"outu1={str(output_file).replace('_merged.fastq.gz', '_unmerged_R1.fastq.gz')}",
        f"outu2={str(output_file).replace('_merged.fastq.gz', '_unmerged_R2.fastq.gz')}",
    ]
    logging.debug("Running: %s", " ".join(command))
    subprocess.run(command, check=True)


def process_pair(r1_file: Path, r2_file: Path, out_dir: Path, *, min_length: int, max_length: int) -> None:
    """Merge a single R1/R2 pair and produce a trimmed merged file."""
    sample = Path(r1_file).name.split("_")[0]
    merged_path = out_dir / f"{sample}_merged.fastq.gz"
    trimmed_path = out_dir / f"{sample}_merged_trimmed.fastq.gz"

    if trimmed_path.exists():
        logging.info("Trimmed exists, skipping: %s", trimmed_path.name)
        return

    if merged_path.exists():
        logging.info("Merged exists, trimming only: %s", merged_path.name)
        trim_and_filter_sequences(merged_path, trimmed_path, min_length=min_length, max_length=max_length)
    else:
        logging.info("Merging: %s + %s -> %s", r1_file.name, r2_file.name, merged_path.name)
        run_bbmerge(r1_file, r2_file, merged_path)
        trim_and_filter_sequences(merged_path, trimmed_path, min_length=min_length, max_length=max_length)


def cmd_merge_trim(folder_path: Path, output_folder: Path, *, min_length: int, max_length: int, workers: int) -> None:
    output_folder.mkdir(parents=True, exist_ok=True)

    # Find pairs
    pairs: List[Tuple[Path, Path]] = []
    for fn in os.listdir(folder_path):
        if fn.endswith(".fastq.gz") and "_R1" in fn:
            r1 = folder_path / fn
            r2 = folder_path / fn.replace("_R1", "_R2")
            if not r2.exists():
                logging.warning("Missing R2 for %s", r1.name)
                continue
            pairs.append((r1, r2))

    logging.info("Found %d R1/R2 pairs.", len(pairs))
    if not pairs:
        return

    max_workers = min(max(1, workers), len(pairs))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [
            ex.submit(process_pair, r1, r2, output_folder, min_length=min_length, max_length=max_length)
            for (r1, r2) in pairs
        ]
        concurrent.futures.wait(futures)

    # Ensure trimming for any pre-existing merged files
    merged_files = [output_folder / f for f in os.listdir(output_folder) if f.endswith("_merged.fastq.gz")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = []
        for merged in merged_files:
            trimmed = Path(str(merged).replace("_merged.fastq.gz", "_merged_trimmed.fastq.gz"))
            if not trimmed.exists():
                futs.append(ex.submit(trim_and_filter_sequences, merged, trimmed, min_length=min_length, max_length=max_length))
        concurrent.futures.wait(futs)


# -------------------------
# UMI Consensus
# -------------------------

def consensus_base(bases: List[str], qualities: List[int]) -> str:
    """Pick the base with highest count; tie-break by highest quality base observed."""
    counts = Counter(bases).most_common(2)
    if len(counts) > 1 and counts[0][1] == counts[1][1]:
        # Tie: choose the base with highest quality among A/G/C/T
        for base, q in sorted(zip(bases, qualities), key=lambda x: x[1], reverse=True):
            if base in "AGCT":
                return base
    return counts[0][0]


def gc_filter(read: str, _quality: List[int], *, min_gc: int = 20, max_gc: int = 80) -> bool:
    gc = 100 * (read.count("G") + read.count("C")) / max(1, len(read))
    return min_gc <= gc <= max_gc


def n_content_filter(read: str, _quality: List[int], *, max_n_percentage: int = 10) -> bool:
    n_pct = 100 * read.count("N") / max(1, len(read))
    return n_pct <= max_n_percentage


def quality_filter(_read: str, quality: List[int], *, min_quality: int = 20) -> bool:
    avg_q = sum(quality) / max(1, len(quality))
    return avg_q >= min_quality


def filter_reads(reads: List[str], qualities: List[List[int]], filters) -> Tuple[List[str], List[List[int]]]:
    kept_r, kept_q = [], []
    for r, q in zip(reads, qualities):
        if all(f(r, q) for f in filters):
            kept_r.append(r)
            kept_q.append(q)
    return kept_r, kept_q


def group_reads_by_umi(file_path: Path, *, umi_length: int) -> dict:
    """Group reads by combined 5' and 3' UMI of given length.

Expects gzipped FASTQ.
"""
    groups = defaultdict(lambda: {"sequences": [], "qualities": []})
    with gzip.open(file_path, "rt") as fh:
        for record in SeqIO.parse(fh, "fastq"):
            seq = str(record.seq)
            umi5 = seq[:umi_length]
            umi3 = seq[-umi_length:]
            combined = f"{umi5}_{umi3}"
            trimmed_seq = seq[umi_length:-umi_length]
            trimmed_qual = record.letter_annotations["phred_quality"][umi_length:-umi_length]
            groups[combined]["sequences"].append(trimmed_seq)
            groups[combined]["qualities"].append(trimmed_qual)
    return groups


def build_consensus_for_group(umi: str, data: dict, filters) -> SeqRecord | None:
    seqs = data["sequences"]
    quals = data["qualities"]
    seqs, quals = filter_reads(seqs, quals, filters)
    if not seqs:
        return None

    cons_bases: List[str] = []
    cons_quals: List[int] = []
    read_len = len(seqs[0])
    for i in range(read_len):
        col_bases = [s[i] for s in seqs]
        col_quals = [q[i] for q in quals]
        cons_bases.append(consensus_base(col_bases, col_quals))
        cons_quals.append(max(col_quals))

    return SeqRecord(Seq("".join(cons_bases)), id=umi, description="consensus sequence", letter_annotations={"phred_quality": cons_quals})


def process_merged_file(merged_file: Path, output_folder: Path, *, umi_length: int, min_gc: int, max_gc: int, max_n: int, min_q: int) -> None:
    out_path = output_folder / (merged_file.name.replace(".fastq.gz", "_consensus.fastq.gz"))
    if out_path.exists():
        logging.info("Consensus exists, skipping: %s", out_path.name)
        return

    groups = group_reads_by_umi(merged_file, umi_length=umi_length)
    filters = [
        lambda r, q: gc_filter(r, q, min_gc=min_gc, max_gc=max_gc),
        lambda r, q: n_content_filter(r, q, max_n_percentage=max_n),
        lambda r, q: quality_filter(r, q, min_quality=min_q),
    ]

    consensus_reads: List[SeqRecord] = []
    for umi, data in groups.items():
        rec = build_consensus_for_group(umi, data, filters)
        if rec is not None:
            consensus_reads.append(rec)

    if consensus_reads:
        with gzip.open(out_path, "wt") as fh:
            SeqIO.write(consensus_reads, fh, "fastq")
        logging.info("Wrote consensus: %s (%d reads)", out_path.name, len(consensus_reads))
    else:
        logging.warning("No consensus reads generated for %s", merged_file.name)


def cmd_consensus(input_folder: Path, output_folder: Path, *, umi_length: int, workers: int, min_gc: int, max_gc: int, max_n: int, min_q: int) -> None:
    output_folder.mkdir(parents=True, exist_ok=True)
    merged = [input_folder / f for f in os.listdir(input_folder) if f.endswith("trimmed.fastq.gz")]
    logging.info("Found %d merged/trimmed files.", len(merged))

    max_workers = min(max(1, workers), len(merged) if merged else 1)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [
            ex.submit(
                process_merged_file,
                mf,
                output_folder,
                umi_length=umi_length,
                min_gc=min_gc,
                max_gc=max_gc,
                max_n=max_n,
                min_q=min_q,
            )
            for mf in merged
        ]
        concurrent.futures.wait(futures)


# -------------------------
# CLI
# -------------------------

def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="FASTQ processing pipeline (merge/trim and UMI consensus).")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")


    sub = p.add_subparsers(dest="cmd", required=True)

    # merge-trim
    pm = sub.add_parser("merge-trim", help="Merge paired FASTQs with BBMerge and trim to fixed length.")
    pm.add_argument("folder", type=Path, help="Folder containing paired R1/R2 .fastq.gz files (expects *_R1*.fastq.gz and *_R2*.fastq.gz)")
    pm.add_argument("out", type=Path, help="Output folder for merged/trimmed files")
    pm.add_argument("--min-length", type=int, default=384, help="Minimum read length to keep (default: 384)")
    pm.add_argument("--max-length", type=int, default=384, help="Trim to this length (default: 384)")
    pm.add_argument("--workers", type=int, default=8, help="Parallel workers (default: 8)")

    # consensus
    pc = sub.add_parser("consensus", help="Build UMI consensus reads from merged/trimmed FASTQ files.")
    pc.add_argument("input", type=Path, help="Folder with *_merged_trimmed.fastq.gz files")
    pc.add_argument("out", type=Path, help="Output folder for *_consensus.fastq.gz files")
    pc.add_argument("--umi-length", type=int, default=9, help="UMI length at each end (default: 9)")
    pc.add_argument("--workers", type=int, default=8, help="Parallel workers (default: 8)")
    pc.add_argument("--min-gc", type=int, default=20, help="Minimum GC%% (default: 20)")
    pc.add_argument("--max-gc", type=int, default=80, help="Maximum GC%% (default: 80)")
    pc.add_argument("--max-n", type=int, default=10, help="Maximum N%% (default: 10)")
    pc.add_argument("--min-q", type=int, default=20, help="Minimum average Phred quality (default: 20)")

    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    init_logger(args.log_level)

    if args.cmd == "merge-trim":
        cmd_merge_trim(args.folder, args.out, min_length=args.min_length, max_length=args.max_length, workers=args.workers)
    elif args.cmd == "consensus":
        cmd_consensus(args.input, args.out, umi_length=args.umi_length, workers=args.workers, min_gc=args.min_gc, max_gc=args.max_gc, max_n=args.max_n, min_q=args.min_q)
    else:
        raise SystemExit("Unknown command")  # shouldn't happen

    return 0


if __name__ == "__main__":
    sys.exit(main())
