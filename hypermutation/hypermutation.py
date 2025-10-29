#!/usr/bin/env python3
"""
hypermutation.py

Extracts reads from FASTQ files that carry specific codon mutations defined
in an enrichment CSV. For each (site, codon, wildtype, aa_mutation) in the CSV
matching a given sample, the script writes matching reads to a per-mutation
FASTA file under the output directory.

Requirements:
    biopython
    pandas
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass(frozen=True)
class MutationRow:
    site: int
    codon: str
    wildtype: str
    aa_mutation: str

    @classmethod
    def from_series(cls, s: pd.Series) -> """Construct from a pandas Series with required fields.""" :
        return cls(
            site=int(s["site"]),
            codon=str(s["codon"]),
            wildtype=str(s["wildtype"]),
            aa_mutation=str(s.get("aa_mutation", "")),
        )


def load_enrich_df(path: Path) -> pd.DataFrame:
    """Load and validate the enrichment CSV.

    Expected columns: site, codon, wildtype, aa_mutation, sample_post
    """
    df = pd.read_csv(path)
    required = {"site", "codon", "wildtype", "sample_post"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"enrich CSV is missing required columns: {sorted(missing)}")
    # Optional: normalize types
    df["site"] = pd.to_numeric(df["site"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["site", "codon", "wildtype", "sample_post"]).copy()
    return df


def open_maybe_gzip(path: Path, mode: str = "rt"):
    """Open plain text or gzipped files based on extension."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[no-any-return]
    return open(path, mode, encoding="utf-8")  # type: ignore[no-any-return]


def sanitize_sample_id(name: str) -> str:
    """Derive a base sample identifier: strip trailing '_filtered' if present."""
    return name[:-9] if name.endswith("_filtered") else name


def build_output_name(out_dir: Path, sample_id: str, m: MutationRow) -> Path:
    """Construct the per-mutation FASTA output path."""
    fname = f"{sample_id}_{m.codon}_{m.wildtype}_{m.site}_{m.aa_mutation}.fasta"
    return out_dir / fname


def iter_mutations_for_sample(enrich_df: pd.DataFrame, sample_id_base: str) -> List[MutationRow]:
    """Filter enrichment rows for this sample and convert to MutationRow objects."""
    mask = enrich_df["sample_post"].astype(str).str.contains(sample_id_base, na=False)
    rows = [MutationRow.from_series(s) for _, s in enrich_df.loc[mask].iterrows()]
    return rows


def extract_codon(seq: str, site: int) -> str | None:
    """Return the codon (triplet) at 1-based codon site, or None if out of range."""
    start = (site - 1) * 3
    end = start + 3
    if start < 0 or end > len(seq):
        return None
    return seq[start:end]


def process_fastq(
    fastq_path: Path,
    out_dir: Path,
    sample_id: str,
    sample_mutations: List[MutationRow],
    *, buffer_size: int = 10_000,
    site_min: int = 12,
    site_max: int = 115,
) -> Tuple[int, int]:
    """Scan a FASTQ and write reads carrying specified mutant codons to FASTA files.

    Returns:
        (total_records, matched_records)
    """
    t0 = time.time()
    total = 0
    matched = 0
    out_handles: Dict[Tuple[str, int], object] = {}
    buffers: Dict[Tuple[str, int], List[SeqRecord]] = defaultdict(list)

    # Limit to sites of interest and where codon != wildtype
    muts = [m for m in sample_mutations if site_min <= m.site <= site_max and m.codon != m.wildtype]

    def get_handle(m: MutationRow):
        key = (m.codon, m.site)
        if key not in out_handles:
            out_path = build_output_name(out_dir, sample_id, m)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_handles[key] = open(out_path, "a", encoding="utf-8")
        return out_handles[key]

    try:
        with open_maybe_gzip(fastq_path, "rt") as fh:
            for idx, rec in enumerate(SeqIO.parse(fh, "fastq")):
                total += 1
                seq = str(rec.seq)

                if idx and idx % 10_000 == 0:
                    logging.info("%s: processed %d reads (%.1fs)", fastq_path.name, idx, time.time() - t0)

                for m in muts:
                    codon = extract_codon(seq, m.site)
                    if codon is None:
                        continue
                    if codon == m.codon:
                        buffers[(m.codon, m.site)].append(rec)
                        if len(buffers[(m.codon, m.site)]) >= buffer_size:
                            SeqIO.write(buffers[(m.codon, m.site)], get_handle(m), "fasta")
                            buffers[(m.codon, m.site)].clear()
                        matched += 1
    finally:
        # Flush any remaining buffers and close files
        for key, buf in buffers.items():
            if buf:
                # find a representative MutationRow to build name/handle
                rep = next(m for m in muts if (m.codon, m.site) == key)
                SeqIO.write(buf, get_handle(rep), "fasta")
        for h in out_handles.values():
            try:
                h.close()
            except Exception:
                pass

    logging.info("%s: total=%d matched=%d (%.1fs)", fastq_path.name, total, matched, time.time() - t0)
    return total, matched


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract hypermutation-matching reads to FASTA files.")
    p.add_argument("--enrich-csv", required=True, type=Path, help="CSV with columns: site,codon,wildtype,aa_mutation,sample_post")
    p.add_argument("--out-dir", required=True, type=Path, help="Output directory for FASTA files")
    p.add_argument("fastq", nargs="+", type=Path, help="Input FASTQ files (.gz ok)")
    p.add_argument("--buffer-size", type=int, default=10_000, help="Number of reads to buffer per (codon,site) before writing")
    p.add_argument("--site-min", type=int, default=12, help="Minimum codon site to evaluate (inclusive)")
    p.add_argument("--site-max", type=int, default=115, help="Maximum codon site to evaluate (inclusive)")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(message)s")

    enrich_df = load_enrich_df(args.enrich_csv)

    totals = 0
    matches = 0
    for fq in args.fastq:
        sample_id = sanitize_sample_id(fq.stem)
        sample_rows = iter_mutations_for_sample(enrich_df, sample_id)
        if not sample_rows:
            logging.warning("No enrichment rows matched sample_id=%s; skipping %s", sample_id, fq)
            continue
        t, m = process_fastq(
            fastq_path=fq,
            out_dir=args.out_dir,
            sample_id=sample_id,
            sample_mutations=sample_rows,
            buffer_size=args.buffer_size,
            site_min=args.site_min,
            site_max=args.site_max,
        )
        totals += t
        matches += m

    logging.info("Done. Total reads: %d | Total matches: %d", totals, matches)
    return 0


if __name__ == "__main__":
    sys.exit(main())
