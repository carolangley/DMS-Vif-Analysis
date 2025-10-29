#!/usr/bin/env python3
"""
y40pilotexperiments.py

Compute and visualize site-40 TAG (WT) and non-TAT (Y40X) percentages from
codon-count CSVs, summarizing pre/post selection and input controls.

This refactor converts the original notebook-style script into a robust CLI with:
- explicit arguments (input folder, output paths, site, patterns)
- non-interactive plotting (matplotlib Agg)
- clear logging and error handling
- reproducible outputs

Expected input files:
  consensuseachUMI/filtered/*WT*codoncounts*.csv
  consensuseachUMI/filtered/*Y40X*codoncounts*.csv

By default, the script searches for files matching the above patterns within
the provided folder and aggregates across replicates R1..R3 for pre/post.


Requirements: pandas, numpy, scipy, matplotlib
"""

from __future__ import annotations

import argparse
import glob
import logging
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import ttest_ind

plt.switch_backend("Agg")


def find_files(folder: Path, pattern: str) -> List[Path]:
    return sorted(Path(p) for p in glob.glob(str(folder / pattern)))


def load_site_percent(df: pd.DataFrame, site: int, col: str) -> Optional[float]:
    if "site" not in df.columns or col not in df.columns:
        return None
    row = df.loc[df["site"] == site]
    if row.empty:
        return None
    # total from all codon columns (exclude site + wildtype if present in col 0/1)
    # robustly sum over all non-'site' columns
    codon_cols = [c for c in df.columns if c != "site"]
    total_reads = row[codon_cols].iloc[:, 1:].sum(axis=1).values[0] if len(codon_cols) > 1 else np.nan
    if total_reads == 0 or np.isnan(total_reads):
        return None
    return float(row[col].values[0]) / float(total_reads) * 100.0


def load_site_non_tat_percent(df: pd.DataFrame, site: int) -> Optional[float]:
    if "site" not in df.columns:
        return None
    row = df.loc[df["site"] == site]
    if row.empty:
        return None
    codon_cols = [c for c in df.columns if c != "site"]
    total_reads = row[codon_cols].iloc[:, 1:].sum(axis=1).values[0] if len(codon_cols) > 1 else np.nan
    if total_reads == 0 or np.isnan(total_reads):
        return None
    tat = float(row.get("TAT", pd.Series([0.0])).values[0])
    return float(total_reads - tat) / float(total_reads) * 100.0


def summarize_patterns(files: List[Path], site: int, *, tag: bool) -> Tuple[pd.DataFrame, Optional[float]]:
    """Return (per-file DataFrame with percentages, input_control_percentage)."""
    rows = []
    input_pct: Optional[float] = None
    for fp in files:
        df = pd.read_csv(fp)
        if tag:
            pct = load_site_percent(df, site, col="TAG")
            if "WTinput" in fp.name and pct is not None:
                input_pct = pct
        else:
            pct = load_site_non_tat_percent(df, site)
            if "Y40Xinput" in fp.name and pct is not None:
                input_pct = pct
        if pct is None:
            continue
        rows.append({"file": fp.name, "pct": pct})
    return pd.DataFrame(rows), input_pct


def pick_sample(file_name: str, prefix: str) -> Optional[str]:
    # Extracts '(pre|post)<prefix>R\d+' from a filename
    import re
    m = re.search(rf'((?:pre|post){prefix}R\d+)_', file_name)
    return m.group(1) if m else None


def aggregate_replicates(df: pd.DataFrame, prefix: str) -> Tuple[List[float], List[float]]:
    pre, post = [], []
    for i in range(1, 4):
        pre_pat = f"pre{prefix}R{i}"
        post_pat = f"post{prefix}R{i}"
        pre_rows = df[df["file"].str.contains(pre_pat)]
        post_rows = df[df["file"].str.contains(post_pat)]
        if not pre_rows.empty:
            pre.append(float(pre_rows["pct"].values[0]))
        if not post_rows.empty:
            post.append(float(post_rows["pct"].values[0]))
    return pre, post


def log_change(pre: float, post: float) -> float:
    return np.log((post / 100.0)) - np.log((pre / 100.0))


def relative_percent_change(log_delta: float) -> float:
    return (np.exp(log_delta) - 1.0) * 100.0


def build_plot(wtinput: Optional[float], pre_wt: List[float], post_wt: List[float],
               y40xinput: Optional[float], pre_y40x: List[float], post_y40x: List[float],
               out_pdf: Path) -> None:
    # Compute summary stats
    def avg_sd(vals: List[float]) -> Tuple[float, float]:
        return (float(np.mean(vals)), float(np.std(vals))) if vals else (np.nan, np.nan)

    avg_pre_wt, sd_pre_wt = avg_sd(pre_wt)
    avg_post_wt, sd_post_wt = avg_sd(post_wt)
    avg_pre_y, sd_pre_y = avg_sd(pre_y40x)
    avg_post_y, sd_post_y = avg_sd(post_y40x)

    fig, axes = plt.subplots(1, 2, figsize=(6, 2))

    # Panel 1: TAG%
    axes[0].bar(
        ["WTinput", "Pre-selection", "Post-selection"],
        [wtinput if wtinput is not None else 0.0, avg_pre_wt, avg_post_wt],
        yerr=[0.0, sd_pre_wt, sd_post_wt],
        capsize=5,
        color="lightgrey",
        edgecolor="black",
    )
    axes[0].set_ylabel("TAG Percentage (%)")
    axes[0].set_title("TAG Percentage")
    axes[0].set_ylim(0, 120)

    # Panel 2: Non-TAT%
    axes[1].bar(
        ["Y40Xinput", "Pre-selection", "Post-selection"],
        [y40xinput if y40xinput is not None else 0.0, avg_pre_y, avg_post_y],
        yerr=[0.0, sd_pre_y, sd_post_y],
        capsize=5,
        color="lightgrey",
        edgecolor="black",
    )
    axes[1].set_ylabel("Non-TAT Percentage (%)")
    axes[1].set_title("Non-TAT Percentage")
    axes[1].set_ylim(0, 120)

    plt.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize site-40 TAG and non-TAT percentages and export a figure + CSV.")
    p.add_argument("--folder", type=Path, default=Path("consensuseachUMI/filtered"), help="Folder with *codoncounts*.csv files")
    p.add_argument("--site", type=int, default=40, help="Codon site index to analyze (default 40)")
    p.add_argument("--out-pdf", type=Path, default=Path("figures/Fig1B.pdf"), help="Output PDF path (default figures/Fig1B.pdf)")
    p.add_argument("--out-csv", type=Path, default=Path("figures/Fig1B_summary.csv"), help="Output CSV with summary stats")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(message)s")

    # Locate files
    csv_files_wt = find_files(args.folder, "*WT*codoncounts*.csv")
    csv_files_y40x = find_files(args.folder, "*Y40X*codoncounts*.csv")

    logging.info("Found %d WT and %d Y40X files", len(csv_files_wt), len(csv_files_y40x))

    # Summaries
    agg_wt, wtinput = summarize_patterns(csv_files_wt, args.site, tag=True)
    agg_y, y40xinput = summarize_patterns(csv_files_y40x, args.site, tag=False)

    pre_wt, post_wt = aggregate_replicates(agg_wt, prefix="WT")
    pre_y, post_y = aggregate_replicates(agg_y, prefix="Y40X")

    # Compute logs & relative changes
    def safe_log_change(pre_list: List[float], post_list: List[float]) -> float:
        import numpy as _np
        if not pre_list or not post_list:
            return _np.nan
        return float(_np.log((_np.mean(post_list) / 100.0)) - _np.log((_np.mean(pre_list) / 100.0)))

    def safe_rel_change(ld: float) -> float:
        import numpy as _np
        return float((_np.exp(ld) - 1.0) * 100.0) if not _np.isnan(ld) else _np.nan

    lc_tag_pre_post = safe_log_change(pre_wt, post_wt)
    lc_non_tat_pre_post = safe_log_change(pre_y, post_y)

    lc_tag_post_input = (float(np.log((np.mean(post_wt)/100.0)) - np.log((wtinput/100.0)))
                         if wtinput is not None and post_wt else float("nan"))
    lc_non_tat_post_input = (float(np.log((np.mean(post_y)/100.0)) - np.log((y40xinput/100.0)))
                             if y40xinput is not None and post_y else float("nan"))

    rel_tag_pre_post = safe_rel_change(lc_tag_pre_post)
    rel_tag_post_input = safe_rel_change(lc_tag_post_input)
    rel_non_tat_pre_post = safe_rel_change(lc_non_tat_pre_post)
    rel_non_tat_post_input = safe_rel_change(lc_non_tat_post_input)

    # t-tests
    def safe_ttest(a: List[float], b_val: Optional[float]) -> float:
        import numpy as _np
        if not a or b_val is None:
            return float("nan")
        return float(ttest_ind(a, [b_val]*len(a))[1])

    p_pre_input_tag = safe_ttest(pre_wt, wtinput)
    p_post_input_tag = safe_ttest(post_wt, wtinput)
    p_pre_post_tag = float(ttest_ind(pre_wt, post_wt)[1]) if pre_wt and post_wt else float("nan")

    p_pre_input_non_tat = safe_ttest(pre_y, y40xinput)
    p_post_input_non_tat = safe_ttest(post_y, y40xinput)
    p_pre_post_non_tat = float(ttest_ind(pre_y, post_y)[1]) if pre_y and post_y else float("nan")

    # Plot
    build_plot(wtinput, pre_wt, post_wt, y40xinput, pre_y, post_y, args.out_pdf)

    # Write CSV summary
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame({
        "metric": [
            "avg_pre_tag", "sd_pre_tag", "avg_post_tag", "sd_post_tag",
            "avg_pre_non_tat", "sd_pre_non_tat", "avg_post_non_tat", "sd_post_non_tat",
            "wtinput_tag", "y40xinput_non_tat",
            "log_change_tag_pre_post", "rel_change_tag_pre_post",
            "log_change_tag_post_input", "rel_change_tag_post_input",
            "log_change_non_tat_pre_post", "rel_change_non_tat_pre_post",
            "log_change_non_tat_post_input", "rel_change_non_tat_post_input",
            "p_pre_input_tag", "p_post_input_tag", "p_pre_post_tag",
            "p_pre_input_non_tat", "p_post_input_non_tat", "p_pre_post_non_tat",
        ],
        "value": [
            float(np.mean(pre_wt)) if pre_wt else float("nan"),
            float(np.std(pre_wt)) if pre_wt else float("nan"),
            float(np.mean(post_wt)) if post_wt else float("nan"),
            float(np.std(post_wt)) if post_wt else float("nan"),
            float(np.mean(pre_y)) if pre_y else float("nan"),
            float(np.std(pre_y)) if pre_y else float("nan"),
            float(np.mean(post_y)) if post_y else float("nan"),
            float(np.std(post_y)) if post_y else float("nan"),
            float(wtinput) if wtinput is not None else float("nan"),
            float(y40xinput) if y40xinput is not None else float("nan"),
            lc_tag_pre_post, rel_tag_pre_post,
            lc_tag_post_input, rel_tag_post_input,
            lc_non_tat_pre_post, rel_non_tat_pre_post,
            lc_non_tat_post_input, rel_non_tat_post_input,
            p_pre_input_tag, p_post_input_tag, p_pre_post_tag,
            p_pre_input_non_tat, p_post_input_non_tat, p_pre_post_non_tat,
        ]
    })
    summary.to_csv(args.out_csv, index=False)

    logging.info("Wrote figure to %s and summary to %s", args.out_pdf, args.out_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
