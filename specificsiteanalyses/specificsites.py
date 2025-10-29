#!/usr/bin/env python3
"""
specificsites.py

"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LIGHTPURPLE = (0.87, 0.63, 0.87)
LIGHTGREEN = (0.56, 0.93, 0.56)

AMINO_ACID_ORDER = ['A','V','I','L','M','F','W','P','G','S','T','C','Y','N','Q','D','E','K','R','H']

def _ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)

def _load_enrich(enrich_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(enrich_csv)
    # Defensive: drop obviously missing columns early
    need = {"site","wt_aa","aa","sample_post","log_enrichment_ratio"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{enrich_csv} missing required columns: {missing}")
    return df

def _mean_std_by_aa(filtered: pd.DataFrame) -> pd.DataFrame:
    g = (filtered.groupby(['aa','group'])['log_enrichment_ratio']
                 .agg(['mean','std']).reset_index())
    # bespoke order, keep only AAs present
    g['sort_key'] = g['aa'].apply(lambda x: AMINO_ACID_ORDER.index(x) if x in AMINO_ACID_ORDER else -1)
    g = g[g['sort_key'] >= 0].sort_values(['sort_key','group'])
    return g

def _plot_one_site(
    site: int,
    df: pd.DataFrame,
    samples_1203: List[str],
    samples_LAI: List[str],
    out_dir: Path,
    ylim: Optional[List[float]] = None,
    include_text: bool = False,
    exclude_star: bool = True,
) -> None:
    # Filter rows for site, drop WT matches; optionally drop '*'
    mask = (df['site'] == site) & (df['wt_aa'] != df['aa']) & (df['sample_post'].isin(samples_1203 + samples_LAI))
    if exclude_star and 'aa' in df.columns:
        mask &= (df['aa'] != '*')
    sub = df.loc[mask].copy()
    if sub.empty:
        return

    sub.loc[:, 'group'] = np.where(sub['sample_post'].isin(samples_1203), '1203', 'LAI')
    stats_df = _mean_std_by_aa(sub)
    if stats_df.empty:
        return
    unique_aas = list(stats_df['aa'].unique())
    x_positions = {aa: i for i, aa in enumerate(unique_aas)}

    # Plot
    plt.figure(figsize=(10, 4))
    palette = {'1203': LIGHTPURPLE, 'LAI': LIGHTGREEN}
    markers = {'1203': 'o', 'LAI': 's'}
    zorder = {'1203': 3, 'LAI': 2}

    # error bars first
    for grp, grp_df in stats_df.groupby('group'):
        xs = [x_positions[a] for a in grp_df['aa']]
        plt.errorbar(xs, grp_df['mean'], yerr=grp_df['std'], fmt='none', capsize=10, elinewidth=2,
                     color=palette[grp], zorder=1)

    # points + optional text labels
    for grp, grp_df in stats_df.groupby('group'):
        xs = [x_positions[a] for a in grp_df['aa']]
        plt.scatter(xs, grp_df['mean'], color=palette[grp], marker=markers[grp], s=200,
                    edgecolor='black', label=grp, zorder=zorder[grp])
        if include_text:
            for x, m in zip(xs, grp_df['mean']):
                plt.text(x, m, f"{m:.2f}", fontsize=10, ha='right', va='bottom', zorder=4)

    plt.title(f"Site {site} in 1203 and LAI replicates")
    plt.xlabel('Amino Acid')
    plt.ylabel('Log Enrichment Ratio')
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    else:
        plt.ylim(-8, 8)
    plt.xticks(ticks=range(len(unique_aas)), labels=unique_aas, rotation=45, ha='right')
    plt.axhline(y=0, color='r', linestyle='--', linewidth=1)
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    plt.legend(title='Group', loc='upper right')
    _ensure_outdir(out_dir)
    fname = f"Site_{site}_LAI_vs_1203.pdf"
    plt.tight_layout()
    plt.savefig(out_dir / fname, bbox_inches='tight')
    plt.close()

def _plot_lai_only(
    site: int,
    df: pd.DataFrame,
    samples_LAI: List[str],
    out_dir: Path,
    ylim: Optional[List[float]] = None,
) -> None:
    mask = (df['site'] == site) & (df['wt_aa'] != df['aa']) & (df['sample_post'].isin(samples_LAI))
    sub = df.loc[mask].copy()
    if sub.empty:
        return
    sub.loc[:, 'group'] = 'LAI'
    stats_df = _mean_std_by_aa(sub)
    if stats_df.empty:
        return
    unique_aas = list(stats_df['aa'])
    xs = list(range(len(unique_aas)))

    plt.figure(figsize=(8, 4))
    plt.errorbar(xs, stats_df['mean'], yerr=stats_df['std'], fmt='none', capsize=10, elinewidth=2,
                 color=LIGHTGREEN, zorder=1)
    plt.scatter(xs, stats_df['mean'], color=LIGHTGREEN, marker='s', s=200, edgecolor='black',
                label='LAI', zorder=2)

    plt.title(f"Site {site} in LAI replicates")
    plt.xlabel('Amino Acid')
    plt.ylabel('Log Enrichment Ratio')
    plt.axhline(0, color='r', linestyle='--', linewidth=1)
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    else:
        plt.ylim(-6, 6)
    plt.xticks(ticks=xs, labels=unique_aas, rotation=45, ha='right')
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    _ensure_outdir(out_dir)
    plt.savefig(out_dir / f"Site_{site}_LAI_only.pdf", bbox_inches='tight')
    plt.close()

def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-site AA scatterplots from enrich_df.csv")    
    p.add_argument("--enrich-csv", type=Path, default=Path("enrich_df.csv"), help="Path to enrich_df.csv")    
    p.add_argument("--sites", type=int, nargs="+", required=True, help="One or more site indices to plot")    
    p.add_argument("--samples-1203", nargs="+", default=["post1203R1","post1203R2","post1203R3"], help="Sample names for 1203 replicates")    
    p.add_argument("--samples-LAI", nargs="+", default=["postLAIR1","postLAIR2","postLAIR3"], help="Sample names for LAI replicates")    
    p.add_argument("--outdir", type=Path, default=Path("figures/scatterplots"), help="Output folder")    
    p.add_argument("--ylim", type=float, nargs=2, default=None, help="Y-axis limits, e.g. --ylim -8 8")    
    p.add_argument("--exclude-star", action="store_true", help="Exclude stop '*' from plots (default off)")    
    p.add_argument("--labels", action="store_true", help="Annotate points with mean values")    
    p.add_argument("--lai-only", action="store_true", help="Also generate LAI-only plots")    
    return p.parse_args(list(argv) if argv is not None else None)

def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    df = _load_enrich(args.enrich_csv)

    for s in args.sites:
        _plot_one_site(
            site=s,
            df=df,
            samples_1203=args.samples_1203,
            samples_LAI=args.samples_LAI,
            out_dir=args.outdir,
            ylim=args.ylim,
            include_text=args.labels,
            exclude_star=args.exclude_star,
        )
        if args.lai_only:
            _plot_lai_only(
                site=s,
                df=df,
                samples_LAI=args.samples_LAI,
                out_dir=args.outdir,
                ylim=args.ylim,
            )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
