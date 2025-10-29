#!/usr/bin/env python3
"""
sitelevelanalyses.py

"""

from __future__ import annotations

import argparse, os, re, glob, io, math, json, warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from collections import Counter

import numpy as np
import pandas as pd

from scipy.stats import wilcoxon, spearmanr, pearsonr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Optional: grids in replicate-qc; safe if available
try:
    import seaborn as sns  # only used for violin aesthetics; not required at runtime
except Exception:  # pragma: no cover
    sns = None

# -------------------------
# Paths and constants
# -------------------------
RAW_DIR   = Path("consensuseachUMI/filtered")
OUT_DIR   = Path("consensuseachUMI/filteredcodoncounts")
QC_DIR    = Path("qc_outputs")
FIG_DIR   = Path("figures")
SUPP_DIR  = Path("supplementary")

SITE_MIN, SITE_MAX = 12, 115
UNMUT_SITES_FOR_BG = list(range(1, 12))
USE_WT_CENTERING = True
PSEUDOCOUNT = 0.5
WT_MIN_FREQ = 0.005
DEPTH_TH_FOR_CONSISTENCY = 100
MIN_UNIQUE_VARIANTS_FOR_STATS = 5

# Fallback colors
GREEN  = (0.133, 0.533, 0.200)
PURPLE = (0.660, 0.200, 0.460)

# -------------------------
# Helpers
# -------------------------
def _sample_basename(path: str) -> str:
    return os.path.basename(path).replace("_codoncounts.csv", "")

def _is_expected_sample_name(name: str) -> bool:
    return isinstance(name, str) and re.fullmatch(r"(pre|post)(LAI|1203)R\d+", name) is not None

def _gather_files(raw_dir: Path) -> Dict[str, List[str]]:
    paths = sorted(glob.glob(str(raw_dir / "*_codoncounts.csv")))
    pre, post = [], []
    for p in paths:
        name = _sample_basename(p)
        if not _is_expected_sample_name(name):
            continue
        (pre if name.startswith("pre") else post).append(p)
    return {"pre": pre, "post": post}

def _homolog_from_name(name: str) -> str:
    if "LAI" in name: return "LAI"
    if "1203" in name: return "1203"
    raise ValueError(f"Cannot infer homolog from sample name: {name}")

def _load_variants_with_wt(csv_path: Path) -> Dict[int, set]:
    df = pd.read_csv(csv_path)
    required = {"site","variant","wildtype"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path} missing columns: {missing}")
    out: Dict[int, set] = {}
    for site, g in df.groupby("site"):
        s = set(g["variant"].dropna().astype(str))
        wt = g["wildtype"].dropna().astype(str)
        if not wt.empty:
            s.add(wt.iloc[0])
        out[int(site)] = s
    return out

def _compute_mean_non_wt_count(df_bg: pd.DataFrame, codon_cols: list) -> float:
    total, nonzero = 0.0, 0
    for _, row in df_bg.iterrows():
        wt = row["wildtype"]
        for codon in codon_cols:
            if codon == wt:
                continue
            cnt = float(row[codon])
            total += cnt
            if cnt > 0:
                nonzero += 1
    return (total / nonzero) if nonzero > 0 else 0.0

def _melt_and_filter_one(path: str, threshold: float, variants_dict: Dict[int,set]) -> pd.DataFrame:
    df = pd.read_csv(path)
    codon_cols = [c for c in df.columns if c not in ("site","wildtype")]
    df[codon_cols] = df[codon_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df = df[(df["site"] >= SITE_MIN) & (df["site"] <= SITE_MAX)].copy()

    long_df = df.melt(id_vars=["site","wildtype"], value_vars=codon_cols,
                      var_name="codon", value_name="counts")
    long_df["site"] = long_df["site"].astype(int)
    long_df["codon"] = long_df["codon"].astype(str)
    long_df["wildtype"] = long_df["wildtype"].astype(str)
    long_df["counts"] = long_df["counts"].astype(float)

    keep = (
        (long_df["codon"] == long_df["wildtype"]) |
        (long_df.apply(lambda r: r["codon"] in variants_dict.get(int(r["site"]), set()), axis=1)) |
        (long_df["counts"] > float(threshold))
    )
    long_df["counts"] = np.where(keep, long_df["counts"], 0.0)
    return long_df[["site","wildtype","codon","counts"]]

def _long_to_wide(df_long: pd.DataFrame) -> pd.DataFrame:
    wide = df_long.pivot_table(index=["site","wildtype"], columns="codon", values="counts",
                               aggfunc="sum", fill_value=0.0).reset_index().sort_values("site")
    return wide[["site","wildtype"] + [c for c in wide.columns if c not in ("site","wildtype")]]

# -------------------------
# Stage 1: prepare
# -------------------------
def cmd_prepare(args):
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)

    files = _gather_files(RAW_DIR)

    # Background thresholds
    background_thresholds = {}
    def _calc_bg(path: str) -> float:
        df = pd.read_csv(path)
        codon_cols = [c for c in df.columns if c not in ("site","wildtype")]
        df[codon_cols] = df[codon_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        df_bg = df[df["site"].isin(UNMUT_SITES_FOR_BG)].copy()
        return _compute_mean_non_wt_count(df_bg, codon_cols)

    for path in files["pre"] + files["post"]:
        name = _sample_basename(path)
        background_thresholds[name] = 0.0 if name.startswith("post") else _calc_bg(path)

    thr_df = pd.DataFrame({"sample": list(background_thresholds.keys()),
                           "threshold": list(background_thresholds.values())}).sort_values("sample")
    thr_df.to_csv(QC_DIR / "background_thresholds.csv", index=False)

    # Variant dicts
    lai_v = _load_variants_with_wt(Path("LAIvariantproportions_with_wildtype.csv"))
    d1203_v = _load_variants_with_wt(Path("1203variantproportions_with_wildtype.csv"))

    def _write_filtered(file_list: List[str], variants: Dict[int,set]):
        for path in file_list:
            name = _sample_basename(path)
            thr = background_thresholds.get(name, 0.0)
            df_long = _melt_and_filter_one(path, thr, variants)
            df_wide = _long_to_wide(df_long)
            outp = OUT_DIR / f"{name}_codoncounts.csv"
            df_wide.to_csv(outp, index=False)

    # Split by homolog
    pre_files_by_h = {"LAI":[], "1203":[]}
    post_files_by_h = {"LAI":[], "1203":[]}
    for path in files["pre"]:
        pre_files_by_h[_homolog_from_name(_sample_basename(path))].append(path)
    for path in files["post"]:
        post_files_by_h[_homolog_from_name(_sample_basename(path))].append(path)

    # Filter per-rep and merged PRE
    _write_filtered(pre_files_by_h["LAI"],  lai_v)
    _write_filtered(pre_files_by_h["1203"], d1203_v)
    _write_filtered(post_files_by_h["LAI"],  lai_v)
    _write_filtered(post_files_by_h["1203"], d1203_v)

    # Pooled PRE (merged) per homolog
    def _write_pre_merged(file_list, variants, homolog):
        if not file_list: return
        merged_long = pd.concat(
            [_melt_and_filter_one(p, background_thresholds.get(_sample_basename(p), 0.0), variants)
             for p in file_list],
            ignore_index=True
        ).groupby(["site","wildtype","codon"], as_index=False)["counts"].sum()
        df_wide = _long_to_wide(merged_long)
        df_wide.to_csv(OUT_DIR / f"mergedpre{homolog}_codoncounts.csv", index=False)

    _write_pre_merged(pre_files_by_h["LAI"],  lai_v,  "LAI")
    _write_pre_merged(pre_files_by_h["1203"], d1203_v, "1203")

# -------------------------
# Stage 2: enrich
# -------------------------
def _aa_of(codon: str) -> str:
    _c2a = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    c = str(codon).upper().replace("U","T")
    return _c2a.get(c, "X")

def _wide_to_long(path: Path, sample_name: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    codon_cols = [c for c in df.columns if c not in ("site","wildtype")]
    df[codon_cols] = df[codon_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df = df[(df["site"] >= SITE_MIN) & (df["site"] <= SITE_MAX)].copy()

    melted = df.melt(id_vars=["site","wildtype"], value_vars=codon_cols,
                     var_name="codon", value_name="counts")
    melted["counts"] = melted["counts"].astype(float)
    melted["counts_w_p"] = melted["counts"] + PSEUDOCOUNT
    melted["site_sum_w_p"] = melted.groupby("site")["counts_w_p"].transform("sum")
    melted["freqs_w_p"] = melted["counts_w_p"] / melted["site_sum_w_p"]
    melted["is_wt"] = melted["codon"] == melted["wildtype"]

    wt_freq = melted[melted["is_wt"]].set_index("site")["freqs_w_p"]
    melted["wt_freq_w_p"] = melted["site"].map(wt_freq).fillna(0.0)

    melted["wt_aa"] = melted["wildtype"].apply(_aa_of)
    melted["aa"]    = melted["codon"].apply(_aa_of)
    melted["sample"] = sample_name
    return melted[["site","wildtype","codon","counts","counts_w_p","freqs_w_p",
                   "wt_freq_w_p","wt_aa","aa","sample"]]

def cmd_enrich(args):
    pre_LAI  = OUT_DIR / "mergedpreLAI_codoncounts.csv"
    pre_1203 = OUT_DIR / "mergedpre1203_codoncounts.csv"
    post_paths = sorted(glob.glob(str(OUT_DIR / "post*R*_codoncounts.csv")))
    if not pre_LAI.exists() or not pre_1203.exists() or not post_paths:
        raise SystemExit("Required filtered codoncounts missing; run 'prepare' first.")

    rows = [
        _wide_to_long(pre_LAI,  "mergedpreLAI"),
        _wide_to_long(pre_1203, "mergedpre1203"),
    ] + [_wide_to_long(Path(p), _sample_basename(p)) for p in post_paths]
    counts_df = pd.concat(rows, ignore_index=True)
    counts_df["homolog"] = counts_df["sample"].apply(lambda s: "LAI" if "LAI" in s else ("1203" if "1203" in s else "unknown"))
    counts_df["site"] = counts_df["site"].astype(int)
    counts_df.to_csv("counts_df.csv", index=False)

    enrich_parts = []
    for path in post_paths:
        post_name = _sample_basename(path)
        strain = "LAI" if "LAI" in post_name else "1203"
        pre_name = f"mergedpre{strain}"

        pre = counts_df[counts_df["sample"] == pre_name].copy()
        post = counts_df[counts_df["sample"] == post_name].copy()
        merged = pre.merge(post, on=["site","wildtype","codon","wt_aa","aa"], suffixes=("_pre","_post"))

        keep = (merged["counts_pre"] > 0) & ((merged["counts_pre"] > 15) | (merged["counts_post"] > 15))
        merged = merged[keep].copy()

        merged["mut_ratio"] = merged["freqs_w_p_post"] / merged["freqs_w_p_pre"].replace(0, np.nan)
        merged["wt_ratio"]  = merged["wt_freq_w_p_post"] / merged["wt_freq_w_p_pre"].replace(0, np.nan)
        merged["enrichment_ratio"]    = merged["mut_ratio"] / merged["wt_ratio"]
        merged["log_enrichment_ratio"] = np.log2(merged["enrichment_ratio"].replace(0, np.nan))
        merged["log2_enrichment_relWT"] = merged["log_enrichment_ratio"] - np.log2(merged["wt_ratio"].replace(0, np.nan))

        merged["wt_ok_for_centering"] = (merged["wt_freq_w_p_pre"] >= WT_MIN_FREQ) & (merged["wt_freq_w_p_post"] >= WT_MIN_FREQ)

        merged["sample_post"] = post_name
        merged["homolog_pre"] = strain
        enrich_parts.append(merged)

    enrich_df = pd.concat(enrich_parts, ignore_index=True)
    enrich_df["site"] = enrich_df["site"].astype(int)
    enrich_df.to_csv("enrich_df.csv", index=False)

# -------------------------
# Stage 3: replicate-qc
# -------------------------
def cmd_replicate_qc(args):
    df = pd.read_csv("enrich_df.csv")
    SCORE_COL = "log2_enrichment_relWT" if USE_WT_CENTERING else "log_enrichment_ratio"

    df_sc = df.copy()
    df_sc = df_sc[(df_sc["aa"] != df_sc["wt_aa"]) & (df_sc["aa"] != "*")]
    df_sc = df_sc[df_sc["site"].between(SITE_MIN, SITE_MAX)]
    df_sc = df_sc[df_sc["wt_ok_for_centering"]] if USE_WT_CENTERING else df_sc

    df_sc["sample_group"] = df_sc["sample_post"].apply(lambda s: "LAI" if "LAI" in s else "1203")
    codon_long = (
        df_sc.assign(key=lambda d: d["site"].astype(str) + "|" + d["codon"])
             [["sample_post","sample_group","key",SCORE_COL]]
             .rename(columns={SCORE_COL: "score"})
    )
    P_codon = codon_long.pivot(index="key", columns="sample_post", values="score")

    site_long = (
        df_sc.groupby(["sample_post","sample_group","site"], as_index=False)[SCORE_COL]
             .median()
             .assign(key=lambda d: d["site"].astype(str))
             [["sample_post","sample_group","key",SCORE_COL]]
             .rename(columns={SCORE_COL: "score"})
    )
    P_site = site_long.pivot(index="key", columns="sample_post", values="score")

    def _pairs(cols, filt):
        arr = [c for c in cols if filt in c]
        return [(arr[i], arr[j]) for i in range(len(arr)) for j in range(i+1, len(arr))]

    samples_all = sorted(P_codon.columns)
    pair_lai = _pairs(samples_all, "LAI")
    pair_1203 = _pairs(samples_all, "1203")
    pair_inter = [(a,b) for a in samples_all for b in samples_all if ("LAI" in a and "1203" in b)]

    def _scatter_grid(P, pair_list, title, fname):
        if not pair_list: return
        n = len(pair_list); ncols = min(3, n); nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
        axes = np.array(axes).ravel()
        for ax, (a,b) in zip(axes, pair_list):
            x, y = P.get(a), P.get(b)
            if x is None or y is None: ax.axis("off"); continue
            valid = x.notna() & y.notna()
            if valid.sum() < 3: ax.axis("off"); continue
            xx, yy = x[valid].astype(float).values, y[valid].astype(float).values
            lo, hi = np.nanpercentile(np.concatenate([xx,yy]), [1, 99])
            pad = 0.05 * (hi - lo if hi > lo else 1.0)
            ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
            ax.scatter(xx, yy, s=10, alpha=0.5, edgecolors="none")
            ax.plot([lo, hi],[lo, hi], ls="--", c="gray")
            rho_s = spearmanr(xx, yy).correlation
            rho_p = np.corrcoef(xx, yy)[0,1]
            ax.set_title(f"{a} vs {b}\nSpearman={rho_s:.2f}, Pearson={rho_p:.2f}, n={valid.sum()}")
        for ax in axes[len(pair_list):]:
            ax.axis("off")
        plt.suptitle(title); plt.tight_layout()
        plt.savefig(FIG_DIR / fname, bbox_inches="tight"); plt.close()

    def _sign_agree(P, pair_list, title, fname):
        rows = []
        for a,b in pair_list:
            x, y = P.get(a), P.get(b)
            if x is None or y is None: continue
            valid = x.notna() & y.notna()
            if valid.sum() < 3: continue
            xx, yy = x[valid].astype(float).values, y[valid].astype(float).values
            pct = 100 * np.mean(np.sign(xx) == np.sign(yy))
            rows.append({"pair": f"{a} vs {b}", "sign_agree_pct": pct, "n": int(valid.sum())})
        if not rows: return
        dfm = pd.DataFrame(rows)
        plt.figure(figsize=(max(6, 1.2*len(rows)), 3.5))
        plt.bar(dfm["pair"], dfm["sign_agree_pct"])
        plt.xticks(rotation=45, ha="right"); plt.ylim(0, 100)
        plt.ylabel("Sign Agreement (%)"); plt.title(title); plt.tight_layout()
        plt.savefig(FIG_DIR / fname, bbox_inches="tight"); plt.close()

    FIG_DIR.mkdir(exist_ok=True, parents=True)
    _scatter_grid(P_codon, pair_lai,   "Codon-level LAI",     "scatter_variant_codon_LAI.pdf")
    _scatter_grid(P_codon, pair_1203,  "Codon-level 1203",    "scatter_variant_codon_1203.pdf")
    _scatter_grid(P_codon, pair_inter, "Codon-level Inter",   "scatter_variant_codon_interstrain.pdf")
    _sign_agree(P_codon, pair_lai,   "Codon Sign Agreement LAI",   "sign_agree_codon_LAI.pdf")
    _sign_agree(P_codon, pair_1203,  "Codon Sign Agreement 1203",  "sign_agree_codon_1203.pdf")
    _sign_agree(P_codon, pair_inter, "Codon Sign Agreement Inter", "sign_agree_codon_interstrain.pdf")

# -------------------------
# Stage 4: site-stats
# -------------------------
def cmd_site_stats(args):
    df = pd.read_csv("enrich_df.csv")
    df["log_enrichment_orig"] = df["log_enrichment_ratio"]
    if USE_WT_CENTERING:
        df = df[df["wt_ok_for_centering"].fillna(False)]
        df["log_enrichment_ratio"] = df["log2_enrichment_relWT"]

    def _rep_from_name(name: str):
        m = re.search(r"(R\d+)$", str(name)); return m.group(1) if m else None
    def _strain_from_name(name: str):
        name = str(name)
        return "LAI" if "LAI" in name else ("1203" if "1203" in name else "unknown")

    df["replicate"]    = df["sample_post"].apply(_rep_from_name)
    df["sample_group"] = df["sample_post"].apply(_strain_from_name)
    df = df[(df["aa"] != df["wt_aa"]) & (df["aa"] != "*")]
    df = df[df["site"].between(SITE_MIN, SITE_MAX)]

    # dynamic replicate sets
    sample_mapping = {grp: sorted(df[df["sample_group"]==grp]["sample_post"].dropna().unique().tolist())
                      for grp in ["LAI","1203"]}
    df_plot = pd.concat([df[(df["sample_group"]==grp) & (df["sample_post"].isin(reps))]
                         for grp, reps in sample_mapping.items()], ignore_index=True)

    # Unique variant counts per site
    unique_counts = (df_plot.groupby(["sample_group","site"])["aa"].nunique()
                            .rename("unique_variant_count").reset_index())
    unique_counts.to_excel("unique_variant_counts_median.xlsx", index=False)
    ok_sites = {
        grp: set(unique_counts.query("sample_group == @grp and unique_variant_count >= @MIN_UNIQUE_VARIANTS_FOR_STATS")["site"])
        for grp in ["LAI","1203"]
    }
    all_sites = list(range(SITE_MIN, SITE_MAX+1))

    # Split by strain
    df_lai = df_plot[df_plot["sample_group"]=="LAI"].copy()
    df_1203 = df_plot[df_plot["sample_group"]=="1203"].copy()

    LAI_median = float(np.nanmedian(df_lai["log_enrichment_ratio"])) if not df_lai.empty else np.nan
    D1203_median = float(np.nanmedian(df_1203["log_enrichment_ratio"])) if not df_1203.empty else np.nan

    def _overall_site_level_median(df_group: pd.DataFrame, ok_site_set: set) -> float:
        site_meds = df_group.groupby("site")["log_enrichment_ratio"].median()
        if ok_site_set:
            site_meds = site_meds[site_meds.index.isin(ok_site_set)]
        return float(np.nanmedian(site_meds)) if len(site_meds) > 0 else np.nan

    LAI_sitelevel_med = _overall_site_level_median(df_lai, ok_sites["LAI"])
    D1203_sitelevel_med = _overall_site_level_median(df_1203, ok_sites["1203"])

    # Per-site Wilcoxon vs global + BH-FDR
    alpha = 0.05; fdr_thresh = 0.10
    def _n_unique_site(df_group): return (df_group.groupby("site")["aa"].nunique()).to_dict()

    def site_wilcoxon_table(df_group: pd.DataFrame, global_median: float, ok_site_set: set, n_unique_map: dict) -> pd.DataFrame:
        rows = []
        for site in all_sites:
            vals = df_group.loc[df_group["site"] == site, "log_enrichment_ratio"].to_numpy()
            nuq = int(n_unique_map.get(site, 0))
            use = (site in ok_site_set) and (vals.size > 0)
            if not use:
                rows.append({"site":site, "site_median":np.nan, "overall_median":global_median,
                             "delta_median":np.nan, "p_less":np.nan, "p_greater":np.nan,
                             "p_two":np.nan, "n_eff":0, "direction":"Equal",
                             "p_sig_nominal_two":False, "n_unique":nuq, "use_for_stats":False})
                continue
            diffs = vals - global_median
            diffs = diffs[~np.isnan(diffs)]; nonzero = diffs[diffs != 0]
            n_eff = nonzero.size
            if n_eff == 0:
                p_less = p_greater = p_two = np.nan
            else:
                p_less    = wilcoxon(diffs, alternative="less", zero_method="wilcox").pvalue
                p_greater = wilcoxon(diffs, alternative="greater", zero_method="wilcox").pvalue
                p_two     = wilcoxon(diffs, alternative="two-sided", zero_method="wilcox").pvalue
            site_median = float(np.nanmedian(vals))
            delta = site_median - global_median
            direction = "Lower" if np.isfinite(delta) and delta < 0 else ("Higher" if np.isfinite(delta) and delta > 0 else "Equal")
            rows.append({"site":site, "site_median":site_median, "overall_median":global_median,
                         "delta_median":delta, "p_less":p_less, "p_greater":p_greater,
                         "p_two":p_two, "n_eff":n_eff, "direction":direction,
                         "p_sig_nominal_two": (p_two < alpha) if np.isfinite(p_two) else False,
                         "n_unique":nuq, "use_for_stats":True})
        out = pd.DataFrame(rows)
        out["q_two_bh"] = np.nan
        m2 = out["use_for_stats"] & out["p_two"].notna()
        if m2.any():
            out.loc[m2, "q_two_bh"] = multipletests(out.loc[m2, "p_two"], method="fdr_bh")[1]
        out["sig_two_bh"] = out["q_two_bh"] <= fdr_thresh

        out["p_dir"] = np.where(out["direction"]=="Lower", out["p_less"],
                         np.where(out["direction"]=="Higher", out["p_greater"], np.nan))
        out["q_dir_bh"] = np.nan
        md = out["use_for_stats"] & out["p_dir"].notna()
        if md.any():
            out.loc[md, "q_dir_bh"] = multipletests(out.loc[md, "p_dir"], method="fdr_bh")[1]
        out["sig_dir_bh"] = out["q_dir_bh"] <= fdr_thresh
        return out

    lai_nu = _n_unique_site(df_lai); d1203_nu = _n_unique_site(df_1203)
    lai_sig = site_wilcoxon_table(df_lai, LAI_median, ok_sites["LAI"], lai_nu);    lai_sig["sample_group"] = "LAI"
    d1203_sig = site_wilcoxon_table(df_1203, D1203_median, ok_sites["1203"], d1203_nu); d1203_sig["sample_group"] = "1203"
    significance_results = pd.concat([lai_sig, d1203_sig], ignore_index=True)
    significance_results.to_csv("site_significance_full.csv", index=False)

    # Replicate-consistent gated lists
    key_cols = ["sample_group","site","direction"]
    sig_two = significance_results[(significance_results["use_for_stats"]) &
                                   (significance_results["sig_two_bh"]) &
                                   (significance_results["direction"]!="Equal")][key_cols].drop_duplicates()
    sig_dir = significance_results[(significance_results["use_for_stats"]) &
                                   (significance_results["sig_dir_bh"]) &
                                   (significance_results["direction"]!="Equal")][key_cols].drop_duplicates()

    # per-rep medians and depth
    gmed = (df_plot.groupby("sample_post")["log_enrichment_ratio"].median().rename("gmed").reset_index())
    site_rep = (df_plot.groupby(["sample_group","sample_post","site"])["log_enrichment_ratio"]
                       .median().rename("site_med").reset_index())
    M = site_rep.merge(gmed, on="sample_post", how="left")
    M["delta"] = M["site_med"] - M["gmed"]
    M["direction"] = np.where(M["delta"]>0, "Higher", np.where(M["delta"]<0, "Lower", "Equal"))
    if {"counts_pre","counts_post"}.issubset(df_plot.columns):
        depth_rep_site = (df_plot.assign(depth=df_plot["counts_pre"] + df_plot["counts_post"])
                                .groupby(["sample_group","sample_post","site"])["depth"].median()
                                .rename("med_depth").reset_index())
        M = M.merge(depth_rep_site, on=["sample_group","sample_post","site"], how="left")
    else:
        M["med_depth"] = np.nan

    ok = M[M["med_depth"] >= DEPTH_TH_FOR_CONSISTENCY]
    rows = []
    for grp, ggrp in ok.groupby("sample_group"):
        for site, g in ggrp.groupby("site"):
            vc = g["direction"].value_counts()
            if vc.get("Higher", 0) >= 2: rows.append((grp, int(site), "Higher"))
            if vc.get("Lower",  0) >= 2: rows.append((grp, int(site), "Lower"))
    consis = pd.DataFrame(rows, columns=key_cols).drop_duplicates()

    def _gate(callset, cons):
        if callset is None or callset.empty or cons is None or cons.empty:
            return pd.DataFrame(columns=key_cols)
        return callset.merge(cons, on=key_cols, how="inner").drop_duplicates()

    both_g = _gate(sig_two.merge(sig_dir, on=key_cols, how="inner"), consis)
    dir_only_g = _gate(sig_dir.merge(sig_two.assign(_in_two=True), on=key_cols, how="left").query("_in_two.isna()")[key_cols], consis)
    two_only_g = _gate(sig_two.merge(sig_dir.assign(_in_dir=True), on=key_cols, how="left").query("_in_dir.isna()")[key_cols], consis)

    def _union_tagged(both, donly, tonly):
        rows = []
        for dfc, tag in [(both,"both"), (donly,"directional_only"), (tonly,"two_sided_only")]:
            if dfc is not None and not dfc.empty:
                t = dfc.copy(); t["support"] = tag; rows.append(t)
        return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=key_cols+["support"])

    summary_gated = _union_tagged(both_g, dir_only_g, two_only_g).sort_values(["sample_group","direction","site"])
    summary_gated.to_csv("replicate_consistent_BH_calls.csv", index=False)

    # Violin plots (matplotlib only to keep deps light)
    FIG_DIR.mkdir(exist_ok=True, parents=True)
    def _violin_one(sub, color, title, out_pdf, median_val, bad_set):
        # Manual violins with matplotlib (no seaborn dependency)
        plt.figure(figsize=(15,5))
        # Build data per site
        sites = sorted(sub["site"].unique().tolist())
        data = [sub.loc[sub["site"]==s, "log_enrichment_ratio"].values for s in sites]
        parts = plt.violinplot(data, positions=np.arange(len(sites)), widths=0.9, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(color); pc.set_edgecolor("black"); pc.set_alpha(0.8)
        if np.isfinite(median_val):
            plt.axhline(median_val, color="red", linestyle="--")
        plt.xticks(np.arange(len(sites))[::5], [str(sites[i]) for i in range(0, len(sites), 5)], rotation=45)
        # mark excluded sites with 'n' near bottom
        y0, y1 = plt.ylim(); y_text = y0 + 0.03*(y1-y0)
        for i, s in enumerate(sites):
            if int(s) in bad_set:
                plt.text(i, y_text, "n", ha="center", va="bottom", fontsize=8)
        plt.xlabel("Site"); plt.ylabel("Log2 enrichment" + (" (WT-centered)" if USE_WT_CENTERING else ""))
        plt.title(title); plt.tight_layout(); plt.savefig(out_pdf, bbox_inches="tight"); plt.close()

    bad_LAI  = set(all_sites) - ok_sites["LAI"]
    bad_1203 = set(all_sites) - ok_sites["1203"]
    _violin_one(df_lai,  GREEN,  "LAI: Log enrichment by site (pooled pre reference)",  FIG_DIR / "sitelevellogenrichment_violin_median_LAI.pdf", LAI_median,  bad_LAI)
    _violin_one(df_1203, PURPLE, "1203: Log enrichment by site (pooled pre reference)", FIG_DIR / "sitelevellogenrichment_violin_median_1203.pdf", D1203_median, bad_1203)

# -------------------------
# Stage 5: conservation
# -------------------------
def _read_fasta_align(path: Path):
    names, seqs = [], []
    cur_name, cur = None, []
    with open(path, "r") as f:
        for line in f:
            t = line.strip()
            if not t: continue
            if t.startswith(">"):
                if cur_name is not None:
                    names.append(cur_name); seqs.append("".join(cur))
                cur_name, cur = t[1:].split()[0], []
            else:
                cur.append(t)
        if cur_name is not None:
            names.append(cur_name); seqs.append("".join(cur))
    return names, seqs

def cmd_conservation(args):
    align_path = Path(args.alignment)
    if not align_path.exists():
        raise SystemExit(f"Missing alignment FASTA: {align_path}")
    names, seqs = _read_fasta_align(align_path)
    if not seqs: raise SystemExit("No sequences in alignment.")
    L = len(seqs[0])
    if any(len(s)!=L for s in seqs):
        raise SystemExit("All alignment sequences must be same length.")

    def pct_vs_modal(col_chars):
        non_gap = [c for c in col_chars if c not in ("-",".")]
        if not non_gap: return np.nan
        modal, _ = Counter(non_gap).most_common(1)[0]
        return 100.0 * sum(c==modal for c in non_gap) / len(non_gap)

    def shannon_entropy(col_chars):
        non_gap = [c for c in col_chars if c not in ("-",".")]
        if not non_gap: return np.nan
        total = len(non_gap); counts = Counter(non_gap)
        return -sum((cnt/total)*math.log(cnt/total) for cnt in counts.values())

    pct_id, Hs = [], []
    for i in range(L):
        col = [s[i] for s in seqs]
        pct_id.append(pct_vs_modal(col)); Hs.append(shannon_entropy(col))

    df_cons = pd.DataFrame({"site": np.arange(1, L+1, dtype=int),
                            "pct_identity_consensus": pct_id,
                            "H_nat_nats": Hs})
    # get per-site medians from site_stats stage
    if Path("site_significance_full.csv").exists():
        sig = pd.read_csv("site_significance_full.csv")
        df_lai_meds  = sig[["site","site_median","sample_group"]].query("sample_group=='LAI'")  .rename(columns={"site_median":"LAI_median"}).drop(columns=["sample_group"]).drop_duplicates("site")
        df_1203_meds = sig[["site","site_median","sample_group"]].query("sample_group=='1203'").rename(columns={"site_median":"1203_median"}).drop(columns=["sample_group"]).drop_duplicates("site")
    else:
        # fallback
        enr = pd.read_csv("enrich_df.csv")
        if USE_WT_CENTERING:
            enr = enr[enr["wt_ok_for_centering"].fillna(False)]
            enr["log_enrichment_ratio"] = enr["log2_enrichment_relWT"]
        df_lai_meds  = enr[enr["sample_post"].str.contains("LAI", na=False)].groupby("site", as_index=False)["log_enrichment_ratio"].median().rename(columns={"log_enrichment_ratio": "LAI_median"})
        df_1203_meds = enr[enr["sample_post"].str.contains("1203", na=False)].groupby("site", as_index=False)["log_enrichment_ratio"].median().rename(columns={"log_enrichment_ratio": "1203_median"})

    df_plot = (df_cons.merge(df_lai_meds, on="site", how="left")
                      .merge(df_1203_meds, on="site", how="left")
                      .query("@SITE_MIN <= site <= @SITE_MAX").sort_values("site"))

    def _smooth(x, window=9):
        s = pd.Series(x, index=df_plot["site"])
        return s.rolling(window=window, center=True, min_periods=1).mean().values

    df_plot["pctid_consensus_smooth"] = _smooth(df_plot["pct_identity_consensus"])
    df_plot["H_nat_nats_smooth"] = _smooth(df_plot["H_nat_nats"])

    def _dual_axis(dfp, strain, color, mode, out_svg):
        ycol = f"{strain}_median"
        cons_id = ("pct_identity_consensus", "pctid_consensus_smooth")
        cons_H = ("H_nat_nats", "H_nat_nats_smooth")

        def _pick(base_col, smooth_col, mode):
            if mode=="raw": return [(base_col, dict(alpha=0.9, ls="-", lw=1.5, label_suffix=""))]
            if mode=="smooth": return [(smooth_col, dict(alpha=1.0, ls="--", lw=1.8, label_suffix=" (smoothed)"))]
            return [(base_col, dict(alpha=0.25, ls="-", lw=1.0, label_suffix=" (raw)")),
                    (smooth_col, dict(alpha=1.0, ls="--", lw=1.8, label_suffix=" (smoothed)"))]

        # % identity panel
        plt.figure(figsize=(11,3.8))
        ax = plt.gca()
        sub = dfp[["site", ycol]].dropna()
        ax.plot(sub["site"], sub[ycol], color=color, label=f"{strain} median enrichment")
        ax.set_xlabel("Vif residue (site)"); ax.set_ylabel("Median log2 enrichment")
        ax2 = ax.twinx()
        for col, st in _pick(*cons_id, mode):
            ax2.plot(dfp["site"], dfp[col], color="black", alpha=st["alpha"], ls=st["ls"], lw=st["lw"],
                     label="% identity"+st["label_suffix"])
        ax2.set_ylabel("Percent identity (%)"); ax2.set_ylim(100, 0)
        h1,l1 = ax.get_legend_handles_labels(); h2,l2 = ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="center left", bbox_to_anchor=(1.01,0.5))
        plt.tight_layout(); plt.savefig(SUPP_DIR / out_svg, bbox_inches="tight"); plt.close()

        # Entropy panel
        plt.figure(figsize=(11,3.8))
        ax = plt.gca()
        ax.plot(sub["site"], sub[ycol], color=color, label=f"{strain} median enrichment")
        ax.set_xlabel("Vif residue (site)"); ax.set_ylabel("Median log2 enrichment")
        ax2 = ax.twinx()
        for col, st in _pick(*cons_H, mode):
            ax2.plot(dfp["site"], dfp[col], color="black", alpha=st["alpha"], ls=st["ls"], lw=st["lw"],
                     label="Shannon entropy"+st["label_suffix"])
        ax2.set_ylabel("Entropy (nats)")
        h1,l1 = ax.get_legend_handles_labels(); h2,l2 = ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="center left", bbox_to_anchor=(1.01,0.5))
        plt.tight_layout(); plt.savefig(SUPP_DIR / out_svg.replace("_pctid_", "_entropy_"), bbox_inches="tight"); plt.close()

    FIG_DIR.mkdir(exist_ok=True, parents=True); SUPP_DIR.mkdir(exist_ok=True, parents=True)
    _dual_axis(df_plot, "LAI",  GREEN,  args.mode, "LAI_median_vs_pctid_{mode}.svg".format(mode=args.mode))
    _dual_axis(df_plot, "1203", PURPLE, args.mode, "1203_median_vs_pctid_{mode}.svg".format(mode=args.mode))

# -------------------------
# Stage 6: cross-strain
# -------------------------
def cmd_cross_strain(args):
    df = pd.read_csv("enrich_df.csv")
    if USE_WT_CENTERING:
        df = df[df["wt_ok_for_centering"].fillna(False)]
        df["log_enrichment_ratio"] = df["log2_enrichment_relWT"]
    df = df[(df["aa"] != df["wt_aa"]) & (df["aa"] != "*")]
    df = df[df["site"].between(SITE_MIN, SITE_MAX)]

    df["sample_group"] = df["sample_post"].apply(lambda s: "LAI" if "LAI" in s else ("1203" if "1203" in s else "unknown"))
    df_lai = df[df["sample_group"]=="LAI"]
    df_1203= df[df["sample_group"]=="1203"]

    lai_meds = df_lai.groupby("site")["log_enrichment_ratio"].median().reset_index(name="median_value_LAI")
    s1203_meds = df_1203.groupby("site")["log_enrichment_ratio"].median().reset_index(name="median_value_1203")
    site_summary = lai_meds.merge(s1203_meds, on="site", how="inner")

    def perform_test(site):
        x = df_lai[df_lai["site"]==site]["log_enrichment_ratio"].values
        y = df_1203[df_1203["site"]==site]["log_enrichment_ratio"].values
        if len(x)>0 and len(y)>0:
            stat, p = mannwhitneyu(x, y, alternative="two-sided")
        else:
            stat, p = np.nan, np.nan
        return pd.Series({"u_stat": stat, "p_value": p})

    test_results = site_summary["site"].apply(perform_test)
    site_summary = pd.concat([site_summary, test_results], axis=1)
    site_summary["significant"] = site_summary["p_value"] < 0.05
    site_summary["direction"] = np.where(site_summary["median_value_LAI"] > site_summary["median_value_1203"],
                                         "Higher in LAI", "Higher in 1203")

    # attach WT amino acids per strain
    wt_lai = (df[df["sample_post"].str.contains("LAI", na=False)]
                .groupby("site")["wt_aa"].first().reset_index().rename(columns={"wt_aa":"wt_aa_LAI"}))
    wt_1203= (df[df["sample_post"].str.contains("1203", na=False)]
                .groupby("site")["wt_aa"].first().reset_index().rename(columns={"wt_aa":"wt_aa_1203"}))
    site_summary = site_summary.merge(wt_lai, on="site", how="left").merge(wt_1203, on="site", how="left")

    sig = site_summary[site_summary["significant"]].copy()
    sig["abs_median_diff"] = np.abs(sig["median_value_LAI"] - sig["median_value_1203"])
    cols = ["site","wt_aa_LAI","wt_aa_1203","median_value_LAI","median_value_1203",
            "abs_median_diff","p_value","direction"]
    sig[cols].to_excel("significant_site_differences.xlsx", index=False)

    # Plot
    FIG_DIR.mkdir(exist_ok=True, parents=True)
    plt.figure(figsize=(4,4))
    for _, r in site_summary.iterrows():
        plt.scatter(r["median_value_LAI"], r["median_value_1203"],
                    s=60, facecolor=("black" if r["significant"] else "white"),
                    edgecolor="black", alpha=0.9)
    plt.xlabel("LAI Site Median"); plt.ylabel("1203 Site Median")
    lims = [min(plt.xlim()[0], plt.ylim()[0]), max(plt.xlim()[1], plt.ylim()[1])]
    plt.plot(lims, lims, "k--", alpha=0.5); plt.xlim(lims); plt.ylim(lims)
    for _, r in site_summary.iterrows():
        if r["significant"]:
            plt.annotate(str(int(r["site"])), (r["median_value_LAI"], r["median_value_1203"]),
                         textcoords="offset points", xytext=(4,4), fontsize=8, color="darkred")
    plt.tight_layout(); plt.savefig(FIG_DIR / "site_median_comparison.pdf", bbox_inches="tight"); plt.close()

# -------------------------
# CLI
# -------------------------
def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Site-level analysis pipeline (publication-ready CLI).")
    p.add_argument("--log-level", default="INFO")
    sub = p.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("prepare", help="Estimate background, filter codoncounts, write pooled PRE tables.")
    s2 = sub.add_parser("enrich", help="Build counts_df and enrich_df (variant enrichment vs pooled pre).")
    s3 = sub.add_parser("replicate-qc", help="Replicate correlation/sign-agreement grids.")
    s4 = sub.add_parser("site-stats", help="Violins + per-site Wilcoxon and BH-FDR; replicate-consistent lists.")
    s5 = sub.add_parser("conservation", help="Dual-axis plots vs consensus %identity / Shannon entropy.")
    s5.add_argument("--alignment", type=Path, default=Path("subtype_ref_protein.fasta"))
    s5.add_argument("--mode", choices=["raw","smooth","overlay"], default="raw", help="Line mode for conservation tracks.")
    s6 = sub.add_parser("cross-strain", help="Mann–Whitney across strains on site medians + plot + Excel.")

    return p.parse_args(list(argv) if argv is not None else None)

def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    if args.cmd == "prepare":       cmd_prepare(args)
    elif args.cmd == "enrich":      cmd_enrich(args)
    elif args.cmd == "replicate-qc":cmd_replicate_qc(args)
    elif args.cmd == "site-stats":  cmd_site_stats(args)
    elif args.cmd == "conservation":cmd_conservation(args)
    elif args.cmd == "cross-strain":cmd_cross_strain(args)
    else:
        raise SystemExit("Unknown command")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
