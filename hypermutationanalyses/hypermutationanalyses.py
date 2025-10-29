#!/usr/bin/env python3
"""
hypermutationanalyses.py

Order of operations:
  1) calculate  -> read FASTA under `--fasta-dir`, compute metrics, write CSVs
  2) analyze    -> read those CSVs + variant-proportion tables, make plots/CSV

Dependencies: pandas, numpy, scipy, matplotlib, biopython, scikit-learn
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import glob
import gzip
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import linregress, spearmanr, pearsonr, ttest_ind
from sklearn.metrics import r2_score

plt.switch_backend("Agg")


# -------------------------
# Shared config
# -------------------------

REF_SEQS: Dict[str, str] = {
    "LAI": "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTT",
    "1203": "ATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTGGACAGGATGAAGATTAAAACATGGAAAAGTTTAGTAAAGCATCATATGTATGTTTCAAAGAAGGCTAGGAGATGGTTTTATAGACATCACTATGAAAGCACTCATCCAAAAATAAGTTCAGAAGTACACATCCCACTAGAGAAGGGTGAATTGGTAGTAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAGGGAGATATAGCACACAAGTAGACCCTGACCTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTT",
}

# Pairwise GG->AG scan (with context rule) used to build allowable variant lists per site
MUTATION_TARGETS = [('GG', 'AG')]

def _check_downstream_context(sequence: str, i: int) -> bool:
    # original rule: skip if +2 is 'C'
    if i + 2 < len(sequence):
        return sequence[i + 2] != 'C'
    return False

def _find_mutation_effects(sequence: str) -> Tuple[List[Tuple[int, str, str, str, str]], List[Tuple[int, str, str, str, str]]]:
    change, no_change = [], []
    for i in range(len(sequence) - 1):
        for orig, mut in MUTATION_TARGETS:
            if sequence[i:i+2] == orig and _check_downstream_context(sequence, i):
                mseq = sequence[:i] + mut + sequence[i+2:]
                codon_start = (i // 3) * 3
                codon = sequence[codon_start:codon_start+3]
                m_codon = mseq[codon_start:codon_start+3]
                if len(codon) == 3 and len(m_codon) == 3:
                    aa = str(Seq(codon).translate())
                    maa = str(Seq(m_codon).translate())
                    tup = (codon_start // 3 + 1, codon, m_codon, aa, maa)
                    (change if aa != maa else no_change).append(tup)
    return change, no_change

def _save_mut_effects(effects: List[Tuple[int, str, str, str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as csvfile:
        w = csv.DictWriter(csvfile, fieldnames=["site", "wildtype", "variant", "wildtype_aa", "variant_aa"])
        w.writeheader()
        for site, wt, var, waa, vaa in effects:
            w.writerow({"site": site, "wildtype": wt, "variant": var, "wildtype_aa": waa, "variant_aa": vaa})


# -------------------------
# Stage 1: CALCULATE
# -------------------------

def _read_fasta_all(fp: Path) -> List[str]:
    seqs = []
    with fp.open("r") as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).strip().upper())
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur).strip().upper())
    return seqs

def _aa_site_to_nt_index(site: int) -> int:
    return 3 * (site - 1)

def _count_specified_variants(sequence: str, specified: pd.DataFrame, exclusion_site_nt: Optional[int]) -> Tuple[int, List[Tuple[int, int, str]]]:
    count = 0
    matches: List[Tuple[int, int, str]] = []
    for _, row in specified.iterrows():
        aa_site = int(row["site"])
        pos = _aa_site_to_nt_index(aa_site)
        if exclusion_site_nt is not None and pos == exclusion_site_nt:
            continue
        var = row["variant"]
        if pos + len(var) <= len(sequence) and sequence[pos:pos+len(var)] == var:
            count += 1
            matches.append((aa_site, pos + 1, var))  # 1-based nt index
    return count, matches

def _count_specific_sites(sequence: str, specified: pd.DataFrame, aa_sites: List[int]) -> int:
    cnt = 0
    for aa_site in aa_sites:
        row = specified[specified["site"] == aa_site]
        if row.empty:
            continue
        pos = _aa_site_to_nt_index(aa_site)
        var = row.iloc[0]["variant"]
        if pos + len(var) <= len(sequence) and sequence[pos:pos+len(var)] == var:
            cnt += 1
    return cnt

def _extract_details_from_fname(filename: str) -> Tuple[str, str, str, str, str, str, int, str]:
    # {sample_type}{sample_name}_filtered_{codon}_{wildtype}_{site}_{aa_mutation}.fasta
    parts = filename.replace(".fasta", "").split("_")
    if len(parts) != 6:
        raise ValueError(f"Unexpected filename format: {filename}")
    if parts[0].startswith("pre"):
        stype = "pre"
    elif parts[0].startswith("post"):
        stype = "post"
    else:
        raise ValueError(f"Filename must start with pre/post: {filename}")
    sample_name = parts[0]
    homolog = "LAI" if "LAI" in sample_name else "1203"
    repl = sample_name.replace(stype, "").replace(homolog, "")
    codon = parts[2]
    wt = parts[3]
    try:
        site = int(parts[4])
    except Exception as e:
        raise ValueError(f"Non-numeric site in filename {filename}") from e
    aa_mut = parts[5]
    return sample_name, stype, homolog, repl, wt, codon, site, aa_mut

def _process_fasta_file(file_name: str, folder: Path, enrich_df: pd.DataFrame, lai_spec: pd.DataFrame, spec_1203: pd.DataFrame):
    specific_sites = [21, 38, 70, 89]
    sample_name, stype, homolog, repl, wt, codon, site, aa_mut = _extract_details_from_fname(file_name)
    refseq = REF_SEQS[homolog]

    e = enrich_df.loc[
        (enrich_df["sample_post"] == sample_name) &
        (enrich_df["site"] == site) &
        (enrich_df["codon"] == codon) &
        (enrich_df["wildtype"] == wt) &
        (enrich_df["aa_mutation"] == aa_mut)
    ]
    if e.empty:
        return None

    specified = lai_spec if homolog == "LAI" else spec_1203

    n_seq = 0
    total_var = 0
    seq_with_var = 0
    seq_with_ge3 = 0
    seq_with_ge1 = 0
    total_specific_sites = 0
    matched: List[Tuple[int,int,str]] = []

    seqs = _read_fasta_all(folder / file_name)
    for s in seqs:
        n_seq += 1
        vcount, vlist = _count_specified_variants(s, specified, _aa_site_to_nt_index(site))
        total_var += vcount
        total_specific_sites += _count_specific_sites(s, specified, specific_sites)
        if vcount > 0:
            seq_with_var += 1
        if vcount >= 3:
            seq_with_ge3 += 1
        if vcount >= 1:
            seq_with_ge1 += 1
        matched.extend(vlist)

    einfo = e.iloc[0].to_dict()
    einfo.update({
        "Replicate": repl,
        "Total Sequences": n_seq,
        "Total Specified Variants": total_var,
        "Total W Site Variants": total_specific_sites,
        "Sequences with 3 or More Variants": seq_with_ge3,
        "Sequences with 1 or More Variants": seq_with_ge1,
        "Normalized with 1 or more variants": (seq_with_ge1 / n_seq) if n_seq > 0 else 0.0,
        "Normalized with 3 or more variants": (seq_with_ge3 / n_seq) if n_seq > 0 else 0.0,
        "Normalized total variants": (total_var / n_seq) if n_seq > 0 else 0.0,
        "Matched Variants": matched,
    })
    return einfo, homolog

def stage_calculate(
    fasta_dir: Path,
    enrich_csv: Path,
    out_dir: Path,
    lai_effects_csv: Path,
    s1203_effects_csv: Path,
    workers: int = 22
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) build GG->AG effects per site (both homologs)
    lai_change, lai_no = _find_mutation_effects(REF_SEQS["LAI"])
    s1203_change, s1203_no = _find_mutation_effects(REF_SEQS["1203"])
    _save_mut_effects(lai_change + lai_no, lai_effects_csv)
    _save_mut_effects(s1203_change + s1203_no, s1203_effects_csv)

    # 2) compute per-FASTA metrics
    enrich_df = pd.read_csv(enrich_csv)
    lai_spec = pd.read_csv(lai_effects_csv)
    spec_1203 = pd.read_csv(s1203_effects_csv)

    lai_rows, s1203_rows = [], []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, workers)) as ex:
        futs = {
            ex.submit(_process_fasta_file, fn, fasta_dir, enrich_df, lai_spec, spec_1203): fn
            for fn in os.listdir(fasta_dir) if fn.endswith(".fasta") and "post" in fn
        }
        for fut in concurrent.futures.as_completed(futs):
            out = fut.result()
            if out is None:
                continue
            einfo, homolog = out
            (lai_rows if homolog == "LAI" else s1203_rows).append(einfo)

    pd.DataFrame(lai_rows).to_csv(out_dir / "LAI_hypermutation_results.csv", index=False)
    pd.DataFrame(s1203_rows).to_csv(out_dir / "1203_hypermutation_results.csv", index=False)


# -------------------------
# Stage 2: ANALYZE
# -------------------------

def _filter_common_codons(df: pd.DataFrame) -> pd.DataFrame:
    # present in all replicates (nunique * 3 / 3 reduces to nunique but preserves original intention)
    rep_thresh = df["Replicate"].nunique()
    codon_counts = df.groupby(["site", "codon"])["Replicate"].nunique().reset_index(name="replicate_count")
    common = codon_counts[codon_counts["replicate_count"] >= rep_thresh][["site", "codon"]]
    return df.merge(common, on=["site", "codon"])

def _more_than_one_bp_diff(wt: str, codon: str) -> bool:
    return sum(1 for a, b in zip(wt, codon) if a != b) > 1

def _merge_with_variant_props(results_df: pd.DataFrame, var_props_df: pd.DataFrame) -> pd.DataFrame:
    return results_df.merge(var_props_df, left_on=["site","codon"], right_on=["site","variant"], how="left", indicator=True, suffixes=("", "_var"))

def _filter_codons(df: pd.DataFrame) -> pd.DataFrame:
    f = _filter_common_codons(df)
    f = f[f["counts_post"] > 0]
    # keep entries present in var-props OR >1bp from wt
    f = f[(f["_merge"] == "both") | (f.apply(lambda r: _more_than_one_bp_diff(r["wildtype"], r["codon"]), axis=1))]
    return f

def _prep_xy(df: pd.DataFrame, y_col: str) -> Tuple[pd.Series, pd.Series]:
    tmp = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log_enrichment_ratio", y_col]).copy()
    y = tmp[y_col]
    if y_col == "Percent with at least 1 A3G Mutation":
        tmp = tmp[tmp[y_col] > 0]
        y = tmp[y_col]
    return tmp["log_enrichment_ratio"], y

def _plot_scatter_fit(df: pd.DataFrame, title: str, out_pdf: Path, y_col: str, logy: bool) -> None:
    plt.figure(figsize=(4,3))
    x, y = _prep_xy(df, y_col)
    if y.empty:
        logging.warning("No points for %s", title)
        return
    plt.scatter(x, y, edgecolors="black", linewidths=0.25)
    if logy:
        slope, intercept, *_ = linregress(x, np.log(y))
        plt.plot(x, np.exp(slope * x + intercept), color="black", linewidth=2)
        plt.yscale("log")
        plt.ylim(10**-1.25, 10**2.25)
        plt.ylabel("% of Reads with at least 1 A3G Mutation")
        r2 = r2_score(np.log(y), slope * x + intercept)
    else:
        slope, intercept, *_ = linregress(x, y)
        plt.plot(x, slope * x + intercept, color="black", linewidth=2)
        plt.ylabel("Total A3G mutations / Total reads")
        r2 = r2_score(y, slope * x + intercept)
    plt.xlim(-5,8)
    plt.xticks(np.arange(-5,9,1))
    plt.xlabel("Log Enrichment Ratio")
    plt.title(title)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    pear, p_p = pearsonr(x, y)
    spear, s_p = spearmanr(x, y)
    logging.info("%s: n=%d, R2=%.3f, Pearson=%.3f (p=%.2e), Spearman=%.3f (p=%.2e)",
                 title, len(x), r2, pear, p_p, spear, s_p)

def _overall_fitness(df: pd.DataFrame, weight: float=1.0, scale_factor: float=10.0) -> pd.DataFrame:
    req = ["site", "Normalized total variants", "log_enrichment_ratio"]
    for c in req:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' for overall fitness")
    out = df.copy()
    out["overall_fitness"] = out["log_enrichment_ratio"] - weight * (scale_factor * out["Normalized total variants"])
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["overall_fitness", "site"])
    return out

def _plot_overall_fitness(df: pd.DataFrame, title: str, out_pdf: Path) -> None:
    if df.empty:
        logging.warning("No data for overall fitness: %s", title); return
    plt.figure(figsize=(6,3))
    plt.scatter(df["site"], df["overall_fitness"], edgecolors="black", linewidths=0.25)
    plt.xlabel("Site"); plt.ylabel("Overall Fitness"); plt.title(title)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

def stage_analyze(
    results_dir: Path,
    varprops_dir: Path,
    out_dir: Path
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Input results
    lai_res = pd.read_csv(results_dir / "LAI_hypermutation_results.csv")
    s1203_res = pd.read_csv(results_dir / "1203_hypermutation_results.csv")
    # Variant proportions (include *_with_wildtype.csv in varprops_dir)
    lai_vp = pd.read_csv(varprops_dir / "LAIvariantproportions_with_wildtype.csv")
    s1203_vp = pd.read_csv(varprops_dir / "1203variantproportions_with_wildtype.csv")

    # Merge + filter
    lai_m = _merge_with_variant_props(lai_res, lai_vp)
    s1203_m = _merge_with_variant_props(s1203_res, s1203_vp)
    lai_f = _filter_codons(lai_m)
    s1203_f = _filter_codons(s1203_m)

    # Derived columns used in the original analysis
    lai_f["Percent with at least 1 A3G Mutation"] = lai_f["Normalized with 1 or more variants"] * 100.0
    s1203_f["Percent with at least 1 A3G Mutation"] = s1203_f["Normalized with 1 or more variants"] * 100.0

    # Plots
    _plot_scatter_fit(lai_f, "LAI Data", out_dir / "LAI_hypermutation.pdf", "Percent with at least 1 A3G Mutation", logy=True)
    _plot_scatter_fit(s1203_f, "1203 Data", out_dir / "1203_hypermutation.pdf", "Percent with at least 1 A3G Mutation", logy=True)

    # Linear (non-log) flavor using 'Total A3G mutations / Total reads'
    if "Normalized total variants" in lai_f.columns:
        _plot_scatter_fit(lai_f.rename(columns={"Normalized total variants": "Total A3G mutations / Total reads"}),
                          "LAI Data (linear)", out_dir / "LAI_hypermutation_linear.pdf",
                          "Total A3G mutations / Total reads", logy=False)
    if "Normalized total variants" in s1203_f.columns:
        _plot_scatter_fit(s1203_f.rename(columns={"Normalized total variants": "Total A3G mutations / Total reads"}),
                          "1203 Data (linear)", out_dir / "1203_hypermutation_linear.pdf",
                          "Total A3G mutations / Total reads", logy=False)

    # Overall fitness (penalize mutation load)
    of_lai = _overall_fitness(lai_f)
    of_1203 = _overall_fitness(s1203_f)
    _plot_overall_fitness(of_lai, "LAI Overall Fitness", out_dir / "LAI_overall_fitness.pdf")
    _plot_overall_fitness(of_1203, "1203 Overall Fitness", out_dir / "1203_overall_fitness.pdf")


# -------------------------
# CLI
# -------------------------

def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Hypermutation pipeline: calculate -> analyze.")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])

    sub = p.add_subparsers(dest="cmd", required=True)

    # calculate
    pc = sub.add_parser("calculate", help="Compute hypermutation metrics from FASTA and enrichment CSV.")
    pc.add_argument("--fasta-dir", type=Path, required=True, help="Folder of per-mutation FASTA files (expects 'post' in filename).")
    pc.add_argument("--enrich-csv", type=Path, required=True, help="enrich_df.csv path (with sample_post/site/codon/wildtype/aa_mutation).")
    pc.add_argument("--out-dir", type=Path, default=Path("."), help="Where to write output CSVs.")
    pc.add_argument("--lai-effects-csv", type=Path, default=Path("LAI_GGtoAGmutation_effects.csv"))
    pc.add_argument("--s1203-effects-csv", type=Path, default=Path("1203_GGtoAGmutation_effects.csv"))
    pc.add_argument("--workers", type=int, default=22)

    # analyze
    pa = sub.add_parser("analyze", help="Generate figures and summaries from results + variant proportions.")
    pa.add_argument("--results-dir", type=Path, default=Path("."), help="Folder with LAI/1203_hypermutation_results.csv (from calculate).")
    pa.add_argument("--varprops-dir", type=Path, default=Path("."), help="Folder with LAIvariantproportions_with_wildtype.csv and 1203variantproportions_with_wildtype.csv.")
    pa.add_argument("--out-dir", type=Path, default=Path("figures"), help="Folder to write figures/CSVs.")

    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(message)s")

    if args.cmd == "calculate":
        stage_calculate(args.fasta_dir, args.enrich_csv, args.out_dir, args.lai_effects_csv, args.s1203_effects_csv, workers=args.workers)
    elif args.cmd == "analyze":
        stage_analyze(args.results_dir, args.varprops_dir, args.out_dir)
    else:
        raise SystemExit("Unknown command")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
