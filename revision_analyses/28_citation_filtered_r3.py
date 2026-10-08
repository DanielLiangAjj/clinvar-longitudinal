#!/usr/bin/env python
"""
Revision-3 analysis #28 — recompute ALL citation statistics and figures with the
variant-specific (methodology-filtered) counting method.

Addresses Reviewer 2, ROUND-2 Major comment #6:
  The round-1 revision switched the TEXT to filtered citation counts (SHERLOC 28492532 and
  ACMG/AMP-2015 25741868 excluded; mean 1.58 -> 0.50 among cited variants; ~12% of all
  variants retain a genuine citation) but Figures 9g/9h, 11b and 12a/12b were still rendered
  from RAW Num_PMIDs. The reviewer correctly read the resulting contradiction off Figure 9g
  (raw domain means 0.80-1.66, medians 1) versus the filtered text (0.50; ~12%).

This script makes the figures and text consistent by computing per-variant FILTERED counts
from the per-variant PMID lists and regenerating every citation display:
  Figure 9  (module renderer figure10_disease_domain):
     panel g -> stacked composition of variant-specific citations per variant by domain
                (median is 0 in every domain and the panel now says so honestly)
     panel h -> adds a third bar: mean variant-specific refs alongside raw refs/submissions
  Figure 11 (module renderer figure4): panel b -> paired raw vs filtered distribution
  Figure 12 (module renderer figure7): panels a,b -> paired raw vs filtered series
All new figure files go to figures/r3/ — nothing from round 2 is overwritten.

Text numbers for the manuscript edit are written to outputs/citation_r3_numbers.json.

Data: variant_pmids_latest.csv (per-variant PMID lists, latest snapshot) + summary CSV
      (Gene, First_Date, Num_PMIDs for all 3,678,303 variants). Variants with no PMID row
      have a filtered count of 0.
"""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
REV  = DATA / "revision_analyses"
EXT  = REV / "external_dependencies" / "ClinVar Research Project"
OUT  = REV / "outputs"
FIGR3 = REV / "figures" / "r3"
FIGR3.mkdir(parents=True, exist_ok=True)

SUMMARY = DATA / "summary_20260219_024805.csv"
PMIDS   = OUT / "variant_pmids_latest.csv"
CACHE_AGG = Path("/Users/danielliang/Downloads/ClinVar_Research_Project/figures/.cache/summary_agg.json")

# The two assertion-method papers named in the manuscript text (primary filter).
EXCL = {"28492532", "25741868"}   # SHERLOC (Nykamp 2017), ACMG/AMP 2015 (Richards)

sys.path.insert(0, str(EXT / "scripts"))
import clinvar_visualizations as V  # noqa: E402  (edited copy with filtered-aware branches)

DOMAINS = json.load(open(EXT / "data" / "disease_domains.json"))
CANCER = set(DOMAINS["cancer_genes"])
COMMON = set(DOMAINS["common_disease_genes"])
GENE2DOM = {g: d for d, gs in DOMAINS["disease_domains"].items() for g in gs}
DOMAIN_ORDER = ["Hereditary Cancer", "Cardiovascular", "Neurological",
                "Renal", "Ophthalmologic", "Other/Unclassified"]


def cls4(n):
    if n == 0: return "0"
    if n == 1: return "1"
    if n <= 5: return "2-5"
    return ">=6"


def bin6(n):
    if n == 0: return "0"
    if n == 1: return "1"
    if n <= 5: return "2-5"
    if n <= 10: return "6-10"
    if n <= 20: return "11-20"
    return ">20"


def main():
    # ---- 1. per-variant filtered counts from the PMID lists --------------------
    print("Streaming per-variant PMID lists ...", flush=True)
    filt, raw_file = {}, {}
    with open(PMIDS) as fh:
        for row in csv.DictReader(fh):
            s = set(p for p in row["PMIDs"].split(";") if p)
            filt[row["Variation_ID"]] = len(s - EXCL)
            raw_file[row["Variation_ID"]] = len(s)
    n_cited = len(filt)
    print(f"  variants with >=1 listed PMID: {n_cited:,}")

    # ---- 2. summary join over ALL variants ------------------------------------
    print("Loading summary ...", flush=True)
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "First_Date", "Num_PMIDs"],
                     dtype={"Variation_ID": str, "Gene": str, "First_Date": str})
    df["Num_PMIDs"] = pd.to_numeric(df["Num_PMIDs"], errors="coerce").fillna(0).astype(int)
    df["filt"] = df["Variation_ID"].map(filt).fillna(0).astype(int)
    df["domain"] = df["Gene"].map(GENE2DOM).fillna("Other/Unclassified")
    total = len(df)
    print(f"  total variants: {total:,}")

    res = {"exclusion_pmids": sorted(EXCL), "n_variants_total": total,
           "n_variants_with_any_listed_pmid_summary": int((df["Num_PMIDs"] >= 1).sum()),
           "n_variants_with_any_listed_pmid_extraction": n_cited}

    # headline text numbers
    cited_mask = df["Variation_ID"].isin(filt)          # the extraction set (n=2,249,743)
    res["mean_raw_over_cited"] = round(float(np.mean([raw_file[v] for v in raw_file])), 3)
    res["mean_filtered_over_cited"] = round(float(df.loc[cited_mask, "filt"].mean()), 3)
    res["mean_filtered_over_all"] = round(float(df["filt"].mean()), 3)
    n_ge1_filt = int((df["filt"] >= 1).sum())
    res["n_ge1_filtered"] = n_ge1_filt
    res["pct_ge1_filtered_of_all"] = round(100 * n_ge1_filt / total, 2)
    res["pct_ge1_raw_of_all_summary"] = round(100 * res["n_variants_with_any_listed_pmid_summary"] / total, 2)
    res["total_refs_raw_summary"] = int(df["Num_PMIDs"].sum())
    res["total_refs_filtered"] = int(df["filt"].sum())

    # ---- 3. figure-11 panel b bins + figure-12 sums ----------------------------
    pmid_bins_filtered = df["filt"].map(bin6).value_counts().to_dict()
    fsum_by_fd = df.groupby("First_Date")["filt"].sum().to_dict()

    # ---- 4. per-domain composition + subgroup means ----------------------------
    by_domain = {}
    for d in DOMAIN_ORDER:
        sub = df[df["domain"] == d]
        if not len(sub):
            continue
        comp = sub["filt"].map(cls4).value_counts(normalize=True).mul(100).round(2).to_dict()
        cited = sub.loc[sub["filt"] >= 1, "filt"]
        by_domain[d] = {
            "n": int(len(sub)),
            "cit_classes": {k: float(comp.get(k, 0.0)) for k in ["0", "1", "2-5", ">=6"]},
            "pct_ge1_filtered": round(100 * float((sub["filt"] >= 1).mean()), 1),
            "mean_filtered_all": round(float(sub["filt"].mean()), 3),
            "mean_filtered_cited": round(float(cited.mean()), 3) if len(cited) else 0.0,
            "median_filtered_cited": float(cited.median()) if len(cited) else 0.0,
            "q25_cited": float(cited.quantile(.25)) if len(cited) else 0.0,
            "q75_cited": float(cited.quantile(.75)) if len(cited) else 0.0,
            "mean_raw_all": round(float(sub["Num_PMIDs"].mean()), 3),
        }
    res["by_domain"] = by_domain

    groups = {
        "cancer":    df["Gene"].isin(CANCER),
        "noncancer": ~df["Gene"].isin(CANCER),
        "common":    df["Gene"].isin(COMMON),
        "rare":      ~df["Gene"].isin(COMMON),
    }
    res["subgroups"] = {}
    for name, mask in groups.items():
        sub = df[mask]
        res["subgroups"][name] = {
            "n": int(len(sub)),
            "mean_raw_all": round(float(sub["Num_PMIDs"].mean()), 3),
            "mean_filtered_all": round(float(sub["filt"].mean()), 3),
            "pct_ge1_filtered": round(100 * float((sub["filt"] >= 1).mean()), 1),
        }

    json.dump(res, open(OUT / "citation_r3_numbers.json", "w"), indent=2)
    print(json.dumps({k: v for k, v in res.items() if not isinstance(v, dict)}, indent=2))
    print("\nsubgroups (raw should reproduce panel h: 1.49/0.90/0.96/2.32):")
    for name, s in res["subgroups"].items():
        print(f"  {name:10s} n={s['n']:>9,}  raw={s['mean_raw_all']:>6.3f}  filt={s['mean_filtered_all']:>6.3f}  %>=1 filt={s['pct_ge1_filtered']}")
    print("\nby domain:")
    for d, s in by_domain.items():
        print(f"  {d:20s} raw_mean={s['mean_raw_all']:>6.3f}  filt_mean={s['mean_filtered_all']:>6.3f}  %>=1={s['pct_ge1_filtered']:>5.1f}  classes={s['cit_classes']}")

    # ---- 5. regenerate Figure 11 (module figure4) and Figure 12 (figure7) ------
    agg = json.load(open(CACHE_AGG))
    agg["pmid_bins_filtered"] = pmid_bins_filtered
    agg["pmid_sum_by_first_date_filtered"] = {k: int(v) for k, v in fsum_by_fd.items()}
    V.figure4(agg, {}, FIGR3)
    V.figure7(agg, FIGR3)

    # ---- 6. regenerate Figure 9 (module figure10_disease_domain) ---------------
    domain = json.load(open(EXT / "data" / "domain_summary.json"))
    for d, s in by_domain.items():
        if d in domain.get("by_domain", {}):
            domain["by_domain"][d]["cit_classes"] = s["cit_classes"]
    for sg in ["cancer", "noncancer"]:
        domain["cancer_vs_noncancer"][sg]["mean_pmids_filtered"] = res["subgroups"][sg]["mean_filtered_all"]
    for sg in ["rare", "common"]:
        domain["rare_vs_common"][sg]["mean_pmids_filtered"] = res["subgroups"][sg]["mean_filtered_all"]
    V.figure10_disease_domain(domain, FIGR3)

    print("\nfigures written to", FIGR3)
    for n in ["figure10_disease_domain_analysis", "figure4_evidence_linkage_completeness",
              "figure7_pubmed_citation_growth"]:
        p = FIGR3 / f"{n}.png"
        print(f"  {p.name}  {'OK ' + str(p.stat().st_size // 1024) + ' KB' if p.exists() else 'MISSING'}")


if __name__ == "__main__":
    main()
