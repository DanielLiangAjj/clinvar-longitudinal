#!/usr/bin/env python
"""
Revision analysis #5 — multivariable logistic regression with BASELINE covariates.

Fixes the collider problem in the first-pass model: R2 asked us to adjust for "number of
submitters" and "review status", but the summary CSV holds their FINAL-snapshot values, which
are consequences of reclassification (a currently-'conflicting' status is an outcome). Here we
use each variant's values AT FIRST APPEARANCE (from 04_history_pass baseline_covariates.csv),
so the model estimates the effect of disease domain / variant class adjusting for the state the
variant STARTED in — a defensible confounder set.

Outcome: Classification_Changed (ever reclassified).
Reports adjusted ORs (95% CI) for the cancer/domain effect and all covariates.
Output: outputs/stats_multivariable_baseline.csv
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
SUMMARY = DATA / "summary_20260219_024805.csv"
DOMAINS = json.load(open(PROJ / "disease_domains.json"))

def vclass(s):
    s = str(s).lower()
    for k in ["missense", "nonsense", "frameshift", "splice", "synonymous", "inframe"]:
        if k in s: return k
    return "other"

def run(df, formula, tag):
    m = smf.logit(formula, data=df).fit(disp=0, maxiter=200)
    params, conf, pv = m.params, m.conf_int(), m.pvalues
    rows = [{"model": tag, "term": t, "OR": round(float(np.exp(params[t])), 3),
             "OR_lo": round(float(np.exp(conf.loc[t, 0])), 3),
             "OR_hi": round(float(np.exp(conf.loc[t, 1])), 3),
             "p": f"{pv[t]:.2e}"} for t in params.index]
    print(f"\n=== {tag}  (n={int(m.nobs):,}, pseudo-R2={m.prsquared:.4f}) ===")
    for r in rows:
        if r["term"] != "Intercept":
            print(f"  {r['term']:52s} OR={r['OR']:>8.3f}  [{r['OR_lo']:.3f}, {r['OR_hi']:.3f}]  p={r['p']}")
    return rows

def main():
    base = pd.read_csv(OUT / "baseline_covariates.csv",
                       dtype={"Variation_ID": str, "base_num_subs": "Int64"})
    print(f"baseline rows: {len(base):,}")
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "First_Date",
                     "Classification_Changed", "Molecular_Consequence", "Start", "Stop"],
                     dtype={"Variation_ID": str, "Gene": str, "First_Date": str})
    df = df.merge(base[["Variation_ID", "base_num_subs", "base_review_tier"]],
                  on="Variation_ID", how="inner")
    print(f"merged rows: {len(df):,}")

    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    df["log_base_subs"] = np.log10(df["base_num_subs"].fillna(1).clip(lower=1).astype(float))
    yr = df["First_Date"].str.slice(0, 4).astype(float)
    mo = df["First_Date"].str.slice(5, 7).astype(float)
    df["follow_up_yr"] = (2025 + 7/12) - (yr + mo/12)
    df["era"] = pd.cut(yr, [0, 2014.99, 2019.99, 3000],
                       labels=["pre2015", "RCV_2015_2019", "VCV_2020plus"])
    gene_span = (df.groupby("Gene").apply(
        lambda g: g["Stop"].astype(float).max() - g["Start"].astype(float).min(), include_groups=False))
    df["log_gene_span"] = np.log10(df["Gene"].map(gene_span).clip(lower=1).fillna(1) + 1)
    df["vclass"] = df["Molecular_Consequence"].map(vclass)

    cancer = set(DOMAINS["cancer_genes"])
    gene2dom = {g: d for d, gs in DOMAINS["disease_domains"].items() for g in gs}
    df["is_cancer"] = df["Gene"].isin(cancer).astype(int)
    df["domain"] = df["Gene"].map(gene2dom).fillna("Other/Unclassified")

    keep = ["changed", "log_base_subs", "follow_up_yr", "log_gene_span", "vclass",
            "base_review_tier", "is_cancer", "domain", "era"]
    df = df[keep].dropna()
    print(f"modeling rows: {len(df):,}")

    all_rows = []
    # Model 1: cancer binary, adjusted for baseline confounders
    f1 = ("changed ~ is_cancer + log_base_subs + follow_up_yr + log_gene_span "
          "+ C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
          "+ C(base_review_tier, Treatment('single'))")
    all_rows += run(df, f1, "M1_cancer_adjusted_baseline")
    # Model 2: disease domain (vs Other/Unclassified), adjusted
    f2 = ("changed ~ C(domain, Treatment('Other/Unclassified')) + log_base_subs + follow_up_yr "
          "+ log_gene_span + C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
          "+ C(base_review_tier, Treatment('single'))")
    all_rows += run(df, f2, "M2_domain_adjusted_baseline")

    pd.DataFrame(all_rows).to_csv(OUT / "stats_multivariable_baseline.csv", index=False)
    print(f"\nwritten: {OUT/'stats_multivariable_baseline.csv'}")

if __name__ == "__main__":
    main()
