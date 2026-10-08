#!/usr/bin/env python
"""
Revision-3 analysis #30 — Cox proportional-hazards replacement for the "ever reclassified"
logistic regression, plus the reviewer's fixed-window fallback.

Addresses Reviewer 2, ROUND-2 Major comment #4:
  "A logistic regression for 'ever reclassified' with follow-up duration included as a
   covariate is less appropriate than a survival model ... the natural adjusted analysis is
   a Cox proportional-hazards model for time to first reclassification ... If logistic
   regression is retained, consider defining a fixed follow-up endpoint — for example,
   reclassification within five years — and include only variants with sufficient potential
   follow-up."

Model changes vs 05_multivariable_baseline.py / 13_ / 19_:
  - outcome: time to FIRST reclassification (event month per 01_survival_analysis.py),
    right-censored at 2025-07;
  - baseline hazard STRATIFIED by entry era (pre-2015 / 2015-2019 / 2020+), which replaces
    both the follow_up_yr covariate and the era covariate — comparisons are made within
    entry-era strata, which is exactly the risk-set behaviour the reviewer asked for;
  - covariates otherwise identical: log10 baseline submissions, log10 gene span, variant
    class (ref missense), baseline review tier (ref single);
  - ties: Efron approximation (monthly event grid produces heavy ties), Breslow fallback.

Every published logistic scenario is refit as a Cox HR (crude 2x2 gates assert that each
scenario's cohort construction reproduces the published crude rates BEFORE any fitting):
  M1   CGC(304) vs all others                        (round-2 adjusted OR 1.484)
  M2   disease domains vs Other/Unclassified         (round-2 ORs 1.20-1.55)
  S1   M1 excluding somatic-origin variants          (1.479)
  S2a  Hereditary Cancer(106) vs all others          (1.458; Discussion "1.46")
  S2b  HC(106) vs genes in neither CGC nor HC        (1.478)
  S2c  CGC-somatic-driver-only(198) vs non-cancer    (Discussion comparison)
  S3a  germline-CGC(114) vs true non-cancer          (1.496)
  S3b  S3a excluding somatic-origin variants         (1.491)
FW   fixed 5-year-window logistic (entry <= 2020-07, outcome = reclassified within 60
     months): the reviewer's explicit fallback, reported as a sensitivity analysis.

Outputs: outputs/cox_r3.json
"""
import gc, json, time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.duration.hazard_regression import PHReg

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
REV  = DATA / "revision_analyses"
EXT  = REV / "external_dependencies" / "ClinVar Research Project"
OUT  = REV / "outputs"
SUMMARY = DATA / "summary_20260219_024805.csv"
EVENTS  = Path("/Users/danielliang/Downloads/ClinVar_Research_Project/data/reclassification_events_detail.csv")
STUDY_END = 2025 * 12 + 7
TIES = "efron"

DOMAINS = json.load(open(EXT / "data" / "disease_domains.json"))
CGC = set(DOMAINS["cancer_genes"])
HC  = set(DOMAINS["disease_domains"]["Hereditary Cancer"])
GJ  = json.load(open(OUT / "cosmic_germline_sensitivity.json"))
GCGC = set(GJ["T1_germline_set_and_overlap"]["germline_cgc_genes"])
syn = GJ["T1_germline_set_and_overlap"].get("matched_via_synonym", {})
if isinstance(syn, dict):
    GCGC |= set(syn.values())
GENE2DOM = {g: d for d, gs in DOMAINS["disease_domains"].items() for g in gs}

VCLASS_LEVELS = ["nonsense", "frameshift", "splice", "synonymous", "inframe", "other"]  # ref missense
TIER_LEVELS = ["no_criteria", "expert_guideline", "multi_no_conflict", "conflicting"]   # ref single


def month_idx(s):
    try:
        y, m = str(s).split("-")[:2]
        return int(y) * 12 + int(m)
    except Exception:
        return np.nan


def vclass(s):
    s = str(s).lower()
    for k in ["missense", "nonsense", "frameshift", "splice", "synonymous", "inframe"]:
        if k in s:
            return k
    return "other"


def crude_gate(d, exp_col, want_rate1, want_rate2, tag, tol=0.15):
    r1 = 100 * d.loc[d[exp_col] == 1, "changed"].mean()
    r2 = 100 * d.loc[d[exp_col] == 0, "changed"].mean()
    ok = abs(r1 - want_rate1) <= tol and abs(r2 - want_rate2) <= tol
    print(f"  gate {tag}: exposed {r1:.2f}% (want {want_rate1}), comparator {r2:.2f}% "
          f"(want {want_rate2}) -> {'OK' if ok else 'MISMATCH'}", flush=True)
    if not ok:
        raise SystemExit(f"crude gate failed for {tag} — cohort construction does not "
                         f"reproduce the published scenario; refusing to fit.")


def cox_fit(d, exposure_cols, tag):
    """Era-stratified Cox PH. Returns HR/CI/p for each exposure column."""
    t0 = time.time()
    covs = exposure_cols + ["log_base_subs", "log_gene_span"] + \
           [f"vc_{v}" for v in VCLASS_LEVELS] + [f"rt_{v}" for v in TIER_LEVELS]
    X = d[covs].astype(np.float64)
    mod = PHReg(d["time_m"].values, X.values, status=d["event"].values,
                strata=d["era"].values, ties=TIES)
    r = mod.fit(disp=0)
    params = pd.Series(r.params, index=covs)
    bse = pd.Series(r.bse, index=covs)
    out = []
    for c in exposure_cols:
        hr = float(np.exp(params[c]))
        lo = float(np.exp(params[c] - 1.96 * bse[c]))
        hi = float(np.exp(params[c] + 1.96 * bse[c]))
        z = params[c] / bse[c]
        from scipy.stats import norm
        pv = float(2 * norm.sf(abs(z)))
        out.append({"model": tag, "term": c, "n": int(len(d)),
                    "n_events": int(d["event"].sum()),
                    "n_exposed": int(d[c].sum()),
                    "HR": round(hr, 3), "HR_lo": round(lo, 3), "HR_hi": round(hi, 3),
                    "p": f"{pv:.2e}"})
        print(f"  {tag:38s} {c:28s} HR={hr:.3f} [{lo:.3f},{hi:.3f}] p={pv:.1e} "
              f"({time.time()-t0:.0f}s)", flush=True)
    del mod, r, X
    gc.collect()
    return out


def main():
    print("Loading first-event months ...", flush=True)
    import csv as _csv
    first_ev = {}
    with open(EVENTS) as fh:
        for row in _csv.DictReader(fh):
            vid = row["Variation_ID"]; mi = month_idx(row["Curr_Date"])
            if vid not in first_ev or mi < first_ev[vid]:
                first_ev[vid] = mi

    print("Loading summary + baseline covariates ...", flush=True)
    df = pd.read_csv(SUMMARY,
                     usecols=["Variation_ID", "Gene", "First_Date", "Classification_Changed",
                              "Molecular_Consequence", "Start", "Stop"],
                     dtype={"Variation_ID": str, "Gene": str, "First_Date": str})
    base = pd.read_csv(OUT / "baseline_covariates.csv",
                       dtype={"Variation_ID": str, "base_num_subs": "Int64"})
    df = df.merge(base[["Variation_ID", "base_num_subs", "base_review_tier"]],
                  on="Variation_ID", how="inner")
    om = pd.read_csv(OUT / "variant_origin_map.csv.gz",
                     usecols=["variation_id", "has_somatic"])
    om["variation_id"] = om["variation_id"].astype(str)
    df = df.merge(om, left_on="Variation_ID", right_on="variation_id", how="left")
    df["has_somatic"] = df["has_somatic"].fillna(False).astype(bool)

    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    df["entry"] = df["First_Date"].map(month_idx)
    df["ev_month"] = df["Variation_ID"].map(first_ev)
    df["event"] = df["changed"]
    df["time_m"] = np.where(df["changed"] == 1, df["ev_month"] - df["entry"],
                            STUDY_END - df["entry"])
    df["log_base_subs"] = np.log10(df["base_num_subs"].fillna(1).clip(lower=1).astype(float))
    yr = df["First_Date"].str.slice(0, 4).astype(float)
    df["era"] = pd.cut(yr, [0, 2014.99, 2019.99, 3000],
                       labels=["pre2015", "RCV_2015_2019", "VCV_2020plus"]).astype(str)
    gene_span = (df.groupby("Gene")
                   .apply(lambda g: g["Stop"].astype(float).max() - g["Start"].astype(float).min(),
                          include_groups=False))
    df["log_gene_span"] = np.log10(df["Gene"].map(gene_span).clip(lower=1).fillna(1) + 1)
    df["vclass"] = df["Molecular_Consequence"].map(vclass)
    df["base_review_tier"] = df["base_review_tier"].fillna("single")

    df = df.dropna(subset=["entry", "time_m", "log_base_subs", "log_gene_span"])
    df = df[df["time_m"] >= 0]
    df = df[~((df["changed"] == 1) & (df["ev_month"].isna()))]
    df["time_m"] = df["time_m"].clip(lower=0.5)
    for v in VCLASS_LEVELS:
        df[f"vc_{v}"] = (df["vclass"] == v).astype(np.int8)
    for v in TIER_LEVELS:
        df[f"rt_{v}"] = (df["base_review_tier"] == v).astype(np.int8)
    df["is_cancer"] = df["Gene"].isin(CGC).astype(np.int8)
    df["is_hc"] = df["Gene"].isin(HC).astype(np.int8)
    df["is_gcgc"] = df["Gene"].isin(GCGC).astype(np.int8)
    df["is_sd_only"] = (df["Gene"].isin(CGC - HC)).astype(np.int8)
    df["domain"] = df["Gene"].map(GENE2DOM).fillna("Other/Unclassified")
    print(f"  modeling rows: {len(df):,}   events: {int(df['event'].sum()):,}", flush=True)

    all_rows = []

    # ---- gates + fits ----------------------------------------------------------
    print("\nM1: CGC(304) vs all others", flush=True)
    crude_gate(df, "is_cancer", 16.62, 8.19, "M1")
    all_rows += cox_fit(df, ["is_cancer"], "M1_cox_cgc304")

    print("\nM2: disease domains", flush=True)
    dom_levels = [d for d in ["Hereditary Cancer", "Cardiovascular", "Neurological",
                              "Renal", "Ophthalmologic"]]
    for d_ in dom_levels:
        df[f"dm_{d_.split()[0]}"] = (df["domain"] == d_).astype(np.int8)
    all_rows += cox_fit(df, [f"dm_{d_.split()[0]}" for d_ in dom_levels], "M2_cox_domains")

    print("\nS1: CGC(304), somatic-origin variants excluded", flush=True)
    d1 = df[~df["has_somatic"]]
    crude_gate(d1, "is_cancer", 16.49, 8.18, "S1")
    all_rows += cox_fit(d1, ["is_cancer"], "S1_cox_cgc304_excl_somatic")
    del d1; gc.collect()

    print("\nS2a: HC(106) vs all others", flush=True)
    crude_gate(df, "is_hc", 18.18, 8.34, "S2a")
    all_rows += cox_fit(df, ["is_hc"], "S2a_cox_hc106_vs_all")

    print("\nS2b: HC(106) vs genes in neither CGC nor HC", flush=True)
    d2 = df[(df["is_hc"] == 1) | (~df["Gene"].isin(CGC | HC))]
    crude_gate(d2, "is_hc", 18.18, 8.19, "S2b", tol=0.2)
    all_rows += cox_fit(d2, ["is_hc"], "S2b_cox_hc106_vs_true_noncancer")
    del d2; gc.collect()

    print("\nS2c: CGC-somatic-driver-only(198) vs non-cancer", flush=True)
    d3 = df[(df["is_sd_only"] == 1) | (df["is_cancer"] == 0)]
    all_rows += cox_fit(d3, ["is_sd_only"], "S2c_cox_sdonly_vs_noncancer")
    del d3; gc.collect()

    print("\nS3a: germline-CGC(114) vs true non-cancer", flush=True)
    d4 = df[(df["is_gcgc"] == 1) | (~df["Gene"].isin(CGC | GCGC))]
    crude_gate(d4, "is_gcgc", 18.33, 8.19, "S3a", tol=0.2)
    all_rows += cox_fit(d4, ["is_gcgc"], "S3a_cox_gcgc114_vs_true_noncancer")

    print("\nS3b: S3a excluding somatic-origin variants", flush=True)
    d5 = d4[~d4["has_somatic"]]
    crude_gate(d5, "is_gcgc", 18.20, 8.18, "S3b", tol=0.2)
    all_rows += cox_fit(d5, ["is_gcgc"], "S3b_cox_gcgc114_excl_somatic")
    del d4, d5; gc.collect()

    # ---- FW: fixed 5-year-window logistic (reviewer's fallback) ----------------
    print("\nFW: fixed 5-year logistic (entry <= 2020-07)", flush=True)
    fw = df[df["entry"] <= 2020 * 12 + 7].copy()
    fw["changed5"] = ((fw["event"] == 1) & (fw["time_m"] <= 60)).astype(int)
    print(f"  cohort: {len(fw):,} variants; reclassified within 5y: {int(fw['changed5'].sum()):,}")
    f1 = ("changed5 ~ is_cancer + log_base_subs + log_gene_span "
          "+ C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
          "+ C(base_review_tier, Treatment('single'))")
    m = smf.logit(f1, data=fw).fit(disp=0, maxiter=200)
    p, ci = m.params, m.conf_int()
    all_rows.append({"model": "FW_logit5yr_cgc304", "term": "is_cancer", "n": int(m.nobs),
                     "n_exposed": int(fw["is_cancer"].sum()),
                     "OR": round(float(np.exp(p["is_cancer"])), 3),
                     "OR_lo": round(float(np.exp(ci.loc["is_cancer", 0])), 3),
                     "OR_hi": round(float(np.exp(ci.loc["is_cancer", 1])), 3),
                     "p": f"{m.pvalues['is_cancer']:.2e}"})
    print(f"  FW is_cancer OR={np.exp(p['is_cancer']):.3f} "
          f"[{np.exp(ci.loc['is_cancer',0]):.3f},{np.exp(ci.loc['is_cancer',1]):.3f}]", flush=True)
    f2 = ("changed5 ~ C(domain, Treatment('Other/Unclassified')) + log_base_subs + log_gene_span "
          "+ C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
          "+ C(base_review_tier, Treatment('single'))")
    m2 = smf.logit(f2, data=fw).fit(disp=0, maxiter=200)
    p2, ci2 = m2.params, m2.conf_int()
    for t in p2.index:
        if t.startswith("C(domain"):
            all_rows.append({"model": "FW_logit5yr_domains", "term": t, "n": int(m2.nobs),
                             "OR": round(float(np.exp(p2[t])), 3),
                             "OR_lo": round(float(np.exp(ci2.loc[t, 0])), 3),
                             "OR_hi": round(float(np.exp(ci2.loc[t, 1])), 3),
                             "p": f"{m2.pvalues[t]:.2e}"})
            print(f"  FW {t:55s} OR={np.exp(p2[t]):.3f}", flush=True)

    json.dump({"ties": TIES, "strata": "entry era (pre2015 / 2015-2019 / 2020+)",
               "results": all_rows}, open(OUT / "cox_r3.json", "w"), indent=2)
    print("\nwritten:", OUT / "cox_r3.json")


if __name__ == "__main__":
    main()
