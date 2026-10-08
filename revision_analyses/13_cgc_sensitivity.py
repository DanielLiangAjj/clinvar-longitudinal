#!/usr/bin/env python
"""
Revision analysis #13 — COSMIC CGC composition audit + gene-set sensitivity (T4)
and adjusted ORs under each gene-set x origin definition (T5).

Co-author's objection: COSMIC's Cancer Gene Census is largely a SOMATIC-driver list, so
it is the wrong exposure for a GERMLINE reclassification comparison. This script
(a) quantifies how much of the 304-gene CGC set sits outside the curated
'Hereditary Cancer' domain, (b) re-runs the crude contrast with Hereditary Cancer as the
exposure, and (c) re-fits the adjusted model of 05_multivariable_baseline.py under every
gene-set x origin-filter combination.

Methodology mirrored EXACTLY from:
  - two_by_two()            -> 02_statistics.py:31-49
  - baseline covariates     -> 05_multivariable_baseline.py:48-78
      log_base_subs (05:59), follow_up_yr (05:60-62), era (05:63-64),
      log_gene_span (05:65-67), vclass (05:28-32,68), base_review_tier (from
      outputs/baseline_covariates.csv, built by 04_history_pass.py)
  - model formula M1       -> 05_multivariable_baseline.py:82-84
NOTE: baseline covariates come from outputs/baseline_covariates.csv, which already
exists; the 14GB classification_history CSV is therefore NOT needed here.

Outputs: outputs/cgc_sensitivity.json
"""
import gc, json, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

DROPBOX = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/outputs")
SUMMARY = DROPBOX / "summary_20260219_024805.csv"
ORIGIN_MAP = OUT / "variant_origin_map.csv.gz"


def two_by_two(x1, n1, x2, n2, label1, label2):
    """Verbatim port of 02_statistics.py:31-49."""
    a, b = x1, n1 - x1
    c, d = x2, n2 - x2
    p1, p2 = a / n1, c / n2
    rr = p1 / p2
    se_lnrr = np.sqrt(b / (a * n1) + d / (c * n2))
    rr_lo, rr_hi = rr * np.exp(-1.96 * se_lnrr), rr * np.exp(1.96 * se_lnrr)
    orr = (a * d) / (b * c)
    se_lnor = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    or_lo, or_hi = orr * np.exp(-1.96 * se_lnor), orr * np.exp(1.96 * se_lnor)
    chi2, p, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=False)
    return {
        "group1": label1, "n1": int(n1), "changed1": int(a), "rate1_pct": round(100 * p1, 3),
        "group2": label2, "n2": int(n2), "changed2": int(c), "rate2_pct": round(100 * p2, 3),
        "risk_ratio": round(rr, 3), "rr_95ci": [round(rr_lo, 3), round(rr_hi, 3)],
        "odds_ratio": round(orr, 3), "or_95ci": [round(or_lo, 3), round(or_hi, 3)],
        "chi2": round(chi2, 1), "p_value": (f"{p:.2e}" if p > 0 else "<1e-300"),
    }


def vclass(s):
    """05_multivariable_baseline.py:28-32"""
    s = str(s).lower()
    for k in ["missense", "nonsense", "frameshift", "splice", "synonymous", "inframe"]:
        if k in s:
            return k
    return "other"


def fit_adjusted(d, exposure, tag):
    """Re-fit M1 from 05_multivariable_baseline.py:82-84 with `exposure` in place of
    is_cancer. Returns the exposure's adjusted OR + 95% CI."""
    f = (f"changed ~ {exposure} + log_base_subs + follow_up_yr + log_gene_span "
         "+ C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
         "+ C(base_review_tier, Treatment('single'))")
    m = smf.logit(f, data=d).fit(disp=0, maxiter=200)
    p, ci, pv = m.params, m.conf_int(), m.pvalues
    out = {
        "model": tag, "term": exposure, "n": int(m.nobs),
        "n_exposed": int(d[exposure].sum()),
        "pseudo_r2": round(float(m.prsquared), 4),
        "OR": round(float(np.exp(p[exposure])), 3),
        "OR_lo": round(float(np.exp(ci.loc[exposure, 0])), 3),
        "OR_hi": round(float(np.exp(ci.loc[exposure, 1])), 3),
        "p": f"{pv[exposure]:.2e}",
    }
    print(f"  {tag:44s} n={out['n']:>9,} exposed={out['n_exposed']:>8,} "
          f"OR={out['OR']:.3f} [{out['OR_lo']:.3f}, {out['OR_hi']:.3f}] p={out['p']}",
          flush=True)
    del m
    gc.collect()
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOMAINS = json.load(open(PROJ / "disease_domains.json"))
    cgc = set(DOMAINS["cancer_genes"])
    hc = set(DOMAINS["disease_domains"]["Hereditary Cancer"])

    # ---------------- T4a: CGC composition ----------------
    comp = {
        "n_cgc": len(cgc),
        "n_hereditary_cancer": len(hc),
        "n_cgc_and_hc": len(cgc & hc),
        "n_cgc_not_hc": len(cgc - hc),
        "n_hc_not_cgc": len(hc - cgc),
        "pct_of_cgc_outside_hereditary_cancer": round(100 * len(cgc - hc) / len(cgc), 1),
        "hereditary_cancer_is_subset_of_cgc": hc.issubset(cgc),
        "examples_cgc_not_hc": sorted(cgc - hc)[:25],
    }
    print("=== T4a: CGC COMPOSITION ===")
    for k, v in comp.items():
        print(f"  {k:44s} {v}")

    # ---------------- data ----------------
    print("\nloading origin map + summary ...", flush=True)
    om = pd.read_csv(ORIGIN_MAP, dtype={"variation_id": str}, keep_default_na=False)
    for c in ("has_somatic", "has_origin", "strict_germline"):
        if om[c].dtype == object:
            om[c] = om[c].map({"True": True, "False": False}).astype(bool)

    # crude analysis uses the full summary (mirrors 02_statistics.py:53-64)
    crude = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "Classification_Changed"],
                        dtype={"Variation_ID": str, "Gene": str})
    crude["changed"] = crude["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    crude = crude.merge(om[["variation_id", "has_somatic"]], left_on="Variation_ID",
                        right_on="variation_id", how="left")
    crude["has_somatic"] = crude["has_somatic"].fillna(False).astype(bool)
    crude["is_cgc"] = crude["Gene"].isin(cgc)
    crude["is_hc"] = crude["Gene"].isin(hc)
    crude["is_cgc_not_hc"] = crude["is_cgc"] & ~crude["is_hc"]
    print(f"  crude rows: {len(crude):,}")

    def grp(d, mask):
        s = d[mask]
        return int(s["changed"].sum()), int(len(s))

    # ---------------- T4b: crude contrasts ----------------
    crude_res = {}
    # (i) all variants, HC as exposure (straight substitution: comparator keeps the
    #     198 somatic-driver CGC genes, per the literal spec)
    x1, n1 = grp(crude, crude["is_hc"]); x2, n2 = grp(crude, ~crude["is_hc"])
    crude_res["hc_vs_all_others_all_variants"] = two_by_two(
        x1, n1, x2, n2, "Hereditary Cancer (106)", "All other genes")

    # (ii) HC exposure, excluding somatic-origin variants
    m = ~crude["has_somatic"]
    x1, n1 = grp(crude, crude["is_hc"] & m); x2, n2 = grp(crude, ~crude["is_hc"] & m)
    crude_res["hc_vs_all_others_exclude_somatic"] = two_by_two(
        x1, n1, x2, n2, "Hereditary Cancer, somatic excluded", "All other genes, somatic excluded")

    # (iii) supplementary — drop the 198 CGC-not-HC genes from the comparator entirely
    keep = ~crude["is_cgc_not_hc"]
    x1, n1 = grp(crude, crude["is_hc"] & keep); x2, n2 = grp(crude, ~crude["is_cgc"] & keep)
    crude_res["hc_vs_noncancer_drop_somatic_drivers"] = two_by_two(
        x1, n1, x2, n2, "Hereditary Cancer (106)", "Genes in neither CGC nor HC")

    keep2 = keep & ~crude["has_somatic"]
    x1, n1 = grp(crude, crude["is_hc"] & keep2); x2, n2 = grp(crude, ~crude["is_cgc"] & keep2)
    crude_res["hc_vs_noncancer_drop_drivers_exclude_somatic"] = two_by_two(
        x1, n1, x2, n2, "Hereditary Cancer, somatic excl.", "Neither CGC nor HC, somatic excl.")

    # (iv) the decomposition that answers the co-author directly:
    #      the 198 somatic-driver-only CGC genes vs non-cancer
    x1, n1 = grp(crude, crude["is_cgc_not_hc"]); x2, n2 = grp(crude, ~crude["is_cgc"])
    crude_res["cgc_not_hc_vs_noncancer"] = two_by_two(
        x1, n1, x2, n2, "CGC somatic-driver-only (198)", "Non-cancer (not in CGC)")

    x1, n1 = grp(crude, crude["is_hc"]); x2, n2 = grp(crude, ~crude["is_cgc"])
    crude_res["hc_vs_noncancer"] = two_by_two(
        x1, n1, x2, n2, "Hereditary Cancer (106)", "Non-cancer (not in CGC)")

    print("\n=== T4b: CRUDE CONTRASTS ===")
    print(f"{'contrast':46s} {'n1':>9} {'rate1':>7} {'n2':>10} {'rate2':>7} {'RR':>6} {'OR':>6}")
    for k, v in crude_res.items():
        print(f"{k:46s} {v['n1']:>9,} {v['rate1_pct']:>6.2f}% {v['n2']:>10,} "
              f"{v['rate2_pct']:>6.2f}% {v['risk_ratio']:>6.3f} {v['odds_ratio']:>6.3f}")

    del crude
    gc.collect()

    # ---------------- T5: adjusted models ----------------
    print("\n=== T5: ADJUSTED MODELS (mirrors 05_multivariable_baseline.py M1) ===", flush=True)
    base = pd.read_csv(OUT / "baseline_covariates.csv",
                       dtype={"Variation_ID": str, "base_num_subs": "Int64"})
    print(f"  baseline rows: {len(base):,}")
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "First_Date",
                                       "Classification_Changed", "Molecular_Consequence",
                                       "Start", "Stop"],
                     dtype={"Variation_ID": str, "Gene": str, "First_Date": str})
    df = df.merge(base[["Variation_ID", "base_num_subs", "base_review_tier"]],
                  on="Variation_ID", how="inner")
    print(f"  merged rows: {len(df):,}")
    del base
    gc.collect()

    # --- covariates, verbatim from 05_multivariable_baseline.py:58-68 ---
    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    df["log_base_subs"] = np.log10(df["base_num_subs"].fillna(1).clip(lower=1).astype(float))
    yr = df["First_Date"].str.slice(0, 4).astype(float)
    mo = df["First_Date"].str.slice(5, 7).astype(float)
    df["follow_up_yr"] = (2025 + 7 / 12) - (yr + mo / 12)
    df["era"] = pd.cut(yr, [0, 2014.99, 2019.99, 3000],
                       labels=["pre2015", "RCV_2015_2019", "VCV_2020plus"])
    gene_span = (df.groupby("Gene").apply(
        lambda g: g["Stop"].astype(float).max() - g["Start"].astype(float).min(),
        include_groups=False))
    df["log_gene_span"] = np.log10(df["Gene"].map(gene_span).clip(lower=1).fillna(1) + 1)
    df["vclass"] = df["Molecular_Consequence"].map(vclass)

    df["is_cancer"] = df["Gene"].isin(cgc).astype(int)      # 05:72
    df["is_hc"] = df["Gene"].isin(hc).astype(int)
    df["is_cgc_not_hc"] = (df["Gene"].isin(cgc) & ~df["Gene"].isin(hc)).astype(int)

    df = df.merge(om[["variation_id", "has_somatic", "has_origin", "strict_germline"]],
                  left_on="Variation_ID", right_on="variation_id", how="left")
    for c in ("has_somatic", "has_origin", "strict_germline"):
        df[c] = df[c].fillna(False).astype(bool)

    keep = ["changed", "log_base_subs", "follow_up_yr", "log_gene_span", "vclass",
            "base_review_tier", "is_cancer", "is_hc", "is_cgc_not_hc", "era",
            "has_somatic", "has_origin", "strict_germline"]
    df = df[keep].dropna()
    print(f"  modeling rows: {len(df):,}", flush=True)

    adj = []
    # (a) baseline — must reproduce OR = 1.484 [1.470, 1.499]
    adj.append(fit_adjusted(df, "is_cancer", "a_baseline_CGC_all_variants"))
    published_or = 1.484
    got = adj[0]["OR"]
    print(f"    -> baseline gate: got OR={got}, published OR={published_or} "
          f"({'MATCH' if abs(got - published_or) < 0.002 else 'MISMATCH'})")

    # (b) CGC exposure, excluding somatic-origin variants
    adj.append(fit_adjusted(df[~df["has_somatic"]], "is_cancer",
                            "b_CGC_exclude_somatic"))
    # (c) Hereditary-Cancer exposure, all variants
    adj.append(fit_adjusted(df, "is_hc", "c_HC_all_variants"))
    # (d) Hereditary-Cancer exposure, excluding somatic-origin variants
    adj.append(fit_adjusted(df[~df["has_somatic"]], "is_hc",
                            "d_HC_exclude_somatic"))
    # supplementary: HC vs true non-cancer (drop the 198 driver-only genes)
    adj.append(fit_adjusted(df[df["is_cgc_not_hc"] == 0], "is_hc",
                            "e_HC_vs_noncancer_drop_drivers"))
    adj.append(fit_adjusted(df[(df["is_cgc_not_hc"] == 0) & ~df["has_somatic"]], "is_hc",
                            "f_HC_vs_noncancer_drop_drivers_excl_somatic"))
    # supplementary: the driver-only genes as their own exposure
    adj.append(fit_adjusted(df[df["is_hc"] == 0], "is_cgc_not_hc",
                            "g_CGCdriver_only_vs_noncancer"))

    res = {"cgc_composition": comp, "crude": crude_res, "adjusted": adj,
           "notes": {
               "baseline_covariates_source": str(OUT / "baseline_covariates.csv"),
               "classification_history_needed": False,
               "covariates_built": ["log_base_subs", "follow_up_yr", "log_gene_span",
                                    "vclass", "base_review_tier", "era"],
               "covariates_dropped": [],
               "model_formula": "05_multivariable_baseline.py:82-84 (M1), exposure swapped",
           }}
    json.dump(res, open(OUT / "cgc_sensitivity.json", "w"), indent=2)
    print(f"\nwritten: {OUT/'cgc_sensitivity.json'}")


if __name__ == "__main__":
    main()
