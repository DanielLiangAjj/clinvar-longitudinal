#!/usr/bin/env python
"""
Revision analysis #2 — Statistical testing + effect sizes + multivariable adjustment.

Addresses Reviewer 2, Major point #4:
  "Statistical testing is needed for proportion comparisons throughout... These comparisons
   require appropriate tests and effect estimates, ideally with models adjusting for
   confounders such as number of submitters, gene size, variant class, review status,
   time in ClinVar, and submission era."

Part A — Univariate: for each headline comparison, a 2x2 table of reclassified/not with
         chi-square test, risk ratio (Katz 95% CI) and odds ratio (95% CI).
Part B — Multivariable logistic regression: P(reclassified) ~ the confounders R2 named,
         reporting adjusted ORs with 95% CIs.

Outputs: outputs/stats_univariate.json, outputs/stats_multivariable.csv
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
SUMMARY = DATA / "summary_20260219_024805.csv"
DOMAINS = json.load(open(PROJ / "disease_domains.json"))

def two_by_two(x1, n1, x2, n2, label1, label2):
    """x=changed, n=total in each group. Returns RR, OR, chi2 p with CIs."""
    a, b = x1, n1 - x1        # group1: changed, unchanged
    c, d = x2, n2 - x2        # group2
    p1, p2 = a / n1, c / n2
    rr = p1 / p2
    se_lnrr = np.sqrt(b / (a * n1) + d / (c * n2))
    rr_lo, rr_hi = rr * np.exp(-1.96 * se_lnrr), rr * np.exp(1.96 * se_lnrr)
    orr = (a * d) / (b * c)
    se_lnor = np.sqrt(1/a + 1/b + 1/c + 1/d)
    or_lo, or_hi = orr * np.exp(-1.96 * se_lnor), orr * np.exp(1.96 * se_lnor)
    chi2, p, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=False)
    return {
        "group1": label1, "n1": int(n1), "changed1": int(a), "rate1_pct": round(100*p1, 3),
        "group2": label2, "n2": int(n2), "changed2": int(c), "rate2_pct": round(100*p2, 3),
        "risk_ratio": round(rr, 3), "rr_95ci": [round(rr_lo, 3), round(rr_hi, 3)],
        "odds_ratio": round(orr, 3), "or_95ci": [round(or_lo, 3), round(or_hi, 3)],
        "chi2": round(chi2, 1), "p_value": (f"{p:.2e}" if p > 0 else "<1e-300"),
    }

def main():
    print("Loading summary ...", flush=True)
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "First_Date",
                     "Classification_Changed", "Num_Submissions", "Molecular_Consequence",
                     "Review_Status", "Variant_Type", "Start", "Stop", "HGVS_p", "Num_PMIDs"],
                     dtype={"Gene": str, "First_Date": str})
    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    N = len(df); print(f"  {N:,} variants, {df['changed'].sum():,} changed", flush=True)

    cancer = set(DOMAINS["cancer_genes"]); common = set(DOMAINS["common_disease_genes"])
    gene2dom = {g: dom for dom, gs in DOMAINS["disease_domains"].items() for g in gs}
    df["is_cancer"] = df["Gene"].isin(cancer)
    df["is_common"] = df["Gene"].isin(common)
    df["domain"] = df["Gene"].map(gene2dom).fillna("Other/Unclassified")

    # ---------- Part A: univariate ----------
    uni = {}
    def grp(mask):  # (changed, total) for a boolean mask
        sub = df[mask]; return int(sub["changed"].sum()), int(len(sub))
    xc, nc = grp(df["is_cancer"]); xnc, nnc = grp(~df["is_cancer"])
    uni["cancer_vs_noncancer"] = two_by_two(xc, nc, xnc, nnc, "Cancer (COSMIC CGC)", "Non-cancer")
    xco, nco = grp(df["is_common"]); xr, nr = grp(~df["is_common"])
    uni["common_vs_rare"] = two_by_two(xco, nco, xr, nr, "Higher-prevalence Mendelian", "Rare disease")
    # each domain vs Other/Unclassified
    xo, no = grp(df["domain"] == "Other/Unclassified")
    for dom in DOMAINS["disease_domains"]:
        xd, nd = grp(df["domain"] == dom)
        uni[f"domain_{dom}_vs_other"] = two_by_two(xd, nd, xo, no, dom, "Other/Unclassified")
    json.dump(uni, open(OUT / "stats_univariate.json", "w"), indent=2)

    print("\n=== PART A: univariate effect sizes ===")
    for k, v in uni.items():
        print(f"{k}:")
        print(f"   {v['group1']} {v['rate1_pct']}% vs {v['group2']} {v['rate2_pct']}%  "
              f"RR={v['risk_ratio']} {v['rr_95ci']}  OR={v['odds_ratio']} {v['or_95ci']}  p={v['p_value']}")

    # ---------- HGVS_p completeness fix (R2 minor) ----------
    mc = df["Molecular_Consequence"].fillna("").str.lower()
    coding = mc.str.contains("missense|nonsense|synonymous|frameshift|inframe|stop|start|initiator|coding")
    has_p = df["HGVS_p"].notna() & (df["HGVS_p"].astype(str).str.strip() != "")
    hgvsp = {
        "hgvs_p_overall_pct": round(100 * has_p.mean(), 2),
        "hgvs_p_among_coding_pct": round(100 * has_p[coding].mean(), 2),
        "n_coding": int(coding.sum()), "n_total": int(N),
        "note": "coding = molecular consequence indicating a protein consequence"
    }
    json.dump(hgvsp, open(OUT / "hgvs_p_completeness.json", "w"), indent=2)
    print(f"\n=== HGVS_p completeness (R2 minor) ===\n   overall {hgvsp['hgvs_p_overall_pct']}%  "
          f"among coding variants {hgvsp['hgvs_p_among_coding_pct']}%  (n_coding={hgvsp['n_coding']:,})")

    # ---------- Part B: multivariable logistic regression ----------
    print("\n=== PART B: multivariable logistic regression ===", flush=True)
    import statsmodels.formula.api as smf
    d = df.copy()
    # confounders R2 named:
    d["log_subs"] = np.log10(d["Num_Submissions"].clip(lower=1))              # number of submitters
    yr = d["First_Date"].str.slice(0, 4).astype(float)
    d["follow_up_yr"] = (2025 + 7/12) - (yr + (d["First_Date"].str.slice(5,7).astype(float))/12)  # time in ClinVar
    d["era"] = pd.cut(yr, [0, 2014.99, 2019.99, 3000], labels=["pre2015", "RCV_2015_2019", "VCV_2020plus"])
    # gene size proxy = genomic span (bp), log10
    span = (d["Stop"].astype(float) - d["Start"].astype(float)).abs()
    gene_span = d.groupby("Gene").apply(lambda g: (g["Stop"].astype(float).max() - g["Start"].astype(float).min()))
    d["log_gene_span"] = np.log10(d["Gene"].map(gene_span).clip(lower=1) + 1)
    # variant class (collapsed)
    def vclass(s):
        s = str(s).lower()
        for key in ["missense", "nonsense", "frameshift", "splice", "synonymous", "inframe"]:
            if key in s: return key
        return "other"
    d["vclass"] = d["Molecular_Consequence"].map(vclass)
    # review tier (collapsed)
    def rtier(s):
        s = str(s).lower()
        if "conflicting" in s: return "conflicting"
        if "expert panel" in s or "practice guideline" in s: return "expert_guideline"
        if "multiple submitters" in s: return "multi_no_conflict"
        if "single submitter" in s: return "single"
        return "no_criteria"
    d["rtier"] = d["Review_Status"].map(rtier)

    cols = ["changed", "log_subs", "follow_up_yr", "era", "log_gene_span", "vclass", "rtier", "is_cancer"]
    d = d[cols].dropna()
    d["is_cancer"] = d["is_cancer"].astype(int)
    print(f"  modeling {len(d):,} variants", flush=True)
    formula = ("changed ~ log_subs + follow_up_yr + log_gene_span + is_cancer "
               "+ C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
               "+ C(rtier, Treatment('single'))")
    model = smf.logit(formula, data=d).fit(disp=1, maxiter=100)
    params = model.params; conf = model.conf_int(); pvals = model.pvalues
    rows = []
    for name in params.index:
        rows.append({"term": name, "OR": round(float(np.exp(params[name])), 3),
                     "OR_lo": round(float(np.exp(conf.loc[name, 0])), 3),
                     "OR_hi": round(float(np.exp(conf.loc[name, 1])), 3),
                     "p": (f"{pvals[name]:.2e}")})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "stats_multivariable.csv", index=False)
    print(out.to_string(index=False))
    print(f"\n  pseudo-R2 (McFadden): {model.prsquared:.4f}  |  n={int(model.nobs):,}")

if __name__ == "__main__":
    main()
