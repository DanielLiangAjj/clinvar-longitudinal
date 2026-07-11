#!/usr/bin/env python
"""
Revision analysis #7 — PubMed citation distribution by domain (fixes Fig 9g) + skew summary.

Addresses Reviewer 2 minor (Fig 9g):
  "should display variability around mean PubMed citations per variant, such as standard
   deviation or interquartile range, or preferably distribution plots because citation counts
   are likely highly skewed."

We replace the mean-only bar with per-domain distribution summaries (median, IQR, mean, SD,
% with >=1 citation) and a box/strip figure on a log scale.

NOTE on R2 major #6 (filtering guideline vs variant-specific citations): the raw per-variant
PMID *lists* are not in the local exports (only Num_PMIDs counts). Filtering out guideline
PMIDs (e.g., ACMG/AMP 2015 = 25741868, SHERLOC = 28492532) requires a variant->PMID export
from the server checkpoints. This script quantifies the skew now; guideline filtering is
flagged as needing that one additional file.

Outputs: outputs/citation_distribution_by_domain.json, figures/fig9g_distribution.png/.pdf
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
FIG  = DATA / "revision_analyses" / "figures"
SUMMARY = DATA / "summary_20260219_024805.csv"
DOMAINS = json.load(open(PROJ / "disease_domains.json"))

def main():
    df = pd.read_csv(SUMMARY, usecols=["Gene", "Num_PMIDs"], dtype={"Gene": str})
    df["Num_PMIDs"] = pd.to_numeric(df["Num_PMIDs"], errors="coerce").fillna(0)
    cancer = set(DOMAINS["cancer_genes"])
    gene2dom = {g: d for d, gs in DOMAINS["disease_domains"].items() for g in gs}
    df["domain"] = df["Gene"].map(gene2dom).fillna("Other/Unclassified")
    df["cancer"] = np.where(df["Gene"].isin(cancer), "Cancer (COSMIC)", "Non-cancer")

    def stats(x):
        x = x.values
        return {"n": int(len(x)), "mean": round(float(x.mean()), 3), "sd": round(float(x.std()), 3),
                "median": float(np.median(x)), "q25": float(np.percentile(x, 25)),
                "q75": float(np.percentile(x, 75)), "q95": float(np.percentile(x, 95)),
                "max": float(x.max()), "pct_ge1": round(100*float((x >= 1).mean()), 1)}

    res = {"by_domain": {}, "by_cancer": {}}
    order = ["Hereditary Cancer", "Cardiovascular", "Neurological", "Renal", "Ophthalmologic", "Other/Unclassified"]
    for d in order:
        sub = df[df["domain"] == d]["Num_PMIDs"]
        if len(sub): res["by_domain"][d] = stats(sub)
    for c in ["Cancer (COSMIC)", "Non-cancer"]:
        res["by_cancer"][c] = stats(df[df["cancer"] == c]["Num_PMIDs"])
    json.dump(res, open(OUT / "citation_distribution_by_domain.json", "w"), indent=2)

    print("=== PubMed citations per variant — distribution (highly skewed) ===")
    print(f"{'domain':22s}{'mean':>7s}{'SD':>7s}{'median':>8s}{'IQR':>12s}{'%>=1':>7s}")
    for d, s in res["by_domain"].items():
        print(f"{d:22s}{s['mean']:>7.2f}{s['sd']:>7.2f}{s['median']:>8.0f}   [{s['q25']:.0f}-{s['q75']:.0f}]      {s['pct_ge1']:>5.1f}")

    # figure: box on symlog, by domain
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    data = [df[df["domain"] == d]["Num_PMIDs"].values for d in order if (df["domain"] == d).any()]
    labels = [d.replace(" ", "\n") for d in order if (df["domain"] == d).any()]
    bp = ax.boxplot(data, labels=labels, showfliers=False, patch_artist=True, whis=(5, 95))
    for patch in bp["boxes"]:
        patch.set_facecolor("#c6dbef"); patch.set_edgecolor("#2166ac")
    # overlay mean as red diamond
    means = [d.mean() for d in data]
    ax.plot(range(1, len(means)+1), means, "D", color="#b2182b", ms=6, label="mean", zorder=5)
    ax.set_yscale("symlog"); ax.set_ylabel("PubMed citations per variant (symlog)")
    ax.set_title("Fig 9g (revised): citation distribution by domain\n(box = median/IQR/5-95%, red = mean)")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.2)
    plt.tight_layout(); fig.savefig(FIG/"fig9g_distribution.png", dpi=300); fig.savefig(FIG/"fig9g_distribution.pdf")
    print("\nfigure written: figures/fig9g_distribution.png")

if __name__ == "__main__":
    main()
