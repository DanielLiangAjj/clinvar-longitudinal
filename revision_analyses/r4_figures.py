#!/usr/bin/env python
"""R4 (editorial round): regenerate the figures flagged by the editor (items 7-9)
with the font/overlap fixes patched into the renderers. Outputs to figures/r4/."""
import json, sys
from pathlib import Path
REV = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses")
EXT = REV / "external_dependencies" / "ClinVar Research Project"
OUT = REV / "outputs"
R4 = REV / "figures" / "r4"; R4.mkdir(parents=True, exist_ok=True)
CACHE = Path("/Users/danielliang/Downloads/ClinVar_Research_Project/figures/.cache/summary_agg.json")
DATA2 = Path("/Users/danielliang/Downloads/ClinVar_Research_Project/data")
sys.path.insert(0, str(EXT / "scripts"))
import clinvar_visualizations as V

agg = json.load(open(CACHE))
cit = json.load(open(OUT / "citation_r3_numbers.json"))
# filtered keys for figures 11/12 (module figure4/figure7) — recompute sums from saved numbers? they
# were computed in 28_; reuse the r3 approach: filtered bins/sums were not persisted, so rebuild from
# the per-variant join is expensive. Instead reuse the r3-rendered inputs: re-run the cheap parts.
import csv
from collections import Counter
EXCL = {"28492532", "25741868"}
filt = {}
with open(OUT / "variant_pmids_latest.csv") as fh:
    for row in csv.DictReader(fh):
        s = set(p for p in row["PMIDs"].split(";") if p)
        filt[row["Variation_ID"]] = len(s - EXCL)
import pandas as pd
df = pd.read_csv("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/summary_20260219_024805.csv",
                 usecols=["Variation_ID", "First_Date"], dtype=str)
df["f"] = df["Variation_ID"].map(filt).fillna(0).astype(int)
def bin6(n):
    if n == 0: return "0"
    if n == 1: return "1"
    if n <= 5: return "2-5"
    if n <= 10: return "6-10"
    if n <= 20: return "11-20"
    return ">20"
agg["pmid_bins_filtered"] = df["f"].map(bin6).value_counts().to_dict()
agg["pmid_sum_by_first_date_filtered"] = {k: int(v) for k, v in df.groupby("First_Date")["f"].sum().items()}
print("filtered aggregates rebuilt")

V.figure4(agg, {}, R4)          # manuscript Figure 11
V.figure7(agg, R4)              # manuscript Figure 12

domain = json.load(open(EXT / "data" / "domain_summary.json"))
for d, s in cit["by_domain"].items():
    if d in domain.get("by_domain", {}):
        domain["by_domain"][d]["cit_classes"] = s["cit_classes"]
for sg in ["cancer", "noncancer"]:
    domain["cancer_vs_noncancer"][sg]["mean_pmids_filtered"] = cit["subgroups"][sg]["mean_filtered_all"]
for sg in ["rare", "common"]:
    domain["rare_vs_common"][sg]["mean_pmids_filtered"] = cit["subgroups"][sg]["mean_filtered_all"]
V.figure10_disease_domain(domain, R4)   # manuscript Figure 9

stats = json.load(open(DATA2 / "statistics_20260219_024805.json"))
V.figure6(agg, stats, R4)       # manuscript Figure 13
acmg = json.load(open(DATA2 / "acmg_criteria_detailed.json"))
V.figure9(acmg, R4)             # manuscript Figure 14
drivers = json.load(open(DATA2 / "reclassification_drivers.json"))
V.figure11_reclassification_drivers(drivers, R4)   # manuscript Figure 6
print("done; files in", R4)
for p in sorted(R4.glob("*.png")): print("  ", p.name, p.stat().st_size // 1024, "KB")

# Additional file 1 Fig. S3 (batch-event decomposition), labels folded into tick labels
V.figure13_batch_event_decomposition(str(EXT / "data"), R4)

# Figures 8 and 15: italic gene symbols; alias rows in Fig 8b/c marked "(alt. name)"
DATA_LOCAL = "/Users/danielliang/Downloads/ClinVar_Research_Project/data"
V.figure14_submitter_landscape(DATA_LOCAL, R4)      # writes figure14 (Fig 7) + figure15 deep dive (Fig 8)
V.figure16_case_study_trajectories(DATA_LOCAL, R4)  # Fig 15
