#!/usr/bin/env python
"""
Revision analysis #3 — Tiered / clinically-stratified reclassification rate.

Addresses (editor-essential) Reviewer 1 major + Reviewer 2 major #1:
  R1: concern the headline 9.3% conflates administrative events with true reclassification.
  R2: "counts any change between the nine categories ... including adjacent changes such as
       P>LP, LB>B, and VUS>conflicting. This may inflate the 9.3% ... Major discordance should
       be analyzed separately from less consequential disagreements."

We decompose the 342,057 reclassified variants (9.3%) by the MOST clinically significant
change each variant ever underwent, using ClinVar's own conflict-severity definitions.

Outputs: outputs/tiered_reclassification.json, figures/tiered_reclassification.png/.pdf
"""
import csv, json, collections
from pathlib import Path
import numpy as np

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
FIG  = DATA / "revision_analyses" / "figures"
TOTAL_VARIANTS = 3678303

# severity label for each (prev,curr) transition, from the project's own definitions
sevmap = {}
with open(PROJ / "conflict_severity_detail.csv") as fh:
    for row in csv.DictReader(fh):
        sevmap[(row["Prev_Category"], row["Curr_Category"])] = row["severity"]

# clinical ordering (most severe wins when a variant has several kinds of change)
rank = {"Medically significant": 4, "Involves composite category": 3,
        "Non-medically significant": 2, "Adjacent (P↔LP)": 1, "Adjacent (B↔LB)": 1}
inv = {4: "Medically significant (crosses actionability: P/LP↔VUS, P/LP↔B/LB)",
       3: "Involves composite/uncertain category (VUS↔Conflicting, ↔P/LP, ↔B/LB)",
       2: "Non-medically significant (VUS↔LB, VUS↔B)",
       1: "Within-tier adjacent only (P↔LP or B↔LB) — ClinVar-concordant"}

var_rank = collections.defaultdict(int)
var_has_medsig = set(); var_has_vusconf = set()
var_only_adjacent = collections.defaultdict(lambda: True)
with open(PROJ / "reclassification_events_detail.csv") as fh:
    for row in csv.DictReader(fh):
        vid = row["Variation_ID"]; key = (row["Prev_Category"], row["Curr_Category"])
        s = sevmap.get(key)
        if s is None:
            continue
        var_rank[vid] = max(var_rank[vid], rank.get(s, 0))
        if s == "Medically significant": var_has_medsig.add(vid)
        if key in (("VUS", "Conflicting"), ("Conflicting", "VUS")): var_has_vusconf.add(vid)
        if s not in ("Adjacent (P↔LP)", "Adjacent (B↔LB)"): var_only_adjacent[vid] = False

nvar = len(var_rank)
by = collections.Counter(var_rank.values())
only_adj = sum(1 for v in var_rank if var_only_adjacent[v])

res = {
    "total_variants": TOTAL_VARIANTS,
    "reclassified_any_change": {"n": nvar, "pct_of_all": round(100*nvar/TOTAL_VARIANTS, 2)},
    "partition_by_most_severe_change": {
        inv[r]: {"n": by[r], "pct_of_all": round(100*by[r]/TOTAL_VARIANTS, 2),
                 "pct_of_reclassified": round(100*by[r]/nvar, 1)}
        for r in sorted(by, reverse=True)},
    "variants_with_any_medically_significant_change": {
        "n": len(var_has_medsig), "pct_of_all": round(100*len(var_has_medsig)/TOTAL_VARIANTS, 2)},
    "variants_within_tier_adjacent_only": {
        "n": only_adj, "pct_of_reclassified": round(100*only_adj/nvar, 1),
        "note": "arguably should NOT count as reclassification (ClinVar treats P/LP and B/LB as concordant)"},
    "variants_involving_VUS_conflicting": {
        "n": len(var_has_vusconf), "pct_of_reclassified": round(100*len(var_has_vusconf)/nvar, 1)},
}
json.dump(res, open(OUT / "tiered_reclassification.json", "w"), indent=2, ensure_ascii=False)
print(json.dumps(res, indent=2, ensure_ascii=False))

# ---- figure: waterfall from 9.3% down to the clinically meaningful core ----
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
labels = ["Any category\nchange\n(9.30%)",
          "Within-tier\nadjacent only\n(excluded)",
          "Non-medically\nsignificant",
          "Composite/\nuncertain",
          "Medically\nsignificant\n(0.68%)"]
med = by[4]; comp = by[3]; nonmed = by[2]; adj = by[1]
pct = lambda n: 100*n/TOTAL_VARIANTS
vals = [pct(nvar), pct(adj), pct(nonmed), pct(comp), pct(med)]
colors = ["#9aa5b1", "#d9d9d9", "#c6dbef", "#fdae6b", "#b2182b"]
fig, ax = plt.subplots(figsize=(8.6, 4.6))
bars = ax.bar(labels, vals, color=colors, edgecolor="white")
for b, n in zip(bars, [nvar, adj, nonmed, comp, med]):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.05,
            f"{n:,}\n({100*n/TOTAL_VARIANTS:.2f}%)", ha="center", va="bottom", fontsize=8)
ax.set_ylabel("% of all 3.68M ClinVar variants")
ax.set_ylim(0, 10.2)
ax.set_title("Decomposition of the 9.3% reclassification rate by clinical significance")
ax.grid(axis="y", alpha=0.25)
plt.tight_layout()
fig.savefig(FIG / "tiered_reclassification.png", dpi=300)
fig.savefig(FIG / "tiered_reclassification.pdf")
print("\nFigure written to figures/tiered_reclassification.png")
