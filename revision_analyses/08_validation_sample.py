#!/usr/bin/env python
"""
Revision analysis #8 — build a stratified random sample for MANUAL validation of the
reclassification-mechanism heuristic.

Addresses Reviewer 2, Major point #3:
  "the inferred mechanisms are not directly observed. The authors should validate the heuristic
   against a manually reviewed sample of SCV/VCV histories and report misclassification rates or
   uncertainty intervals for the estimated mechanism proportions."

Produces a review spreadsheet: a reproducible (fixed-seed) stratified random sample of events
across every mechanism label, each with a direct ClinVar link to the variant's full submission
history, the heuristic's label + signals, and blank columns for the human adjudicator. Once the
co-authors fill in `true_mechanism` and `agree`, a companion scorer computes per-class
precision/recall and a Wilson CI on overall accuracy.

Output: outputs/validation_sample_for_manual_review.csv
"""
import csv, random, collections
from pathlib import Path

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
EVENTS = PROJ / "reclassification_events_detail.csv"
PER_LABEL = 6            # ~50 total across 8 mechanism labels (manual-review burden)
SEED = 20260710

random.seed(SEED)
by_label = collections.defaultdict(list)
with open(EVENTS) as fh:
    for row in csv.DictReader(fh):
        by_label[row["Driver_Label"]].append(row)

sample = []
for label, rows in sorted(by_label.items()):
    k = min(PER_LABEL, len(rows))
    sample.extend(random.sample(rows, k))
    print(f"  {label:24s} population={len(rows):>7,}  sampled={k}")

random.shuffle(sample)   # blind the reviewer to label ordering
cols = ["review_id", "clinvar_link", "Variation_ID", "Gene", "Prev_Date", "Curr_Date",
        "Prev_Category", "Curr_Category", "Prev_Num_Subs", "Curr_Num_Subs",
        "heuristic_label", "heuristic_signals",
        "true_mechanism_(fill_in)", "agree_Y_N_(fill_in)", "notes_(fill_in)"]
with open(OUT / "validation_sample_for_manual_review.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for i, r in enumerate(sample, 1):
        vid = r["Variation_ID"]
        w.writerow([f"R{i:04d}",
                    f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/",
                    vid, r["Gene"], r["Prev_Date"], r["Curr_Date"],
                    r["Prev_Category"], r["Curr_Category"], r["Prev_Num_Subs"], r["Curr_Num_Subs"],
                    r["Driver_Label"], r["Signals"], "", "", ""])
print(f"\nWrote {len(sample)} events to outputs/validation_sample_for_manual_review.csv")
print("Reviewers: open each clinvar_link, inspect the submission history, and fill the last 3 columns.")
print("Then run 09_validation_score.py to get per-class precision/recall + accuracy CI.")
