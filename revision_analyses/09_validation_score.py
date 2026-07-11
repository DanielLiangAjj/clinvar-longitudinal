#!/usr/bin/env python
"""
Revision analysis #9 — score the manually-reviewed validation sample (run AFTER co-authors
fill in `true_mechanism_(fill_in)` and `agree_Y_N_(fill_in)` in the sample CSV).

Computes: overall accuracy with Wilson 95% CI, per-class precision/recall/F1, and a confusion
matrix — the "misclassification rates / uncertainty intervals" Reviewer 2 asked for. These
numbers let us attach uncertainty to the headline mechanism proportions (e.g., reinterpretation
33.2% +/- validated error).
"""
import csv, math, collections
from pathlib import Path

OUT = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/outputs")
SAMPLE = OUT / "validation_sample_for_manual_review.csv"

def wilson(k, n, z=1.96):
    if n == 0: return (0, 0, 0)
    p = k / n; denom = 1 + z*z/n
    center = (p + z*z/(2*n)) / denom
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return p, center-half, center+half

def main():
    rows = list(csv.DictReader(open(SAMPLE)))
    filled = [r for r in rows if r.get("true_mechanism_(fill_in)", "").strip()]
    if not filled:
        print(f"No rows adjudicated yet in {SAMPLE.name}.")
        print("Fill 'true_mechanism_(fill_in)' (and optionally 'agree_Y_N_(fill_in)') then re-run.")
        return
    n = len(filled)
    correct = sum(1 for r in filled if r["heuristic_label"].strip() == r["true_mechanism_(fill_in)"].strip())
    p, lo, hi = wilson(correct, n)
    print(f"Adjudicated: {n} / {len(rows)}")
    print(f"Overall accuracy: {correct}/{n} = {100*p:.1f}%  (Wilson 95% CI {100*lo:.1f}-{100*hi:.1f})")

    labels = sorted(set([r["heuristic_label"].strip() for r in filled] +
                        [r["true_mechanism_(fill_in)"].strip() for r in filled]))
    cm = collections.defaultdict(lambda: collections.Counter())
    for r in filled:
        cm[r["true_mechanism_(fill_in)"].strip()][r["heuristic_label"].strip()] += 1
    print("\nPer-class precision / recall / F1:")
    for lab in labels:
        tp = cm[lab][lab]
        fp = sum(cm[t][lab] for t in labels if t != lab)
        fn = sum(cm[lab][pred] for pred in labels if pred != lab)
        prec = tp/(tp+fp) if tp+fp else float("nan")
        rec = tp/(tp+fn) if tp+fn else float("nan")
        f1 = 2*prec*rec/(prec+rec) if prec and rec and (prec+rec) else float("nan")
        print(f"  {lab:24s} P={prec:5.2f}  R={rec:5.2f}  F1={f1:5.2f}  (n_true={sum(cm[lab].values())})")

if __name__ == "__main__":
    main()
