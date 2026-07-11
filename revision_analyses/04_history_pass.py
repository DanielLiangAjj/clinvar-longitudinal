#!/usr/bin/env python
"""
Revision analysis #4 — single streaming pass over the 14 GB classification-history file.

Produces two things needed downstream:
  (A) baseline_covariates.csv  — per variant, the values AT FIRST APPEARANCE
      (baseline Num_Submissions, review tier, category, date). Needed so the multivariable
      model (R2 #4) adjusts for BASELINE state, not post-outcome (final-snapshot) values.
  (B) monthly_timeseries.csv   — per snapshot month, aggregate counts for the interrupted
      time-series (R2 #5): total variants, multi-submitter variants, conflicting, concordant,
      and category composition (VUS / P-LP / B-LB / Conflicting).

The file is grouped by Variation_ID with dates ascending within each group; we stream it
variant-by-variant (O(1) memory) and dedup repeated (variant, month) rows.
"""
import csv, collections
from pathlib import Path

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
OUT  = DATA / "revision_analyses" / "outputs"
HIST = DATA / "classification_history_20260219_024805.csv"

def review_tier(s):
    s = s.lower()
    if "conflicting" in s: return "conflicting"
    if "expert panel" in s or "practice guideline" in s: return "expert_guideline"
    if "multiple submitter" in s: return "multi_no_conflict"
    if "single submitter" in s: return "single"
    return "no_criteria"

# monthly counters
m_total = collections.Counter()
m_multi = collections.Counter()
m_conflicting = collections.Counter()   # review status conflicting
m_concordant = collections.Counter()    # multi-submitter, no conflict / expert / guideline
m_cat = collections.defaultdict(collections.Counter)   # month -> Counter(category)

def flush_baseline(bw, vid, gene, row):
    """row = (date, category, review_status, num_subs) at first appearance."""
    date, cat, rs, ns = row
    bw.writerow([vid, gene, date, cat, review_tier(rs), ns])

def main():
    n_rows = 0; n_variants = 0
    cur_vid = None; cur_gene = None; cur_min_date = None; cur_base = None
    seen_month_for_var = set()   # (dedup within current variant) months already counted

    bfh = open(OUT / "baseline_covariates.csv", "w", newline="")
    bw = csv.writer(bfh)
    bw.writerow(["Variation_ID", "Gene", "base_date", "base_category", "base_review_tier", "base_num_subs"])

    with open(HIST, newline="") as fh:
        r = csv.reader(fh)
        header = next(r)
        idx = {name: i for i, name in enumerate(header)}
        iV, iG, iD, iCat, iRS, iNS = (idx["Variation_ID"], idx["Gene"], idx["Snapshot_Date"],
                                      idx["Category"], idx["Review_Status"], idx["Num_Submissions"])
        for row in r:
            n_rows += 1
            vid = row[iV]
            if vid != cur_vid:
                # close previous variant: write its baseline
                if cur_vid is not None:
                    flush_baseline(bw, cur_vid, cur_gene, cur_base)
                    n_variants += 1
                cur_vid = vid; cur_gene = row[iG]
                cur_min_date = row[iD]; cur_base = (row[iD], row[iCat], row[iRS], row[iNS])
                seen_month_for_var = set()
            else:
                # track earliest-date row as baseline (robust to any local disorder)
                if row[iD] < cur_min_date:
                    cur_min_date = row[iD]; cur_base = (row[iD], row[iCat], row[iRS], row[iNS])

            date = row[iD]
            if date in seen_month_for_var:
                continue                      # dedup repeated (variant, month)
            seen_month_for_var.add(date)

            m_total[date] += 1
            try:
                ns = int(row[iNS])
            except ValueError:
                ns = 1
            rs = row[iRS].lower()
            cat = row[iCat]
            m_cat[date][cat] += 1
            if ns > 1 or "multiple submitter" in rs:
                m_multi[date] += 1
                if "no conflict" in rs or "expert panel" in rs or "practice guideline" in rs:
                    m_concordant[date] += 1
            if "conflicting" in rs:
                m_conflicting[date] += 1

            if n_rows % 20_000_000 == 0:
                print(f"  ...{n_rows:,} rows, {n_variants:,} variants", flush=True)

        # last variant
        if cur_vid is not None:
            flush_baseline(bw, cur_vid, cur_gene, cur_base)
            n_variants += 1
    bfh.close()

    # write monthly time series
    months = sorted(m_total)
    with open(OUT / "monthly_timeseries.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["month", "n_variants", "n_multi", "n_conflicting_reviewstatus", "n_concordant_multi",
                    "n_VUS", "n_PLP", "n_BLB", "n_Conflicting_cat",
                    "pct_conflicting_of_all", "pct_concordant_of_multi", "pct_VUS"])
        for mo in months:
            tot = m_total[mo]; multi = m_multi[mo]
            cats = m_cat[mo]
            vus = cats.get("VUS", 0)
            plp = cats.get("P", 0) + cats.get("LP", 0) + cats.get("P/LP", 0)
            blb = cats.get("B", 0) + cats.get("LB", 0) + cats.get("B/LB", 0)
            confcat = cats.get("Conflicting", 0)
            w.writerow([mo, tot, multi, m_conflicting[mo], m_concordant[mo],
                        vus, plp, blb, confcat,
                        round(100*m_conflicting[mo]/tot, 3) if tot else 0,
                        round(100*m_concordant[mo]/multi, 3) if multi else 0,
                        round(100*vus/tot, 3) if tot else 0])
    print(f"\nDONE. rows={n_rows:,}  variants={n_variants:,}  months={len(months)}")
    print(f"  baseline_covariates.csv and monthly_timeseries.csv written to {OUT}")

if __name__ == "__main__":
    main()
