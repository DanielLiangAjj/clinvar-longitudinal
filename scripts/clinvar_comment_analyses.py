#!/usr/bin/env python3
"""
ClinVar Comment-Driven Analyses
================================
Addresses reviewer comments #16 and #17 from Yiming Luo,
plus extracts P/LP↔B/LB reversal variant details.

Analysis 1 (Comment #16): Time-to-reclassification statistics
Analysis 2 (Comment #17): Batch event decomposition
Analysis 3: P/LP↔B/LB reversal variant extraction with full history

Usage:
    python3 clinvar_comment_analyses.py \
        --summary-csv  /path/to/summary_*.csv \
        --events-csv   /path/to/reclassification_events_detail.csv \
        --output-dir   /path/to/output
"""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

# ── Constants ────────────────────────────────────────────────────────────

P_SIDE = {"P", "LP", "P/LP"}
B_SIDE = {"B", "LB", "B/LB"}

BATCH_DATES = {"2013-12", "2014-01", "2014-02", "2018-06", "2020-01"}


def _months_between(d1: str, d2: str) -> float:
    """Return months between two YYYY-MM date strings."""
    y1, m1 = int(d1[:4]), int(d1[5:7])
    y2, m2 = int(d2[:4]), int(d2[5:7])
    return (y2 - y1) * 12 + (m2 - m1)


def _years_between(d1: str, d2: str) -> float:
    return _months_between(d1, d2) / 12.0


# ═════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Time-to-Reclassification (Comment #16)
# ═════════════════════════════════════════════════════════════════════════

def analysis_time_to_reclass(summary_path: Path, events_path: Path,
                              out_dir: Path):
    """Compute time from first ClinVar appearance to first reclassification."""
    print("\n" + "=" * 70)
    print("ANALYSIS 1: Time-to-Reclassification (Comment #16)")
    print("=" * 70)

    # Step 1: Load first appearance date from summary CSV
    print("  Loading variant first-appearance dates from summary CSV ...")
    first_date = {}  # variant_id -> first_date
    first_cat = {}   # variant_id -> first_category
    with open(summary_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            vid = row["Variation_ID"]
            changed = row.get("Classification_Changed", "False")
            if changed == "True":
                first_date[vid] = row["First_Date"]
                first_cat[vid] = row.get("First_Classification", "")

    print(f"    Reclassified variants with first dates: {len(first_date):,}")

    # Step 2: Find earliest reclassification event per variant
    print("  Loading first reclassification event from events CSV ...")
    first_event = {}  # variant_id -> (curr_date, prev_cat, curr_cat)
    with open(events_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            vid = row["Variation_ID"]
            curr_date = row["Curr_Date"]
            if vid not in first_event or curr_date < first_event[vid][0]:
                first_event[vid] = (
                    curr_date,
                    row["Prev_Category"],
                    row["Curr_Category"],
                )

    print(f"    Variants with reclassification events: {len(first_event):,}")

    # Step 3: Compute time-to-first-reclassification
    times_all = []
    times_by_transition = defaultdict(list)

    for vid, (event_date, prev, curr) in first_event.items():
        if vid not in first_date:
            continue
        fd = first_date[vid]
        if fd >= event_date:
            continue  # skip if dates don't make sense

        years = _years_between(fd, event_date)
        if years < 0:
            continue

        times_all.append(years)
        transition_key = f"{prev}→{curr}"
        times_by_transition[transition_key].append(years)

    times_all = np.array(times_all)
    print(f"\n    Variants with valid time-to-reclassification: {len(times_all):,}")

    # Step 4: Compute statistics
    results = {
        "overall": {
            "n": len(times_all),
            "mean_years": float(np.mean(times_all)),
            "median_years": float(np.median(times_all)),
            "q25_years": float(np.percentile(times_all, 25)),
            "q75_years": float(np.percentile(times_all, 75)),
            "min_years": float(np.min(times_all)),
            "max_years": float(np.max(times_all)),
            "pct_within_1yr": float(np.mean(times_all <= 1.0) * 100),
            "pct_within_2yr": float(np.mean(times_all <= 2.0) * 100),
            "pct_within_5yr": float(np.mean(times_all <= 5.0) * 100),
        },
        "by_transition": {},
    }

    print(f"\n  ── Overall Statistics ──")
    print(f"    Mean:   {results['overall']['mean_years']:.2f} years")
    print(f"    Median: {results['overall']['median_years']:.2f} years")
    print(f"    IQR:    {results['overall']['q25_years']:.2f} – "
          f"{results['overall']['q75_years']:.2f} years")
    print(f"    Within 1 year:  {results['overall']['pct_within_1yr']:.1f}%")
    print(f"    Within 2 years: {results['overall']['pct_within_2yr']:.1f}%")
    print(f"    Within 5 years: {results['overall']['pct_within_5yr']:.1f}%")

    # Top transitions by frequency
    top_transitions = sorted(times_by_transition.items(),
                             key=lambda x: -len(x[1]))[:20]

    print(f"\n  ── By Transition Type (top 20) ──")
    print(f"    {'Transition':<25s} {'N':>8s} {'Mean':>7s} {'Median':>7s} "
          f"{'IQR':>13s}")
    print(f"    {'─' * 25} {'─' * 8} {'─' * 7} {'─' * 7} {'─' * 13}")

    for tkey, tvals in top_transitions:
        arr = np.array(tvals)
        stats = {
            "n": len(arr),
            "mean_years": float(np.mean(arr)),
            "median_years": float(np.median(arr)),
            "q25_years": float(np.percentile(arr, 25)),
            "q75_years": float(np.percentile(arr, 75)),
        }
        results["by_transition"][tkey] = stats
        print(f"    {tkey:<25s} {stats['n']:>8,d} {stats['mean_years']:>6.2f}y "
              f"{stats['median_years']:>6.2f}y "
              f"{stats['q25_years']:.1f}–{stats['q75_years']:.1f}y")

    # Clinically important transitions
    print(f"\n  ── Clinically Important Transitions ──")
    important = {
        "VUS→P/LP": ["VUS→P", "VUS→LP", "VUS→P/LP"],
        "VUS→B/LB": ["VUS→B", "VUS→LB", "VUS→B/LB"],
        "P/LP→B/LB (drastic)": [f"{p}→{b}" for p in P_SIDE for b in B_SIDE],
        "B/LB→P/LP (drastic)": [f"{b}→{p}" for b in B_SIDE for p in P_SIDE],
        "VUS→Conflicting": ["VUS→Conflicting"],
    }
    results["clinically_important"] = {}
    for label, keys in important.items():
        combined = []
        for k in keys:
            if k in times_by_transition:
                combined.extend(times_by_transition[k])
        if combined:
            arr = np.array(combined)
            stats = {
                "n": len(arr),
                "mean_years": float(np.mean(arr)),
                "median_years": float(np.median(arr)),
                "q25_years": float(np.percentile(arr, 25)),
                "q75_years": float(np.percentile(arr, 75)),
            }
            results["clinically_important"][label] = stats
            print(f"    {label:<30s}  n={stats['n']:>7,d}  "
                  f"median={stats['median_years']:.2f}y  "
                  f"IQR={stats['q25_years']:.1f}–{stats['q75_years']:.1f}y")

    # Save results
    out_file = out_dir / "time_to_reclassification.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_file}")

    # Generate distribution data for histogram
    # Bin by 0.5-year intervals
    bins = np.arange(0, max(15.5, np.max(times_all) + 0.5), 0.5)
    hist, bin_edges = np.histogram(times_all, bins=bins)
    hist_data = {
        "bin_edges_years": bin_edges.tolist(),
        "counts": hist.tolist(),
    }
    hist_file = out_dir / "time_to_reclassification_histogram.json"
    with open(hist_file, "w") as f:
        json.dump(hist_data, f, indent=2)
    print(f"  Saved: {hist_file}")

    return results


# ═════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Batch Event Decomposition (Comment #17)
# ═════════════════════════════════════════════════════════════════════════

def analysis_batch_events(events_path: Path, out_dir: Path):
    """Decompose each batch event to determine ClinVar-wide vs lab-specific."""
    print("\n" + "=" * 70)
    print("ANALYSIS 2: Batch Event Decomposition (Comment #17)")
    print("=" * 70)

    # Load all batch events
    batch_events = defaultdict(list)  # date -> list of event dicts
    with open(events_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["Driver_Label"] == "batch_event":
                batch_events[row["Curr_Date"]].append(row)

    print(f"  Batch event dates found: {sorted(batch_events.keys())}")
    print(f"  Total batch event records: "
          f"{sum(len(v) for v in batch_events.values()):,}")

    results = {}
    for date in sorted(batch_events.keys()):
        events = batch_events[date]
        n = len(events)

        # Submission count analysis
        subs_stable = sum(
            1 for e in events
            if e["Prev_Num_Subs"] == e["Curr_Num_Subs"]
        )
        subs_increased = sum(
            1 for e in events
            if int(e["Curr_Num_Subs"]) > int(e["Prev_Num_Subs"])
        )
        subs_decreased = sum(
            1 for e in events
            if int(e["Curr_Num_Subs"]) < int(e["Prev_Num_Subs"])
        )

        # Gene diversity
        genes = Counter(e["Gene"] for e in events if e["Gene"])
        n_unique_genes = len(genes)
        top_10_genes = genes.most_common(10)

        # Classification change patterns
        change_patterns = Counter(
            f"{e['Prev_Category']}→{e['Curr_Category']}" for e in events
        )
        top_changes = change_patterns.most_common(10)

        # Signal analysis
        all_signals = Counter()
        for e in events:
            for sig in e.get("Signals", "").split("|"):
                if sig:
                    all_signals[sig] += 1

        # Determine likely source
        if subs_stable / max(n, 1) > 0.8:
            likely_source = "ClinVar-wide re-processing (submissions unchanged)"
        elif subs_increased / max(n, 1) > 0.5:
            likely_source = "Large-scale lab submission (submissions increased)"
        elif subs_decreased / max(n, 1) > 0.3:
            likely_source = "ClinVar re-aggregation (submissions decreased)"
        else:
            likely_source = "Mixed (combination of factors)"

        # Gene concentration (Herfindahl index)
        if n > 0:
            gene_shares = np.array([c / n for c in genes.values()])
            hhi = float(np.sum(gene_shares ** 2))
        else:
            hhi = 0.0

        batch_result = {
            "date": date,
            "total_events": n,
            "subs_stable": subs_stable,
            "subs_stable_pct": round(subs_stable / max(n, 1) * 100, 1),
            "subs_increased": subs_increased,
            "subs_increased_pct": round(subs_increased / max(n, 1) * 100, 1),
            "subs_decreased": subs_decreased,
            "subs_decreased_pct": round(subs_decreased / max(n, 1) * 100, 1),
            "unique_genes": n_unique_genes,
            "gene_concentration_hhi": round(hhi, 4),
            "top_10_genes": [{"gene": g, "count": c} for g, c in top_10_genes],
            "top_10_changes": [{"pattern": p, "count": c}
                               for p, c in top_changes],
            "signals": dict(all_signals.most_common(10)),
            "likely_source": likely_source,
        }
        results[date] = batch_result

        print(f"\n  ── {date} ({n:,} events) ──")
        print(f"    Subs stable:    {subs_stable:>6,d} ({batch_result['subs_stable_pct']}%)")
        print(f"    Subs increased: {subs_increased:>6,d} ({batch_result['subs_increased_pct']}%)")
        print(f"    Subs decreased: {subs_decreased:>6,d} ({batch_result['subs_decreased_pct']}%)")
        print(f"    Unique genes:   {n_unique_genes:,d}  "
              f"(HHI: {hhi:.4f}  {'concentrated' if hhi > 0.05 else 'dispersed'})")
        print(f"    Top genes:      {', '.join(f'{g}({c})' for g, c in top_10_genes[:5])}")
        print(f"    Top changes:    {', '.join(f'{p}({c})' for p, c in top_changes[:5])}")
        print(f"    → Likely: {likely_source}")

    out_file = out_dir / "batch_event_decomposition.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_file}")

    return results


# ═════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: P/LP ↔ B/LB Reversal Variant Extraction
# ═════════════════════════════════════════════════════════════════════════

def analysis_drastic_reversals(summary_path: Path, events_path: Path,
                                out_dir: Path):
    """Extract all variants with P/LP↔B/LB reversals and their full history."""
    print("\n" + "=" * 70)
    print("ANALYSIS 3: P/LP ↔ B/LB Reversal Extraction")
    print("=" * 70)

    # Step 1: Identify reversal variants from summary CSV
    print("  Scanning summary CSV for reversal variants ...")
    reversal_variants = {}  # vid -> summary row dict
    p_to_b_ids = set()
    b_to_p_ids = set()

    with open(summary_path, "r") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            if row.get("Classification_Changed") != "True":
                continue

            change_type = row.get("Change_Type", "")
            first_cat = _normalize_to_side(row.get("First_Classification", ""))
            curr_cat = _normalize_to_side(row.get("Current_Category", ""))

            # Check if the transition path includes a P-side to B-side
            # or B-side to P-side transition
            path_cats = change_type.split("->") if change_type else []
            has_p_to_b = False
            has_b_to_p = False

            for i in range(len(path_cats) - 1):
                prev_side = _cat_to_side(path_cats[i].strip())
                next_side = _cat_to_side(path_cats[i + 1].strip())
                if prev_side == "P" and next_side == "B":
                    has_p_to_b = True
                elif prev_side == "B" and next_side == "P":
                    has_b_to_p = True

            if has_p_to_b or has_b_to_p:
                vid = row["Variation_ID"]
                reversal_variants[vid] = row
                if has_p_to_b:
                    p_to_b_ids.add(vid)
                if has_b_to_p:
                    b_to_p_ids.add(vid)

    print(f"    Total reversal variants: {len(reversal_variants):,}")
    print(f"    P/LP → B/LB:            {len(p_to_b_ids):,}")
    print(f"    B/LB → P/LP:            {len(b_to_p_ids):,}")
    print(f"    Both directions:         "
          f"{len(p_to_b_ids & b_to_p_ids):,}")

    # Step 2: Get full reclassification event history for reversal variants
    print("  Loading reclassification event history for reversal variants ...")
    reversal_events = defaultdict(list)  # vid -> list of event rows
    with open(events_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            vid = row["Variation_ID"]
            if vid in reversal_variants:
                reversal_events[vid].append(row)

    # Sort events chronologically per variant
    for vid in reversal_events:
        reversal_events[vid].sort(key=lambda x: x["Curr_Date"])

    # Step 3: Write comprehensive output CSV
    out_csv = out_dir / "plp_blb_reversal_variants.csv"
    out_fields = [
        "Variation_ID", "VCV_ID", "Gene", "Variant_Type",
        "Molecular_Consequence", "HGVS_c", "HGVS_p", "dbSNP",
        "Chromosome", "Assembly",
        "First_Classification", "First_Date",
        "Current_Classification", "Current_Category", "Review_Status",
        "Num_Submissions", "Num_PMIDs", "Num_ACMG_Codes",
        "Change_Type", "Num_Changes", "Had_Oscillation",
        "Reversal_Direction",  # P→B, B→P, or Both
        "Num_Reclass_Events", "Event_History",
    ]

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()

        for vid, summary_row in sorted(reversal_variants.items(),
                                        key=lambda x: x[0]):
            direction = []
            if vid in p_to_b_ids:
                direction.append("P→B")
            if vid in b_to_p_ids:
                direction.append("B→P")

            events = reversal_events.get(vid, [])
            event_str = " | ".join(
                f"{e['Prev_Date']}:{e['Prev_Category']}→"
                f"{e['Curr_Date']}:{e['Curr_Category']}"
                f"[{e['Driver_Label']}]"
                for e in events
            )

            out_row = {
                "Variation_ID": vid,
                "VCV_ID": summary_row.get("VCV_ID", ""),
                "Gene": summary_row.get("Gene", ""),
                "Variant_Type": summary_row.get("Variant_Type", ""),
                "Molecular_Consequence": summary_row.get(
                    "Molecular_Consequence", ""),
                "HGVS_c": summary_row.get("HGVS_c", ""),
                "HGVS_p": summary_row.get("HGVS_p", ""),
                "dbSNP": summary_row.get("dbSNP", ""),
                "Chromosome": summary_row.get("Chromosome", ""),
                "Assembly": summary_row.get("Assembly", ""),
                "First_Classification": summary_row.get(
                    "First_Classification", ""),
                "First_Date": summary_row.get("First_Date", ""),
                "Current_Classification": summary_row.get(
                    "Current_Classification", ""),
                "Current_Category": summary_row.get("Current_Category", ""),
                "Review_Status": summary_row.get("Review_Status", ""),
                "Num_Submissions": summary_row.get("Num_Submissions", ""),
                "Num_PMIDs": summary_row.get("Num_PMIDs", ""),
                "Num_ACMG_Codes": summary_row.get("Num_ACMG_Codes", ""),
                "Change_Type": summary_row.get("Change_Type", ""),
                "Num_Changes": summary_row.get("Num_Changes", ""),
                "Had_Oscillation": summary_row.get("Had_Oscillation", ""),
                "Reversal_Direction": " & ".join(direction),
                "Num_Reclass_Events": len(events),
                "Event_History": event_str,
            }
            writer.writerow(out_row)

    print(f"  Saved: {out_csv}")

    # Step 4: Summary statistics
    gene_counter = Counter(
        reversal_variants[vid].get("Gene", "")
        for vid in reversal_variants
    )
    top_genes = gene_counter.most_common(20)

    # Driver distribution for reversal events
    driver_counter = Counter()
    for vid in reversal_variants:
        for e in reversal_events.get(vid, []):
            # Only count the drastic reversal events themselves
            prev_side = _cat_to_side(e["Prev_Category"])
            curr_side = _cat_to_side(e["Curr_Category"])
            if (prev_side == "P" and curr_side == "B") or \
               (prev_side == "B" and curr_side == "P"):
                driver_counter[e["Driver_Label"]] += 1

    summary_stats = {
        "total_reversal_variants": len(reversal_variants),
        "p_to_b_variants": len(p_to_b_ids),
        "b_to_p_variants": len(b_to_p_ids),
        "both_directions": len(p_to_b_ids & b_to_p_ids),
        "top_20_genes": [{"gene": g, "count": c} for g, c in top_genes],
        "drivers_of_drastic_reversals": dict(driver_counter.most_common()),
    }

    out_json = out_dir / "plp_blb_reversal_summary.json"
    with open(out_json, "w") as f:
        json.dump(summary_stats, f, indent=2)
    print(f"  Saved: {out_json}")

    print(f"\n  ── Top 20 Genes with Reversals ──")
    for g, c in top_genes:
        print(f"    {g:<15s} {c:>5,d}")

    print(f"\n  ── Drivers of Drastic Reversal Events ──")
    for d, c in driver_counter.most_common():
        print(f"    {d:<30s} {c:>6,d}")

    return summary_stats


def _cat_to_side(cat: str) -> str:
    """Map a classification category to P-side, B-side, or other."""
    cat = cat.strip()
    if cat in P_SIDE:
        return "P"
    elif cat in B_SIDE:
        return "B"
    else:
        return "O"


def _normalize_to_side(raw: str) -> str:
    """Normalize raw classification string to P/B/O side."""
    r = raw.strip().lower()
    if "pathogenic" in r and "benign" not in r:
        return "P"
    elif "benign" in r and "pathogenic" not in r:
        return "B"
    else:
        return "O"


# ═════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="ClinVar comment-driven analyses (#16, #17, reversals)")
    parser.add_argument("--summary-csv", type=Path, required=True,
                        help="Path to summary_*.csv")
    parser.add_argument("--events-csv", type=Path, required=True,
                        help="Path to reclassification_events_detail.csv")
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="Output directory (default: current)")
    args = parser.parse_args()

    if not args.summary_csv.exists():
        print(f"ERROR: Summary CSV not found: {args.summary_csv}")
        return 1
    if not args.events_csv.exists():
        print(f"ERROR: Events CSV not found: {args.events_csv}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("ClinVar Comment-Driven Analyses")
    print("=" * 70)
    print(f"  Summary CSV:  {args.summary_csv}")
    print(f"  Events CSV:   {args.events_csv}")
    print(f"  Output dir:   {args.output_dir}")

    # Run all three analyses
    time_results = analysis_time_to_reclass(
        args.summary_csv, args.events_csv, args.output_dir)

    batch_results = analysis_batch_events(
        args.events_csv, args.output_dir)

    reversal_results = analysis_drastic_reversals(
        args.summary_csv, args.events_csv, args.output_dir)

    print("\n" + "=" * 70)
    print("All analyses complete.")
    print("=" * 70)
    print(f"  Outputs saved to: {args.output_dir}")
    print(f"    - time_to_reclassification.json")
    print(f"    - time_to_reclassification_histogram.json")
    print(f"    - batch_event_decomposition.json")
    print(f"    - plp_blb_reversal_variants.csv")
    print(f"    - plp_blb_reversal_summary.json")

    return 0


if __name__ == "__main__":
    sys.exit(main())
