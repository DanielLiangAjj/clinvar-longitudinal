#!/usr/bin/env python3
"""
ClinVar Submitter Landscape Analysis
=====================================
Parses ClinVar's submission_summary.txt.gz and organization_summary.txt
to produce submitter-level statistics for the manuscript.

Outputs:
    - submitter_overall.json:   top submitters ranked by total submission count,
                                with institution type, review status, collection method
    - submitter_temporal.json:  yearly submission counts for top N submitters
    - submitter_stats.json:     aggregate statistics (total orgs, submissions by type, etc.)

Data sources (NCBI ClinVar FTP):
    submission_summary.txt.gz   — every current SCV record with submitter name & date
    organization_summary.txt    — org metadata (type, country, submission breakdown)

Usage:
    python3 clinvar_submitter_analysis.py \
        --submission-file /path/to/submission_summary.txt.gz \
        --org-file        /path/to/organization_summary.txt \
        --output-dir      /path/to/output
"""

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


# ── Column indices in submission_summary.txt (0-based, tab-delimited) ──
COL_VARIATION_ID = 0
COL_CLIN_SIG = 1
COL_DATE_EVALUATED = 2
COL_REVIEW_STATUS = 6
COL_COLLECTION_METHOD = 7
COL_SUBMITTER = 9
COL_SCV = 10

# Year parsing helper
MONTH_MAP = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
}


def parse_clinvar_date(date_str: str):
    """Parse ClinVar date format 'Mon DD, YYYY' → year int or None."""
    if not date_str or date_str == "-":
        return None
    parts = date_str.replace(",", "").split()
    if len(parts) == 3 and parts[0] in MONTH_MAP:
        try:
            return int(parts[2])
        except ValueError:
            return None
    # Try 'YYYY-MM-DD' fallback
    try:
        return int(date_str[:4])
    except (ValueError, IndexError):
        return None


def load_org_metadata(org_path: Path) -> dict:
    """Load organization_summary.txt → dict keyed by org name."""
    orgs = {}
    with open(org_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            name = cols[0]
            orgs[name] = {
                "org_id": cols[1],
                "institution_type": cols[2],
                "city": cols[3] if len(cols) > 3 else "",
                "country": cols[5] if len(cols) > 5 else "",
                "total_submissions": int(cols[6]) if cols[6].isdigit() else 0,
                "date_last_submitted": cols[7] if len(cols) > 7 else "",
                "max_review_status": cols[8] if len(cols) > 8 else "",
                "collection_methods": cols[9] if len(cols) > 9 else "",
                "n_clinical_testing": int(cols[12]) if len(cols) > 12 and cols[12].isdigit() else 0,
                "n_research": int(cols[13]) if len(cols) > 13 and cols[13].isdigit() else 0,
                "n_literature": int(cols[14]) if len(cols) > 14 and cols[14].isdigit() else 0,
                "n_curation": int(cols[15]) if len(cols) > 15 and cols[15].isdigit() else 0,
            }
    return orgs


def analyse_submissions(submission_path: Path, org_metadata: dict, output_dir: Path):
    """Main analysis: parse submission_summary.txt.gz and produce statistics."""

    print(f"Parsing {submission_path.name} ...")

    # Accumulators
    submitter_total = Counter()             # submitter → total count
    submitter_by_year = defaultdict(Counter) # year → submitter → count
    submitter_clinsig = defaultdict(Counter)  # submitter → clin_sig → count
    submitter_review = defaultdict(Counter)   # submitter → review_status → count
    submitter_method = defaultdict(Counter)   # submitter → method → count
    submitter_first_year = {}                 # submitter → earliest evaluation year
    year_total = Counter()                    # year → total count

    n_rows = 0
    n_no_date = 0

    opener = gzip.open if str(submission_path).endswith(".gz") else open
    with opener(submission_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= COL_SUBMITTER:
                continue

            submitter = cols[COL_SUBMITTER]
            if not submitter or submitter == "-":
                continue

            n_rows += 1
            submitter_total[submitter] += 1

            # Clinical significance
            clin_sig = cols[COL_CLIN_SIG] if len(cols) > COL_CLIN_SIG else ""
            if clin_sig:
                submitter_clinsig[submitter][clin_sig] += 1

            # Review status
            review = cols[COL_REVIEW_STATUS] if len(cols) > COL_REVIEW_STATUS else ""
            if review:
                submitter_review[submitter][review] += 1

            # Collection method
            method = cols[COL_COLLECTION_METHOD] if len(cols) > COL_COLLECTION_METHOD else ""
            if method:
                submitter_method[submitter][method] += 1

            # Year from DateLastEvaluated
            date_str = cols[COL_DATE_EVALUATED] if len(cols) > COL_DATE_EVALUATED else ""
            year = parse_clinvar_date(date_str)
            if year and 2000 <= year <= 2026:
                submitter_by_year[year][submitter] += 1
                year_total[year] += 1
                if submitter not in submitter_first_year or year < submitter_first_year[submitter]:
                    submitter_first_year[submitter] = year
            else:
                n_no_date += 1

            if n_rows % 1_000_000 == 0:
                print(f"  processed {n_rows:,} rows ...")

    print(f"  Total rows: {n_rows:,} ({n_no_date:,} without valid date)")

    # ── Build overall rankings ────────────────────────────────────────
    top_n = 50
    top_submitters = submitter_total.most_common(top_n)

    overall_list = []
    for name, count in top_submitters:
        meta = org_metadata.get(name, {})
        # Dominant classification
        sig_counts = submitter_clinsig.get(name, {})
        top_sig = max(sig_counts, key=sig_counts.get) if sig_counts else ""
        # Review status distribution
        rev_counts = submitter_review.get(name, {})
        top_review = max(rev_counts, key=rev_counts.get) if rev_counts else ""
        # Method
        meth_counts = submitter_method.get(name, {})
        top_method = max(meth_counts, key=meth_counts.get) if meth_counts else ""

        overall_list.append({
            "rank": len(overall_list) + 1,
            "submitter": name,
            "total_submissions": count,
            "pct_of_total": round(100 * count / n_rows, 2),
            "institution_type": meta.get("institution_type", "unknown"),
            "country": meta.get("country", ""),
            "dominant_classification": top_sig,
            "dominant_review_status": top_review,
            "dominant_method": top_method,
            "first_evaluation_year": submitter_first_year.get(name),
            "classification_breakdown": dict(sig_counts),
        })

    # ── Build temporal data for top 20 ────────────────────────────────
    top20_names = [x["submitter"] for x in overall_list[:20]]
    years = sorted(year_total.keys())

    temporal = {}
    for name in top20_names:
        temporal[name] = {
            str(y): submitter_by_year[y].get(name, 0) for y in years
        }

    # Add "Other" aggregation
    other_by_year = {}
    for y in years:
        top20_sum = sum(submitter_by_year[y].get(n, 0) for n in top20_names)
        other_by_year[str(y)] = year_total[y] - top20_sum
    temporal["_Other"] = other_by_year

    # ── Build aggregate stats ─────────────────────────────────────────
    # Institution type distribution
    type_counter = Counter()
    type_submission_counter = Counter()
    for name, count in submitter_total.items():
        meta = org_metadata.get(name, {})
        itype = meta.get("institution_type", "unknown")
        type_counter[itype] += 1
        type_submission_counter[itype] += count

    # Collection method distribution
    method_total = Counter()
    for name, mc in submitter_method.items():
        for meth, cnt in mc.items():
            method_total[meth] += cnt

    # Review status distribution
    review_total = Counter()
    for name, rc in submitter_review.items():
        for rev, cnt in rc.items():
            review_total[rev] += cnt

    # Concentration: top 10 submitters as % of all submissions
    top10_sum = sum(c for _, c in submitter_total.most_common(10))
    top20_sum = sum(c for _, c in submitter_total.most_common(20))

    stats_dict = {
        "total_submissions": n_rows,
        "total_unique_submitters": len(submitter_total),
        "submissions_without_date": n_no_date,
        "top10_concentration_pct": round(100 * top10_sum / n_rows, 2),
        "top20_concentration_pct": round(100 * top20_sum / n_rows, 2),
        "institution_type_orgs": dict(type_counter.most_common()),
        "institution_type_submissions": dict(type_submission_counter.most_common()),
        "collection_method_totals": dict(method_total.most_common()),
        "review_status_totals": dict(review_total.most_common()),
        "submissions_per_year": {str(y): year_total[y] for y in years},
    }

    # ── Write outputs ─────────────────────────────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(output_dir / "submitter_overall.json", "w") as f:
        json.dump(overall_list, f, indent=2)
    print(f"  Wrote submitter_overall.json ({len(overall_list)} entries)")

    with open(output_dir / "submitter_temporal.json", "w") as f:
        json.dump({"years": [str(y) for y in years], "submitters": temporal}, f, indent=2)
    print(f"  Wrote submitter_temporal.json ({len(temporal)} submitters × {len(years)} years)")

    with open(output_dir / "submitter_stats.json", "w") as f:
        json.dump(stats_dict, f, indent=2)
    print(f"  Wrote submitter_stats.json")

    return overall_list, temporal, stats_dict


def main():
    parser = argparse.ArgumentParser(
        description="ClinVar Submitter Landscape Analysis"
    )
    parser.add_argument("--submission-file", required=True,
                        help="Path to submission_summary.txt.gz from ClinVar FTP")
    parser.add_argument("--org-file", required=True,
                        help="Path to organization_summary.txt from ClinVar FTP")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for output JSON files")
    args = parser.parse_args()

    submission_path = Path(args.submission_file).expanduser()
    org_path = Path(args.org_file).expanduser()
    output_dir = Path(args.output_dir).expanduser()

    if not submission_path.exists():
        print(f"ERROR: {submission_path} not found")
        sys.exit(1)
    if not org_path.exists():
        print(f"ERROR: {org_path} not found")
        sys.exit(1)

    print("=" * 70)
    print("ClinVar Submitter Landscape Analysis")
    print("=" * 70)

    org_metadata = load_org_metadata(org_path)
    print(f"Loaded metadata for {len(org_metadata)} organizations")

    analyse_submissions(submission_path, org_metadata, output_dir)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
