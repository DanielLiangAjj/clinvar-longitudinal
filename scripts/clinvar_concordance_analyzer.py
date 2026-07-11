#!/usr/bin/env python3
"""
ClinVar Inter-Submitter Concordance/Discordance Analyzer (v2 — CSV-based)
=========================================================================

Analyzes inter-submitter concordance and discordance for ClinVar variants
across longitudinal snapshots (2012-2025). Uses the classification_history CSV
produced by clinvar_ftp_only_v5.py, relying on the Review_Status and
Num_Submissions fields to determine concordance.

Concordance logic:
  - Review_Status contains "conflicting" → DISCORDANT
  - Review_Status contains "multiple submitters, no conflicts" → CONCORDANT
  - Review_Status contains "expert panel" or "practice guideline" → CONCORDANT
  - Num_Submissions > 1 but none of the above → EXCLUDED (ambiguous)
  - Num_Submissions <= 1 → SINGLE_SUBMITTER (not applicable)

This approach avoids the per-submission extraction bug in VCV-era checkpoint
files, where .//Interpretation was not searched for per-ClinicalAssertion
classifications.

Input:  classification_history CSV from clinvar_ftp_only_v5.py
Output: concordance_timeseries.json (same schema as v1 for viz compatibility)

Author: Generated for ClinVar Research Project
Date: 2026-02
"""

import sys
import csv
import json
import logging
import argparse
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
from collections import Counter, defaultdict


# ============================================================================
# LOGGING
# ============================================================================

def setup_logging(name: str = "concordance_v2") -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers = []
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(ch)
    return logger


# ============================================================================
# CONCORDANCE CLASSIFICATION
# ============================================================================

def classify_review_status(review_status: str, num_submissions: int) -> str:
    """Classify a variant's concordance state from its Review_Status.

    Returns one of:
        'concordant'        — multiple submitters agree
        'discordant'        — submitters disagree (conflicting interpretations)
        'single_submitter'  — only one submitter, concordance N/A
        'excluded'          — multi-submitter but review status is ambiguous
    """
    if num_submissions <= 1:
        return 'single_submitter'

    rs_lower = review_status.lower().strip()

    # Discordant: conflicting interpretations
    if 'conflicting' in rs_lower:
        return 'discordant'

    # Concordant: multiple submitters agree, expert panel, practice guideline
    if 'multiple submitters, no conflicts' in rs_lower:
        return 'concordant'
    if 'expert panel' in rs_lower:
        return 'concordant'
    if 'practice guideline' in rs_lower:
        return 'concordant'

    # Multi-submitter but classified by "single submitter" review status
    # (can happen when ClinVar aggregates but doesn't update review status)
    # or "no assertion criteria provided" with multiple subs — ambiguous
    return 'excluded'


# ============================================================================
# SNAPSHOT ACCUMULATOR
# ============================================================================

class SnapshotAccumulator:
    """Accumulates concordance stats for a single snapshot date."""

    def __init__(self, date: str):
        self.date = date
        self.total_variants = 0
        self.single_submitter = 0
        self.multi_sub_variants = 0
        self.concordant = 0
        self.discordant = 0
        self.excluded = 0
        self.classification_patterns = Counter()  # for discordant: category string

    def add(self, status: str, category: str):
        self.total_variants += 1
        if status == 'single_submitter':
            self.single_submitter += 1
        elif status == 'concordant':
            self.multi_sub_variants += 1
            self.concordant += 1
        elif status == 'discordant':
            self.multi_sub_variants += 1
            self.discordant += 1
            # Use the ClinVar classification category as the "pattern"
            # For discordant variants, the overall classification is "Conflicting"
            # The actual sub-classifications are not available in the CSV,
            # but the Category field shows the aggregate classification
            self.classification_patterns[category] += 1
        elif status == 'excluded':
            self.multi_sub_variants += 1
            self.excluded += 1

    def to_dict(self) -> Dict[str, Any]:
        if self.multi_sub_variants > 0:
            concordance_rate_pct = round(
                self.concordant / self.multi_sub_variants * 100, 1
            )
            discordance_rate_pct = round(
                self.discordant / self.multi_sub_variants * 100, 1
            )
        else:
            concordance_rate_pct = 0.0
            discordance_rate_pct = 0.0

        # Top discordant classification categories
        top_patterns = [
            {'pattern': pattern, 'count': count}
            for pattern, count in self.classification_patterns.most_common(10)
        ]

        return {
            'date': self.date,
            'total_variants': self.total_variants,
            'single_submitter': self.single_submitter,
            'multi_sub_variants': self.multi_sub_variants,
            'concordant': self.concordant,
            'discordant': self.discordant,
            'excluded': self.excluded,
            'concordance_rate_pct': concordance_rate_pct,
            'discordance_rate_pct': discordance_rate_pct,
            'top_discordance_patterns': top_patterns,
        }


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def run_concordance_analysis(history_csv: Path, output_path: Path,
                              logger: logging.Logger) -> bool:
    """Run concordance analysis by streaming through classification_history CSV.

    Each row represents one variant at one snapshot. We group by Snapshot_Date
    and accumulate concordance statistics.

    The CSV columns used:
        Variation_ID, Snapshot_Date, Category, Review_Status, Num_Submissions
    """
    start_time = time.time()

    logger.info("=" * 70)
    logger.info("ClinVar Inter-Submitter Concordance Analyzer v2 (CSV-based)")
    logger.info("=" * 70)
    logger.info(f"Input CSV: {history_csv}")
    logger.info(f"Output: {output_path}")

    if not history_csv.exists():
        logger.error(f"Input file not found: {history_csv}")
        return False

    # ---- Stream through CSV ----
    logger.info("")
    logger.info("Streaming through classification history CSV...")

    accumulators: Dict[str, SnapshotAccumulator] = {}
    row_count = 0
    error_count = 0

    with open(history_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            row_count += 1

            if row_count % 5_000_000 == 0:
                logger.info(f"  Processed {row_count:,} rows...")

            try:
                snapshot_date = row['Snapshot_Date']
                category = row.get('Category', '')
                review_status = row.get('Review_Status', '')
                num_subs_str = row.get('Num_Submissions', '1')

                try:
                    num_subs = int(num_subs_str)
                except (ValueError, TypeError):
                    num_subs = 1

                # Get or create accumulator for this snapshot
                if snapshot_date not in accumulators:
                    accumulators[snapshot_date] = SnapshotAccumulator(snapshot_date)

                status = classify_review_status(review_status, num_subs)
                accumulators[snapshot_date].add(status, category)

            except Exception as e:
                error_count += 1
                if error_count <= 5:
                    logger.warning(f"  Row {row_count}: {e}")

    logger.info(f"  Total rows processed: {row_count:,}")
    logger.info(f"  Unique snapshots: {len(accumulators)}")
    if error_count > 0:
        logger.warning(f"  Errors: {error_count}")

    if not accumulators:
        logger.error("No data found in CSV!")
        return False

    # ---- Build results sorted by date ----
    logger.info("")
    logger.info("Building output...")

    results = []
    for date in sorted(accumulators.keys()):
        results.append(accumulators[date].to_dict())

    output_data = {
        'generated_at': datetime.now().isoformat(),
        'description': (
            'Inter-submitter concordance analysis across ClinVar longitudinal '
            'snapshots. Uses Review_Status field to determine concordance: '
            '"conflicting interpretations" = discordant, '
            '"multiple submitters, no conflicts" / "expert panel" / '
            '"practice guideline" = concordant.'
        ),
        'method': 'review_status_based',
        'total_rows_processed': row_count,
        'snapshots': results,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2)

    elapsed = time.time() - start_time

    # ---- Summary ----
    logger.info("")
    logger.info("=" * 70)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Snapshots: {len(results)}")
    logger.info(f"Date range: {results[0]['date']} to {results[-1]['date']}")
    logger.info(f"Output: {output_path}")
    logger.info(f"Elapsed: {elapsed:.1f}s")

    logger.info("")
    logger.info("Concordance rate over time (multi-submitter variants):")
    for snap in results:
        if snap['multi_sub_variants'] > 0:
            bar_len = int(snap['concordance_rate_pct'] / 2)
            bar = '#' * bar_len
            logger.info(
                f"  {snap['date']}  conc={snap['concordance_rate_pct']:5.1f}%  "
                f"disc={snap['discordance_rate_pct']:5.1f}%  "
                f"({snap['concordant']:,}/{snap['multi_sub_variants']:,})  {bar}"
            )

    return True


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze inter-submitter concordance/discordance for ClinVar variants "
            "using the classification_history CSV. Uses Review_Status to determine "
            "concordance (v2 approach, avoids VCV checkpoint extraction bug)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s --history-csv classification_history_20260219.csv\n"
            "  %(prog)s --history-csv /path/to/history.csv --output concordance.json\n"
        ),
    )

    parser.add_argument(
        '--history-csv',
        type=str,
        required=True,
        help='Path to classification_history CSV from clinvar_ftp_only_v5.py',
    )
    parser.add_argument(
        '--output',
        type=str,
        default='concordance_timeseries.json',
        help='Path for output JSON (default: concordance_timeseries.json)',
    )

    args = parser.parse_args()

    logger = setup_logging()

    try:
        success = run_concordance_analysis(
            history_csv=Path(args.history_csv),
            output_path=Path(args.output),
            logger=logger,
        )
        return 0 if success else 1

    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 1
    except Exception as e:
        logger.error(f"FATAL ERROR: {e}")
        logger.error(traceback.format_exc())
        return 1


if __name__ == '__main__':
    sys.exit(main())
