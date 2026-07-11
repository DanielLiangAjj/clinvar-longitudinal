#!/usr/bin/env python3
"""
ClinVar Reclassification Driver Classifier
===========================================

Attempts to distinguish schema-driven reclassifications (artifacts of the
RCV→VCV data model transition) from knowledge-driven reclassifications
(genuine updates to variant interpretation based on new evidence).

Three signals are used:
  1. ERA BOUNDARIES — reclassifications at the RCV→VCV transition (2019-12
     to 2020-01) or the early VCV stabilization window (2020-01 to 2020-06)
     are flagged as likely schema-driven.
  2. NUM_SUBMISSIONS CHANGE — if num_submissions changes at the same snapshot
     as the reclassification, it suggests re-aggregation (schema-driven).
     Stable num_submissions suggests genuine reinterpretation.
  3. REVIEW_STATUS CHANGE — if review_status is lost or gains "conflicting"
     at the transition, it suggests re-aggregation. Stable review_status
     suggests knowledge-driven change.

Each reclassification event is labeled:
  - likely_knowledge_driven — no era boundary signal, stable submissions
  - likely_schema_driven   — at era boundary + submission/review change
  - ambiguous              — mixed signals
  - batch_event            — >5% of variants in that snapshot reclassified
                             (suggests systematic data processing change)

Input:
  - classification_history CSV (per-variant longitudinal data)
  - summary CSV (optional, for enrichment with gene/variant info)

Output:
  - reclassification_drivers.json — aggregate statistics
  - reclassification_events_detail.csv — per-event detail (optional)

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
from typing import Dict, List, Optional, Any, Tuple
from collections import Counter, defaultdict


# ============================================================================
# CONSTANTS
# ============================================================================

# Schema eras
ERA_RCV_END = '2019-12'       # Last pure RCV snapshot
ERA_VCV_START = '2020-01'     # First VCV snapshot
ERA_VCV_STABLE = '2020-07'    # VCV format stabilized

# Transition window: snapshots where schema-driven artifacts are most likely
TRANSITION_SNAPSHOTS = {
    '2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06'
}

# The single most critical boundary
PRIMARY_BOUNDARY = ('2019-12', '2020-01')

# Batch threshold: if >5% of tracked variants reclassify in one snapshot,
# flag the entire snapshot as a batch event
BATCH_THRESHOLD_PCT = 5.0


# ============================================================================
# LOGGING
# ============================================================================

def setup_logging(name: str = "reclass_driver") -> logging.Logger:
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
# DRIVER CLASSIFICATION
# ============================================================================

def classify_driver(
    prev_date: str,
    curr_date: str,
    prev_category: str,
    curr_category: str,
    prev_num_subs: int,
    curr_num_subs: int,
    prev_review_status: str,
    curr_review_status: str,
    is_batch_snapshot: bool,
) -> Tuple[str, List[str]]:
    """Classify a single reclassification event.

    Strategy:
      - BATCH EVENTS: snapshots with >5% simultaneous reclassification
        override all other signals.
      - TRANSITION WINDOW (2020-01 to 2020-06): The only period where
        schema-driven vs knowledge-driven can be meaningfully distinguished
        using era boundary, submission count, and review status signals.
      - OUTSIDE TRANSITION: We classify by *observable mechanism* rather
        than pretending to know the underlying cause. Categories:
          new_submitter — num_submissions increased (new lab contributed)
          reinterpretation — num_submissions stable, no conflict status
                            change (existing submitter updated)
          conflict_emerged — gained 'conflicting' review status
          conflict_resolved — lost 'conflicting' review status

    Returns:
        (label, signals) where label is one of:
            'likely_schema_driven', 'likely_knowledge_driven',
            'transition_ambiguous',
            'new_submitter', 'reinterpretation',
            'conflict_emerged', 'conflict_resolved',
            'batch_event'
        and signals is a list of reason strings.
    """
    signals = []

    # ---- Batch event override (checked first) ----
    if is_batch_snapshot:
        return ('batch_event', ['batch_snapshot'])

    # ---- Derived features ----
    at_primary_boundary = (prev_date == PRIMARY_BOUNDARY[0] and
                           curr_date == PRIMARY_BOUNDARY[1])
    at_transition_window = curr_date in TRANSITION_SNAPSHOTS
    subs_changed = (prev_num_subs != curr_num_subs)
    subs_increased = (curr_num_subs > prev_num_subs)

    prev_rs = prev_review_status.lower().strip()
    curr_rs = curr_review_status.lower().strip()
    review_status_changed = (prev_rs != curr_rs)
    review_status_lost = (prev_rs and not curr_rs)
    gained_conflicting = ('conflicting' not in prev_rs and
                          'conflicting' in curr_rs)
    lost_conflicting = ('conflicting' in prev_rs and
                        'conflicting' not in curr_rs)

    # ==================================================================
    # PATH A: TRANSITION WINDOW (2020-01 to 2020-06)
    # Schema vs knowledge scoring only meaningful here.
    # ==================================================================
    if at_transition_window:
        schema_signals = 0
        knowledge_signals = 0

        # Signal 1: Era boundary
        if at_primary_boundary:
            schema_signals += 2
            signals.append('at_rcv_vcv_boundary')
        else:
            schema_signals += 1
            signals.append('in_transition_window')

        # Signal 2: Submission count
        if subs_changed:
            schema_signals += 1
            signals.append('subs_changed_at_transition')
        else:
            knowledge_signals += 1
            signals.append('subs_stable_at_transition')

        # Signal 3: Review status
        if review_status_lost:
            schema_signals += 1
            signals.append('review_status_lost_at_transition')
        elif gained_conflicting:
            schema_signals += 1
            signals.append('gained_conflicting_at_transition')
        elif lost_conflicting:
            knowledge_signals += 1
            signals.append('resolved_conflict_at_transition')
        elif not review_status_changed:
            knowledge_signals += 1
            signals.append('review_status_stable_at_transition')
        else:
            signals.append('review_status_changed_at_transition')

        # Signal 4: "Other" category at transition = schema artifact
        if curr_category == 'Other':
            schema_signals += 1
            signals.append('to_other_at_transition')

        # Final transition classification
        if schema_signals >= 3:
            return ('likely_schema_driven', signals)
        elif schema_signals >= 2 and knowledge_signals == 0:
            return ('likely_schema_driven', signals)
        elif knowledge_signals >= 2 and schema_signals == 0:
            return ('likely_knowledge_driven', signals)
        elif knowledge_signals > schema_signals:
            return ('likely_knowledge_driven', signals)
        elif schema_signals > knowledge_signals:
            return ('likely_schema_driven', signals)
        else:
            return ('transition_ambiguous', signals)

    # ==================================================================
    # PATH B: OUTSIDE TRANSITION (RCV era + VCV-stable era)
    # Classify by observable mechanism — no schema/knowledge inference.
    # ==================================================================
    signals.append('outside_transition')

    # Observable mechanism classification
    if gained_conflicting:
        signals.append('gained_conflicting')
        return ('conflict_emerged', signals)

    if lost_conflicting:
        signals.append('lost_conflicting')
        return ('conflict_resolved', signals)

    if subs_increased:
        signals.append('subs_increased')
        return ('new_submitter', signals)

    # Stable or decreased submissions, no conflict change
    if subs_changed:
        signals.append('subs_decreased')
    else:
        signals.append('subs_stable')
    if review_status_changed:
        signals.append('review_status_changed')
    else:
        signals.append('review_status_stable')

    return ('reinterpretation', signals)


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def run_driver_analysis(history_csv: Path, output_dir: Path,
                         write_detail: bool,
                         logger: logging.Logger) -> bool:
    """Stream through classification_history CSV and classify each
    reclassification event.

    For each variant, we track consecutive snapshots and detect when the
    Category field changes. At each change, we collect the surrounding
    context (num_submissions, review_status, dates) and classify the driver.
    """
    start_time = time.time()

    logger.info("=" * 70)
    logger.info("ClinVar Reclassification Driver Classifier")
    logger.info("=" * 70)
    logger.info(f"Input: {history_csv}")
    logger.info(f"Output dir: {output_dir}")

    if not history_csv.exists():
        logger.error(f"Input not found: {history_csv}")
        return False

    # ---- Pass 1: Count variants per snapshot (for batch detection) ----
    logger.info("")
    logger.info("Pass 1: Counting variants per snapshot for batch detection...")

    snapshot_total = Counter()     # total variants per snapshot
    snapshot_changed = Counter()   # reclassified variants per snapshot

    # Track previous state per variant
    prev_state: Dict[str, Dict] = {}  # variation_id -> {category, num_subs, review_status, date}
    reclassification_events = []       # collected for pass 1

    row_count = 0
    with open(history_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            row_count += 1
            if row_count % 5_000_000 == 0:
                logger.info(f"  {row_count:,} rows...")

            var_id = row['Variation_ID']
            snapshot_date = row['Snapshot_Date']
            category = row.get('Category', '')
            review_status = row.get('Review_Status', '')
            try:
                num_subs = int(row.get('Num_Submissions', '1'))
            except (ValueError, TypeError):
                num_subs = 1

            snapshot_total[snapshot_date] += 1

            if var_id in prev_state:
                prev = prev_state[var_id]
                if category and prev['category'] and category != prev['category']:
                    snapshot_changed[snapshot_date] += 1
                    reclassification_events.append({
                        'variation_id': var_id,
                        'gene': row.get('Gene', ''),
                        'prev_date': prev['date'],
                        'curr_date': snapshot_date,
                        'prev_category': prev['category'],
                        'curr_category': category,
                        'prev_num_subs': prev['num_subs'],
                        'curr_num_subs': num_subs,
                        'prev_review_status': prev['review_status'],
                        'curr_review_status': review_status,
                    })

            prev_state[var_id] = {
                'category': category,
                'num_subs': num_subs,
                'review_status': review_status,
                'date': snapshot_date,
            }

    logger.info(f"  Total rows: {row_count:,}")
    logger.info(f"  Reclassification events detected: {len(reclassification_events):,}")

    # Free memory
    del prev_state

    # ---- Identify batch snapshots ----
    batch_snapshots = set()
    for date in sorted(snapshot_total.keys()):
        total = snapshot_total[date]
        changed = snapshot_changed.get(date, 0)
        if total > 0:
            pct = changed / total * 100
            if pct > BATCH_THRESHOLD_PCT:
                batch_snapshots.add(date)
                logger.info(f"  BATCH snapshot: {date} — {changed:,}/{total:,} = {pct:.1f}%")

    # ---- Pass 2: Classify each reclassification event ----
    logger.info("")
    logger.info("Pass 2: Classifying reclassification events...")

    driver_counts = Counter()       # label -> count
    era_driver_counts = defaultdict(Counter)  # era -> label -> count
    per_snapshot_drivers = defaultdict(Counter)  # date -> label -> count
    change_type_by_driver = defaultdict(Counter)  # label -> change_type -> count
    signal_counts = Counter()       # individual signals
    detail_rows = []

    for event in reclassification_events:
        is_batch = event['curr_date'] in batch_snapshots

        label, signals = classify_driver(
            prev_date=event['prev_date'],
            curr_date=event['curr_date'],
            prev_category=event['prev_category'],
            curr_category=event['curr_category'],
            prev_num_subs=event['prev_num_subs'],
            curr_num_subs=event['curr_num_subs'],
            prev_review_status=event['prev_review_status'],
            curr_review_status=event['curr_review_status'],
            is_batch_snapshot=is_batch,
        )

        driver_counts[label] += 1

        # Determine era
        curr_date = event['curr_date']
        if curr_date <= ERA_RCV_END:
            era = 'RCV (2012-2019)'
        elif curr_date <= '2020-06':
            era = 'Transition (2020-01 to 2020-06)'
        else:
            era = 'VCV-stable (2020-07+)'
        era_driver_counts[era][label] += 1

        per_snapshot_drivers[curr_date][label] += 1

        change_type = f"{event['prev_category']}->{event['curr_category']}"
        change_type_by_driver[label][change_type] += 1

        for sig in signals:
            signal_counts[sig] += 1

        if write_detail:
            detail_rows.append({
                'Variation_ID': event['variation_id'],
                'Gene': event['gene'],
                'Prev_Date': event['prev_date'],
                'Curr_Date': event['curr_date'],
                'Prev_Category': event['prev_category'],
                'Curr_Category': event['curr_category'],
                'Prev_Num_Subs': event['prev_num_subs'],
                'Curr_Num_Subs': event['curr_num_subs'],
                'Driver_Label': label,
                'Signals': '|'.join(signals),
            })

    # ---- Build output JSON ----
    ALL_LABELS = [
        'likely_schema_driven', 'likely_knowledge_driven', 'transition_ambiguous',
        'new_submitter', 'reinterpretation', 'conflict_emerged', 'conflict_resolved',
        'batch_event',
    ]
    # Higher-level groupings for summary
    TRANSITION_LABELS = {'likely_schema_driven', 'likely_knowledge_driven', 'transition_ambiguous'}
    MECHANISM_LABELS = {'new_submitter', 'reinterpretation', 'conflict_emerged', 'conflict_resolved'}

    total_events = sum(driver_counts.values())
    logger.info(f"  Total reclassification events: {total_events:,}")
    for label in ALL_LABELS:
        count = driver_counts.get(label, 0)
        if count > 0:
            pct = count / total_events * 100
            logger.info(f"    {label}: {count:,} ({pct:.1f}%)")

    # Per-snapshot time series
    snapshot_timeseries = []
    for date in sorted(per_snapshot_drivers.keys()):
        entry = {'date': date}
        for label in ALL_LABELS:
            entry[label] = per_snapshot_drivers[date].get(label, 0)
        entry['total_reclassifications'] = sum(per_snapshot_drivers[date].values())
        snapshot_timeseries.append(entry)

    # Per-era summary
    era_summary = {}
    for era in ['RCV (2012-2019)', 'Transition (2020-01 to 2020-06)', 'VCV-stable (2020-07+)']:
        era_total = sum(era_driver_counts[era].values())
        era_entry = {'total_events': era_total}
        for label in ALL_LABELS:
            era_entry[label] = era_driver_counts[era].get(label, 0)
        if era_total > 0:
            for label in ALL_LABELS:
                era_entry[f'{label}_pct'] = round(
                    era_entry[label] / era_total * 100, 1
                )
        # Also add grouped summaries
        era_entry['transition_total'] = sum(
            era_driver_counts[era].get(l, 0) for l in TRANSITION_LABELS)
        era_entry['mechanism_total'] = sum(
            era_driver_counts[era].get(l, 0) for l in MECHANISM_LABELS)
        era_summary[era] = era_entry

    # Top change types per driver
    top_changes_by_driver = {}
    for label, counter in change_type_by_driver.items():
        top_changes_by_driver[label] = [
            {'change_type': ct, 'count': c}
            for ct, c in counter.most_common(15)
        ]

    output_data = {
        'generated_at': datetime.now().isoformat(),
        'description': (
            'Reclassification driver classification (v2). '
            'Transition window (2020-01 to 2020-06): schema-driven vs '
            'knowledge-driven using era boundaries, submission count, and '
            'review status signals. '
            'Outside transition: classified by observable mechanism '
            '(new_submitter, reinterpretation, conflict_emerged, conflict_resolved). '
            'Batch events: snapshots with >5% simultaneous reclassification.'
        ),
        'total_reclassification_events': total_events,
        'all_labels': ALL_LABELS,
        'driver_summary': {
            label: {
                'count': driver_counts.get(label, 0),
                'pct': round(driver_counts.get(label, 0) / total_events * 100, 1)
                       if total_events > 0 else 0.0
            }
            for label in ALL_LABELS
        },
        'grouped_summary': {
            'transition_schema_driven': driver_counts.get('likely_schema_driven', 0),
            'transition_knowledge_driven': driver_counts.get('likely_knowledge_driven', 0),
            'transition_ambiguous': driver_counts.get('transition_ambiguous', 0),
            'mechanism_new_submitter': driver_counts.get('new_submitter', 0),
            'mechanism_reinterpretation': driver_counts.get('reinterpretation', 0),
            'mechanism_conflict_emerged': driver_counts.get('conflict_emerged', 0),
            'mechanism_conflict_resolved': driver_counts.get('conflict_resolved', 0),
            'batch_event': driver_counts.get('batch_event', 0),
        },
        'era_summary': era_summary,
        'signal_frequency': dict(signal_counts.most_common()),
        'top_changes_by_driver': top_changes_by_driver,
        'batch_snapshots': sorted(batch_snapshots),
        'timeseries': snapshot_timeseries,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / 'reclassification_drivers.json'
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2)
    logger.info(f"  Wrote: {json_path}")

    # ---- Optional detail CSV ----
    if write_detail and detail_rows:
        detail_path = output_dir / 'reclassification_events_detail.csv'
        fieldnames = [
            'Variation_ID', 'Gene', 'Prev_Date', 'Curr_Date',
            'Prev_Category', 'Curr_Category', 'Prev_Num_Subs',
            'Curr_Num_Subs', 'Driver_Label', 'Signals'
        ]
        with open(detail_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(detail_rows)
        logger.info(f"  Wrote detail CSV: {detail_path} ({len(detail_rows):,} rows)")

    elapsed = time.time() - start_time
    logger.info("")
    logger.info("=" * 70)
    logger.info(f"DONE in {elapsed:.1f}s")
    logger.info("=" * 70)

    return True


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Classify ClinVar reclassification events as schema-driven vs "
            "knowledge-driven using era boundaries, submission count changes, "
            "and review status changes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s --history-csv classification_history_20260219.csv\n"
            "  %(prog)s --history-csv history.csv --output-dir results/ --detail\n"
        ),
    )

    parser.add_argument(
        '--history-csv', type=str, required=True,
        help='Path to classification_history CSV from clinvar_ftp_only_v5.py',
    )
    parser.add_argument(
        '--output-dir', type=str, default='.',
        help='Directory for output files (default: current dir)',
    )
    parser.add_argument(
        '--detail', action='store_true',
        help='Write per-event detail CSV (can be large)',
    )

    args = parser.parse_args()
    logger = setup_logging()

    try:
        success = run_driver_analysis(
            history_csv=Path(args.history_csv),
            output_dir=Path(args.output_dir),
            write_detail=args.detail,
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
