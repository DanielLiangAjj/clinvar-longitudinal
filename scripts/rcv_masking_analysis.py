#!/usr/bin/env python3
"""
RCV Masking Analysis
====================

Quantifies the impact of the "highest-weight review status" selection used
in the RCV-era parser. When a variant appears in multiple RCV records (one
per condition), the parser selects the record with the highest-weight review
status. This can mask condition-specific conflicts.

This script parses one or more RCV-era XML snapshots WITHOUT collapsing to
the highest-weight record, and reports:
  1. How many variants have multiple RCV records (multiple conditions)
  2. How many have different review statuses across conditions
  3. How many have different classifications across conditions
  4. How many would have a conflict masked by the highest-weight selection
  5. The impact on reclassification detection

Usage:
    # Analyze the last RCV snapshot (December 2019)
    python rcv_masking_analysis.py --snapshot-dir /path/to/ftp_snapshots \
                                    --output /path/to/output/

    # Analyze a specific snapshot
    python rcv_masking_analysis.py --snapshot /path/to/ClinVarFullRelease_2019-12.xml.gz \
                                    --output /path/to/output/

    # Analyze multiple snapshots to show trend over time
    python rcv_masking_analysis.py --snapshot-dir /path/to/ftp_snapshots \
                                    --all-rcv \
                                    --output /path/to/output/

Author: ClinVar Research Project
Date: 2026-03
"""

import sys
import gzip
import json
import csv
import re
import time
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict, Counter
from datetime import datetime


# ============================================================================
# REVIEW STATUS WEIGHT (same as clinvar_ftp_only_v5.py)
# ============================================================================

def _review_weight(review_status: str) -> int:
    """Convert review status to weight for sorting."""
    status = (review_status or '').lower()
    if 'practice guideline' in status:
        return 4
    elif 'expert panel' in status:
        return 3
    elif 'multiple submitters' in status:
        return 2
    elif 'criteria provided' in status:
        return 1
    return 0


def _normalize_classification(classification: str) -> str:
    """Normalize classification to standard categories."""
    if not classification:
        return 'Other'
    c = classification.lower().strip()
    if 'pathogenic/likely' in c or 'pathogenic, likely' in c:
        return 'P/LP'
    if 'benign/likely' in c or 'benign, likely' in c:
        return 'B/LB'
    if 'likely pathogenic' in c:
        return 'LP'
    if 'pathogenic' in c:
        return 'P'
    if 'likely benign' in c:
        return 'LB'
    if 'benign' in c:
        return 'B'
    if 'uncertain' in c:
        return 'VUS'
    if 'conflicting' in c:
        return 'Conflicting'
    return 'Other'


# ============================================================================
# RCV RECORD EXTRACTION (lightweight — only fields needed for masking analysis)
# ============================================================================

def _extract_rcv_record_light(clinvar_set_elem) -> Optional[Dict]:
    """Extract key fields from a ClinVarSet element for masking analysis.

    Returns a lightweight record with only the fields needed:
      - variation_id, rcv_id, classification, review_status, condition,
        num_submissions, gene_symbol
    """
    ref_assertion = clinvar_set_elem.find('.//ReferenceClinVarAssertion')
    if ref_assertion is None:
        return None

    # Get variation ID
    measure_set = ref_assertion.find('MeasureSet')
    genotype_set = None
    if measure_set is None:
        genotype_set = ref_assertion.find('GenotypeSet')
        if genotype_set is None:
            return None
        variation_id = genotype_set.get('ID')
    else:
        variation_id = measure_set.get('ID')

    if not variation_id:
        return None

    record = {
        'variation_id': variation_id,
        'rcv_id': ref_assertion.get('ID', ''),
    }

    # Classification and review status (from ReferenceClinVarAssertion)
    cls_elem = ref_assertion.find('ClinicalSignificance')
    if cls_elem is not None:
        desc = cls_elem.find('Description')
        if desc is not None and desc.text:
            record['classification'] = desc.text.strip()
        review = cls_elem.find('ReviewStatus')
        if review is not None and review.text:
            record['review_status'] = review.text.strip()

    # Condition (from TraitSet)
    trait_set = ref_assertion.find('.//TraitSet')
    if trait_set is not None:
        trait = trait_set.find('.//Trait')
        if trait is not None:
            name_elem = trait.find('Name/ElementValue')
            if name_elem is not None and name_elem.text:
                record['condition'] = name_elem.text.strip()

    # Gene symbol (first found)
    if measure_set is not None:
        measures = measure_set.findall('.//Measure')
    elif genotype_set is not None:
        measures = genotype_set.findall('.//MeasureSet//Measure')
    else:
        measures = []

    for measure in measures:
        for gene in measure.findall('.//MeasureRelationship/Gene'):
            symbol = gene.get('Symbol', '')
            if symbol:
                record['gene_symbol'] = symbol
                break
        if 'gene_symbol' in record:
            break

    # Count individual submissions (ClinVarAssertion elements)
    submissions = clinvar_set_elem.findall('.//ClinVarAssertion')
    record['num_submissions'] = len(submissions)

    # Extract per-submission classifications for within-condition concordance
    sub_classifications = []
    for assertion in submissions:
        sub_cls = assertion.find('.//ClinicalSignificance/Description')
        if sub_cls is not None and sub_cls.text:
            sub_classifications.append(sub_cls.text.strip())
    record['submission_classifications'] = sub_classifications

    return record


# ============================================================================
# SINGLE SNAPSHOT ANALYSIS
# ============================================================================

def analyze_snapshot(snapshot_path: Path) -> Dict:
    """Parse a single RCV snapshot and analyze cross-condition masking.

    Returns a dict with all analysis results for this snapshot.
    """
    print(f"\n{'='*70}")
    print(f"Analyzing: {snapshot_path.name}")
    print(f"{'='*70}")
    t0 = time.time()

    # Phase 1: Parse all RCV records, grouped by variation_id
    variant_records = defaultdict(list)  # variation_id -> list of RCV records
    total_rcv_records = 0
    parse_errors = 0

    print("Phase 1: Parsing RCV records ...")
    with gzip.open(snapshot_path, 'rb') as f:
        context = ET.iterparse(f, events=('start', 'end'))
        root = None

        for event, elem in context:
            if event == 'start':
                if root is None:
                    root = elem
                continue

            if elem.tag == 'ClinVarSet':
                try:
                    record = _extract_rcv_record_light(elem)
                    if record and record.get('variation_id'):
                        variant_records[record['variation_id']].append(record)
                        total_rcv_records += 1

                        if total_rcv_records % 100_000 == 0:
                            elapsed = time.time() - t0
                            rate = total_rcv_records / elapsed
                            print(f"  Parsed {total_rcv_records:,} RCV records "
                                  f"({rate:.0f}/s, {elapsed:.0f}s)")
                except Exception as e:
                    parse_errors += 1

                if root is not None:
                    root.clear()

    parse_time = time.time() - t0
    total_variants = len(variant_records)
    print(f"  Done: {total_rcv_records:,} RCV records -> "
          f"{total_variants:,} unique variants in {parse_time:.1f}s")
    if parse_errors:
        print(f"  Parse errors: {parse_errors:,}")

    # Phase 2: Analyze masking
    print("\nPhase 2: Analyzing cross-condition masking ...")

    # Counters
    single_rcv = 0          # Variants with only one RCV record
    multi_rcv = 0            # Variants with multiple RCV records
    multi_same_class = 0     # Multi-RCV, all same classification
    multi_diff_class = 0     # Multi-RCV, different classifications across conditions
    multi_same_review = 0    # Multi-RCV, all same review status
    multi_diff_review = 0    # Multi-RCV, different review statuses

    # The key masking metrics
    masking_review_upgrades = 0   # Higher-weight status hides lower-weight status
    masking_hides_conflict = 0    # Specifically: concordant hides conflicting
    masking_changes_class = 0     # Selected classification differs from lower-weight
    masking_concordant_hides_discordant = 0  # Within-condition: one concordant, one discordant

    # Detail tracking
    masked_examples = []
    classification_impact = Counter()  # (selected_cat, masked_cat) -> count

    for var_id, records in variant_records.items():
        if len(records) == 1:
            single_rcv += 1
            continue

        multi_rcv += 1

        # Normalize classifications and review statuses
        for r in records:
            r['norm_class'] = _normalize_classification(
                r.get('classification', ''))
            r['weight'] = _review_weight(r.get('review_status', ''))

        # Check classification uniformity
        classes = set(r['norm_class'] for r in records)
        if len(classes) <= 1:
            multi_same_class += 1
        else:
            multi_diff_class += 1

        # Check review status uniformity
        statuses = set(r.get('review_status', '').lower() for r in records)
        if len(statuses) <= 1:
            multi_same_review += 1
        else:
            multi_diff_review += 1

        # Sort by weight (descending) — same as the parser does
        sorted_records = sorted(records, key=lambda x: (
            x['weight'], x.get('date_updated', '')), reverse=True)
        selected = sorted_records[0]
        others = sorted_records[1:]

        # Check if selection masks a lower review status
        for other in others:
            if other['weight'] < selected['weight']:
                masking_review_upgrades += 1

                # Does it specifically hide a "conflicting" status?
                other_rs = other.get('review_status', '').lower()
                selected_rs = selected.get('review_status', '').lower()
                if ('conflicting' in other_rs and
                        'conflicting' not in selected_rs):
                    masking_hides_conflict += 1

                # Does the classification differ?
                if other['norm_class'] != selected['norm_class']:
                    masking_changes_class += 1
                    classification_impact[
                        (selected['norm_class'], other['norm_class'])] += 1

                    if len(masked_examples) < 20:
                        masked_examples.append({
                            'variation_id': var_id,
                            'gene': selected.get('gene_symbol', ''),
                            'selected_rcv': selected.get('rcv_id', ''),
                            'selected_condition': selected.get(
                                'condition', '?')[:80],
                            'selected_classification': selected.get(
                                'classification', ''),
                            'selected_review_status': selected.get(
                                'review_status', ''),
                            'selected_weight': selected['weight'],
                            'masked_rcv': other.get('rcv_id', ''),
                            'masked_condition': other.get(
                                'condition', '?')[:80],
                            'masked_classification': other.get(
                                'classification', ''),
                            'masked_review_status': other.get(
                                'review_status', ''),
                            'masked_weight': other['weight'],
                        })
                break  # Only count once per variant (first masked record)

        # Check within-condition concordance mismatch
        concordant_conditions = []
        discordant_conditions = []
        for r in records:
            sub_classes = r.get('submission_classifications', [])
            if len(sub_classes) < 2:
                continue
            norm_subs = set(_normalize_classification(c) for c in sub_classes)
            norm_subs -= {'Other'}
            if len(norm_subs) <= 1:
                concordant_conditions.append(r)
            else:
                discordant_conditions.append(r)

        if concordant_conditions and discordant_conditions:
            masking_concordant_hides_discordant += 1

    analysis_time = time.time() - t0 - parse_time

    # Phase 3: Compile results
    print("\nPhase 3: Compiling results ...")

    results = {
        'snapshot': snapshot_path.name,
        'parse_time_seconds': round(parse_time, 1),
        'analysis_time_seconds': round(analysis_time, 1),
        'total_rcv_records': total_rcv_records,
        'total_variants': total_variants,
        'parse_errors': parse_errors,

        'single_rcv_variants': single_rcv,
        'multi_rcv_variants': multi_rcv,
        'multi_rcv_pct': round(multi_rcv / max(total_variants, 1) * 100, 2),

        'multi_same_classification': multi_same_class,
        'multi_diff_classification': multi_diff_class,
        'multi_same_review_status': multi_same_review,
        'multi_diff_review_status': multi_diff_review,

        'masking_review_upgrades': masking_review_upgrades,
        'masking_hides_conflict': masking_hides_conflict,
        'masking_changes_classification': masking_changes_class,
        'masking_concordant_hides_discordant':
            masking_concordant_hides_discordant,

        'masking_review_upgrades_pct': round(
            masking_review_upgrades / max(total_variants, 1) * 100, 3),
        'masking_hides_conflict_pct': round(
            masking_hides_conflict / max(total_variants, 1) * 100, 3),
        'masking_changes_classification_pct': round(
            masking_changes_class / max(total_variants, 1) * 100, 3),
        'masking_concordant_hides_discordant_pct': round(
            masking_concordant_hides_discordant / max(total_variants, 1) * 100,
            3),

        'classification_impact': {
            f"{sel}->{masked}": count
            for (sel, masked), count in classification_impact.most_common(30)
        },
        'masked_examples': masked_examples,
    }

    # Print summary
    _print_summary(results)

    return results


def _print_summary(results: Dict):
    """Print a formatted summary of masking analysis results."""
    print()
    print("=" * 70)
    print(f"RCV MASKING ANALYSIS RESULTS: {results['snapshot']}")
    print("=" * 70)

    tv = results['total_variants']
    mr = results['multi_rcv_variants']

    print(f"\n{'OVERVIEW':}")
    print(f"  Total RCV records parsed:       {results['total_rcv_records']:>12,}")
    print(f"  Unique variants:                {tv:>12,}")
    print(f"  Single-RCV variants:            {results['single_rcv_variants']:>12,}")
    print(f"  Multi-RCV variants:             {mr:>12,}  "
          f"({results['multi_rcv_pct']:.1f}%)")

    print(f"\n{'MULTI-RCV BREAKDOWN':}")
    print(f"  Same classification across all: "
          f"{results['multi_same_classification']:>12,}  "
          f"({results['multi_same_classification']/max(mr,1)*100:.1f}%"
          f" of multi-RCV)")
    print(f"  Diff classification across:     "
          f"{results['multi_diff_classification']:>12,}  "
          f"({results['multi_diff_classification']/max(mr,1)*100:.1f}%"
          f" of multi-RCV)")
    print(f"  Same review status across all:  "
          f"{results['multi_same_review_status']:>12,}  "
          f"({results['multi_same_review_status']/max(mr,1)*100:.1f}%"
          f" of multi-RCV)")
    print(f"  Diff review status across:      "
          f"{results['multi_diff_review_status']:>12,}  "
          f"({results['multi_diff_review_status']/max(mr,1)*100:.1f}%"
          f" of multi-RCV)")

    print(f"\n{'MASKING IMPACT (key metrics)':}")
    print(f"  Variants where higher-weight")
    print(f"    status was selected:          "
          f"{results['masking_review_upgrades']:>12,}  "
          f"({results['masking_review_upgrades_pct']:.3f}% of all variants)")
    print(f"  Variants where 'conflicting'")
    print(f"    was hidden by non-conflicting: "
          f"{results['masking_hides_conflict']:>12,}  "
          f"({results['masking_hides_conflict_pct']:.3f}% of all variants)")
    print(f"  Variants where selected class")
    print(f"    differs from masked class:    "
          f"{results['masking_changes_classification']:>12,}  "
          f"({results['masking_changes_classification_pct']:.3f}%"
          f" of all variants)")
    print(f"  Concordant condition hides")
    print(f"    discordant condition:          "
          f"{results['masking_concordant_hides_discordant']:>12,}  "
          f"({results['masking_concordant_hides_discordant_pct']:.3f}%"
          f" of all variants)")

    if results['classification_impact']:
        print(f"\n{'TOP CLASSIFICATION MISMATCHES':}")
        print(f"  {'Selected':<15} {'Masked':<15} {'Count':>8}")
        print(f"  {'-'*15} {'-'*15} {'-'*8}")
        for transition, count in list(
                results['classification_impact'].items())[:15]:
            parts = transition.split('->')
            print(f"  {parts[0]:<15} {parts[1]:<15} {count:>8,}")

    if results['masked_examples']:
        print(f"\n{'EXAMPLE MASKED VARIANTS (first 5)':}")
        for ex in results['masked_examples'][:5]:
            print(f"  Variant {ex['variation_id']} ({ex.get('gene', '?')}):")
            print(f"    SELECTED: {ex['selected_condition']}")
            print(f"      Class: {ex['selected_classification']} "
                  f"| Status: {ex['selected_review_status']} "
                  f"(weight={ex['selected_weight']})")
            print(f"    MASKED:   {ex['masked_condition']}")
            print(f"      Class: {ex['masked_classification']} "
                  f"| Status: {ex['masked_review_status']} "
                  f"(weight={ex['masked_weight']})")
            print()

    print(f"\n{'INTERPRETATION':}")
    mc = results['masking_changes_classification']
    hc = results['masking_hides_conflict']
    print(f"  Out of {tv:,} variants:")
    print(f"  - {mc:,} ({results['masking_changes_classification_pct']:.3f}%)"
          f" may have been assigned a different classification")
    print(f"    due to highest-weight selection (classification-level impact)")
    print(f"  - {hc:,} ({results['masking_hides_conflict_pct']:.3f}%)"
          f" had a condition-specific conflict hidden")
    print(f"    (concordance-level impact)")
    print(f"  This represents an upper bound — reclassification detection")
    print(f"  tracks classification CHANGES between snapshots, so a")
    print(f"  consistently selected condition still produces correct")
    print(f"  change detection.")
    print()
    print(f"  Parse time: {results['parse_time_seconds']:.1f}s | "
          f"Analysis time: {results['analysis_time_seconds']:.1f}s")
    print("=" * 70)


# ============================================================================
# MULTI-SNAPSHOT ANALYSIS
# ============================================================================

def find_rcv_snapshots(snapshot_dir: Path) -> List[Path]:
    """Find all RCV-format XML snapshots in the directory."""
    snapshots = []
    for f in sorted(snapshot_dir.iterdir()):
        if (f.name.endswith('.xml.gz') and
                'ClinVarFullRelease' in f.name and
                f.stat().st_size > 1000):
            snapshots.append(f)
    return snapshots


def analyze_multiple_snapshots(snapshot_dir: Path, output_dir: Path,
                                all_rcv: bool = False,
                                sample_years: List[int] = None):
    """Analyze multiple RCV snapshots and write combined results.

    By default, analyzes a representative sample:
      - 2013-12 (early era)
      - 2016-12 (mid era)
      - 2019-12 (last RCV snapshot, most relevant)

    With --all-rcv, analyzes every RCV snapshot.
    """
    snapshots = find_rcv_snapshots(snapshot_dir)
    if not snapshots:
        print(f"ERROR: No RCV snapshots found in {snapshot_dir}")
        sys.exit(1)

    print(f"Found {len(snapshots)} RCV snapshots")

    if not all_rcv:
        # Select representative snapshots
        if sample_years is None:
            sample_years = [2013, 2016, 2019]

        selected = []
        for year in sample_years:
            # Find the December snapshot for each year, or the latest available
            year_snaps = [s for s in snapshots if f'{year}' in s.name]
            dec_snaps = [s for s in year_snaps if f'{year}-12' in s.name]
            if dec_snaps:
                selected.append(dec_snaps[0])
            elif year_snaps:
                selected.append(year_snaps[-1])

        if not selected:
            print("ERROR: Could not find any snapshots for the requested years")
            sys.exit(1)

        snapshots = selected
        print(f"Selected {len(snapshots)} representative snapshots:")
        for s in snapshots:
            print(f"  {s.name}")

    # Run analysis
    all_results = []
    for snapshot_path in snapshots:
        result = analyze_snapshot(snapshot_path)
        all_results.append(result)

    # Write combined results
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "rcv_masking_analysis.json"

    combined = {
        'generated_at': datetime.now().isoformat(),
        'description': (
            "Analysis of cross-condition masking in RCV-era parsing. "
            "Quantifies how many variants have condition-specific conflicts "
            "or classification differences that are hidden by selecting "
            "the highest-weight review status record."
        ),
        'num_snapshots_analyzed': len(all_results),
        'snapshots': all_results,
        'summary': _build_trend_summary(all_results),
    }

    with open(output_path, 'w') as f:
        json.dump(combined, f, indent=2)
    print(f"\nResults written to: {output_path}")

    # Also write a CSV for easy plotting
    csv_path = output_dir / "rcv_masking_trend.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'snapshot', 'total_variants', 'multi_rcv_variants',
            'multi_rcv_pct', 'masking_review_upgrades',
            'masking_hides_conflict', 'masking_changes_classification',
            'masking_changes_classification_pct',
            'masking_hides_conflict_pct',
        ])
        for r in all_results:
            writer.writerow([
                r['snapshot'], r['total_variants'], r['multi_rcv_variants'],
                r['multi_rcv_pct'], r['masking_review_upgrades'],
                r['masking_hides_conflict'],
                r['masking_changes_classification'],
                r['masking_changes_classification_pct'],
                r['masking_hides_conflict_pct'],
            ])
    print(f"Trend CSV written to: {csv_path}")

    return combined


def _build_trend_summary(results: List[Dict]) -> Dict:
    """Build a trend summary across multiple snapshots."""
    if not results:
        return {}

    latest = results[-1]
    return {
        'latest_snapshot': latest['snapshot'],
        'latest_total_variants': latest['total_variants'],
        'latest_multi_rcv_variants': latest['multi_rcv_variants'],
        'latest_multi_rcv_pct': latest['multi_rcv_pct'],
        'latest_masking_changes_classification':
            latest['masking_changes_classification'],
        'latest_masking_changes_classification_pct':
            latest['masking_changes_classification_pct'],
        'latest_masking_hides_conflict':
            latest['masking_hides_conflict'],
        'latest_masking_hides_conflict_pct':
            latest['masking_hides_conflict_pct'],
        'conclusion': (
            f"In the latest RCV snapshot ({latest['snapshot']}), "
            f"{latest['masking_changes_classification']:,} variants "
            f"({latest['masking_changes_classification_pct']:.3f}% of "
            f"{latest['total_variants']:,} total) may have been assigned "
            f"a different classification due to highest-weight selection. "
            f"{latest['masking_hides_conflict']:,} variants "
            f"({latest['masking_hides_conflict_pct']:.3f}%) had a "
            f"condition-specific conflict hidden. These represent upper "
            f"bounds as reclassification detection tracks changes between "
            f"snapshots rather than absolute classifications."
        ),
    }


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="RCV Masking Analysis — Quantify impact of "
                    "highest-weight review status selection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze the last RCV snapshot only
  python rcv_masking_analysis.py \\
      --snapshot /path/to/ClinVarFullRelease_2019-12.xml.gz \\
      --output ./masking_results/

  # Analyze representative snapshots (2013, 2016, 2019)
  python rcv_masking_analysis.py \\
      --snapshot-dir /path/to/ftp_snapshots/ \\
      --output ./masking_results/

  # Analyze all RCV snapshots
  python rcv_masking_analysis.py \\
      --snapshot-dir /path/to/ftp_snapshots/ \\
      --all-rcv \\
      --output ./masking_results/
        """
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--snapshot', '-s', type=Path,
        help="Path to a single RCV XML snapshot (.xml.gz)")
    group.add_argument(
        '--snapshot-dir', '-d', type=Path,
        help="Directory containing RCV snapshots (ftp_snapshots/)")

    parser.add_argument(
        '--output', '-o', type=Path, default=Path('./masking_results'),
        help="Output directory (default: ./masking_results/)")
    parser.add_argument(
        '--all-rcv', action='store_true',
        help="Analyze all RCV snapshots (default: sample 2013, 2016, 2019)")
    parser.add_argument(
        '--years', nargs='+', type=int, default=None,
        help="Specific years to sample (e.g., --years 2015 2017 2019)")

    args = parser.parse_args()

    print("=" * 70)
    print("RCV MASKING ANALYSIS")
    print("Quantifying cross-condition masking from highest-weight selection")
    print("=" * 70)

    if args.snapshot:
        # Single snapshot mode
        if not args.snapshot.exists():
            print(f"ERROR: Snapshot not found: {args.snapshot}")
            sys.exit(1)

        result = analyze_snapshot(args.snapshot)

        args.output.mkdir(parents=True, exist_ok=True)
        output_path = args.output / "rcv_masking_analysis.json"
        with open(output_path, 'w') as f:
            json.dump({
                'generated_at': datetime.now().isoformat(),
                'description': (
                    "Analysis of cross-condition masking in RCV-era parsing."
                ),
                'num_snapshots_analyzed': 1,
                'snapshots': [result],
                'summary': _build_trend_summary([result]),
            }, f, indent=2)
        print(f"\nResults written to: {output_path}")

    else:
        # Multi-snapshot mode
        if not args.snapshot_dir.exists():
            print(f"ERROR: Snapshot directory not found: {args.snapshot_dir}")
            sys.exit(1)

        analyze_multiple_snapshots(
            args.snapshot_dir, args.output,
            all_rcv=args.all_rcv,
            sample_years=args.years)


if __name__ == '__main__':
    main()
