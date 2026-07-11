#!/usr/bin/env python3
"""
ClinVar Submitter Landscape Analysis — Server-Side Script
===========================================================
Designed to run on the Columbia DBMI server where checkpoint pickle files
and FTP snapshot XML files are available.

IMPORTANT: VCV-era checkpoints (2020+) have empty submitter fields due to
a parsing bug (ClinVarSubmissionID doesn't carry the 'submitter' attribute
in VCV XML — the name lives on ClinVarAccession/@SubmitterName instead).
This script automatically detects this and re-extracts submitter names
from the original XML files for affected snapshots.

Analyses produced:
  A. OVERALL RANKINGS       — top N submitters by total submissions
  B. SNAPSHOT-BASED TEMPORAL — month-by-month submitter counts (148 snapshots)
  C. SUBMITTER x GENE       — which submitters dominate which genes
  D. RECLASSIFICATION LINK  — which submitter triggered each reclassification
  E. SUBMITTER ENTRY/EXIT   — when each submitter first/last appeared
  F. CONFLICT PAIRS         — which submitter pairs disagree most often

Usage (on the server):
    nohup python3 clinvar_submitter_server.py \
        --checkpoint-dir clinvar_output/checkpoints \
        --snapshot-dir   clinvar_cache/ftp_snapshots \
        --output-dir     submitter_analysis > submitter_analysis.log 2>&1 &
    tail -f submitter_analysis.log
"""

import argparse
import gc
import gzip
import json
import pickle
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from xml.etree.ElementTree import iterparse


# ── Classification category mapping (must match clinvar_ftp_only_v5) ──

_CATEGORY_MAP = {
    'pathogenic':            'P',
    'likely pathogenic':     'LP',
    'uncertain significance':'VUS',
    'likely benign':         'LB',
    'benign':                'B',
    'conflicting interpretations of pathogenicity': 'Conflicting',
    'conflicting classifications of pathogenicity': 'Conflicting',
    'drug response':         'Other',
    'risk factor':           'Other',
    'association':           'Other',
    'protective':            'Other',
    'affects':               'Other',
    'not provided':          'Other',
    'other':                 'Other',
}

def categorize(cls_str: str) -> str:
    if not cls_str:
        return 'Other'
    key = cls_str.strip().lower()
    if key in _CATEGORY_MAP:
        return _CATEGORY_MAP[key]
    if 'pathogenic' in key and 'likely' in key:
        return 'LP'
    if 'pathogenic' in key:
        return 'P'
    if 'benign' in key and 'likely' in key:
        return 'LB'
    if 'benign' in key:
        return 'B'
    if 'uncertain' in key or 'vus' in key:
        return 'VUS'
    if 'conflict' in key:
        return 'Conflicting'
    return 'Other'


# ═══════════════════════════════════════════════════════════════════════
# CHECKPOINT LOADING
# ═══════════════════════════════════════════════════════════════════════

def find_checkpoints(checkpoint_dir: Path) -> List[Tuple[Path, str]]:
    """Find and sort all checkpoint files by date."""
    checkpoints = []
    for f in checkpoint_dir.glob("*.checkpoint.pkl.gz"):
        if f.stat().st_size < 100:
            continue
        match = re.search(r'(\d{4}-\d{2})', f.name)
        if match:
            checkpoints.append((f, match.group(1)))
    checkpoints.sort(key=lambda x: x[1])
    return checkpoints


def load_checkpoint(path: Path) -> Optional[Tuple[str, Dict]]:
    """Load a single checkpoint file."""
    try:
        with gzip.open(path, 'rb') as f:
            ckpt = pickle.load(f)
        return ckpt['snapshot_date'], ckpt['data']
    except Exception as e:
        print(f"  WARNING: Failed to load {path.name}: {e}", flush=True)
        return None


# ═══════════════════════════════════════════════════════════════════════
# XML SUBMITTER EXTRACTION (for VCV snapshots with missing data)
# ═══════════════════════════════════════════════════════════════════════

# Pre-compiled regexes for fast line-by-line extraction
_RE_VAR_ID = re.compile(r'<VariationArchive[^>]+VariationID="(\d+)"')
_RE_SUBMITTER = re.compile(r'SubmitterName="([^"]*)"')
_RE_ORG_ID = re.compile(r'OrgID="([^"]*)"')
_RE_SCV_TYPE = re.compile(r'Type="SCV"')

def extract_submitters_from_xml(xml_path: Path) -> Dict[str, List[str]]:
    """
    Fast regex-based scan of a VCV XML file.
    Avoids full XML parsing — just scans lines for VariationID and
    SubmitterName attributes on SCV-type ClinVarAccession elements.

    ~10-20x faster than iterparse for this specific extraction.

    Returns: {variation_id: [submitter_name, ...]}
    """
    result = {}
    opener = gzip.open if str(xml_path).endswith('.gz') else open

    current_var_id = None
    current_subs = []

    with opener(xml_path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            # Detect new VariationArchive
            m_var = _RE_VAR_ID.search(line)
            if m_var:
                # Save previous variant
                if current_var_id and current_subs:
                    result[current_var_id] = current_subs
                current_var_id = m_var.group(1)
                current_subs = []
                continue

            # Detect SCV ClinVarAccession with SubmitterName
            if 'SubmitterName=' in line and _RE_SCV_TYPE.search(line):
                m_sub = _RE_SUBMITTER.search(line)
                if m_sub:
                    name = m_sub.group(1)
                    if name:
                        current_subs.append(name)

    # Save last variant
    if current_var_id and current_subs:
        result[current_var_id] = current_subs

    return result


def find_xml_for_checkpoint(ckpt_path: Path, snapshot_dir: Path) -> Optional[Path]:
    """Find the FTP XML file corresponding to a checkpoint."""
    # Extract date from checkpoint name
    match = re.search(r'(\d{4}-\d{2})', ckpt_path.name)
    if not match:
        return None
    date_str = match.group(1)

    # Try VCV naming patterns
    patterns = [
        f"ClinVarVariationRelease_{date_str}*.xml.gz",
        f"ClinVarVCVRelease_{date_str}*.xml.gz",
        f"ClinVarFullRelease_{date_str}*.xml.gz",
    ]
    for pat in patterns:
        matches = list(snapshot_dir.glob(pat))
        if matches:
            return matches[0]
    return None


def patch_checkpoint_submitters(data: Dict, xml_submitters: Dict) -> int:
    """Patch checkpoint data with submitter names from XML extraction.
    xml_submitters: {var_id: [submitter_name, ...]}
    Returns count of variants patched."""
    patched = 0
    for var_id, record in data.items():
        names = xml_submitters.get(var_id)
        if not names:
            continue

        # Patch all_submitters
        record['all_submitters'] = names
        patched += 1

        # Patch individual submission records if counts match
        existing_subs = record.get('submissions', [])
        if existing_subs and len(existing_subs) == len(names):
            for i, name in enumerate(names):
                existing_subs[i]['submitter'] = name

    return patched


# ═══════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def analyse_submitter_landscape(checkpoint_dir: Path, snapshot_dir: Optional[Path],
                                output_dir: Path):
    checkpoints = find_checkpoints(checkpoint_dir)
    if not checkpoints:
        print("ERROR: No checkpoint files found!")
        sys.exit(1)

    print(f"Found {len(checkpoints)} checkpoints ({checkpoints[0][1]} – {checkpoints[-1][1]})",
          flush=True)

    # Accumulators
    temporal_counts = {}
    snapshot_variant_counts = {}
    submitter_lifespan = {}
    snapshot_active_submitters = {}

    latest_date = None
    latest_submitter_total = Counter()
    latest_submitter_genes = defaultdict(Counter)
    latest_submitter_clinsig = defaultdict(Counter)
    latest_submitter_review = defaultdict(Counter)

    prev_snapshot_submitters = {}
    prev_snapshot_date = None
    prev_snapshot_categories = {}

    reclass_by_new_submitter = Counter()
    reclass_by_existing_submitter = Counter()
    conflict_pair_counts = Counter()

    n_checkpoints_loaded = 0
    n_xml_patched = 0
    total_reclass_events = 0

    for ckpt_path, ckpt_date in checkpoints:
        t0 = time.time()
        result = load_checkpoint(ckpt_path)
        if result is None:
            continue
        snapshot_date, data = result
        n_checkpoints_loaded += 1
        snapshot_variant_counts[snapshot_date] = len(data)

        # ── Detect missing submitter data and patch from XML ──────
        # Sample a few records to check if submitters are populated
        sample_keys = list(data.keys())[:20]
        has_submitters = any(
            bool(data[k].get('all_submitters', []))
            for k in sample_keys
        )

        if not has_submitters and snapshot_dir and snapshot_dir.exists():
            xml_path = find_xml_for_checkpoint(ckpt_path, snapshot_dir)
            if xml_path and xml_path.exists():
                print(f"  [{snapshot_date}] Missing submitter data — re-extracting from {xml_path.name} ...",
                      flush=True)
                xml_subs = extract_submitters_from_xml(xml_path)
                patched = patch_checkpoint_submitters(data, xml_subs)
                print(f"  [{snapshot_date}] Patched {patched:,} / {len(data):,} variants with submitter data",
                      flush=True)
                n_xml_patched += 1
                del xml_subs
            elif not xml_path:
                print(f"  [{snapshot_date}] WARNING: No XML file found for patching", flush=True)

        # ── Per-snapshot submitter counting ────────────────────────
        submitter_count_this_snapshot = Counter()
        active_submitters = set()

        for var_id, record in data.items():
            submitters = record.get('all_submitters', [])
            gene = record.get('gene_symbol', '')
            classification = record.get('classification', '')
            category = categorize(classification)

            for sub in submitters:
                if not sub:
                    continue
                submitter_count_this_snapshot[sub] += 1
                active_submitters.add(sub)
                if sub not in submitter_lifespan:
                    submitter_lifespan[sub] = [snapshot_date, snapshot_date]
                else:
                    submitter_lifespan[sub][1] = snapshot_date

            # ── Reclassification → submitter linkage ──────────────
            if prev_snapshot_date is not None and var_id in prev_snapshot_submitters:
                prev_cat = prev_snapshot_categories.get(var_id, '')
                curr_cat = category
                if (prev_cat and curr_cat and prev_cat != curr_cat
                        and curr_cat != 'Other' and prev_cat != 'Other'):
                    total_reclass_events += 1
                    prev_subs = prev_snapshot_submitters.get(var_id, set())
                    curr_subs = set(s for s in submitters if s)
                    new_subs = curr_subs - prev_subs
                    if new_subs:
                        for ns in new_subs:
                            reclass_by_new_submitter[ns] += 1
                    else:
                        for es in curr_subs:
                            reclass_by_existing_submitter[es] += 1

            # ── Conflict pair detection ───────────────────────────
            if len(submitters) >= 2:
                subs_list = record.get('submissions', [])
                if len(subs_list) >= 2:
                    sub_cls = {}
                    for s in subs_list:
                        sname = s.get('submitter', '')
                        scls = categorize(s.get('classification', ''))
                        if sname and scls and scls != 'Other':
                            sub_cls[sname] = scls
                    snames = list(sub_cls.keys())
                    for i in range(len(snames)):
                        for j in range(i + 1, len(snames)):
                            if sub_cls[snames[i]] != sub_cls[snames[j]]:
                                pair = tuple(sorted([snames[i], snames[j]]))
                                conflict_pair_counts[pair] += 1

        temporal_counts[snapshot_date] = dict(submitter_count_this_snapshot.most_common())
        snapshot_active_submitters[snapshot_date] = len(active_submitters)

        # Update prev_snapshot for reclassification linkage
        prev_snapshot_submitters = {}
        prev_snapshot_categories = {}
        for var_id, record in data.items():
            subs = set(s for s in record.get('all_submitters', []) if s)
            if subs:
                prev_snapshot_submitters[var_id] = subs
            cls = record.get('classification', '')
            if cls:
                prev_snapshot_categories[var_id] = categorize(cls)
        prev_snapshot_date = snapshot_date

        # Capture detailed stats for the latest snapshot
        latest_date = snapshot_date
        latest_submitter_total = Counter()
        latest_submitter_genes = defaultdict(Counter)
        latest_submitter_clinsig = defaultdict(Counter)
        latest_submitter_review = defaultdict(Counter)
        for var_id, record in data.items():
            gene = record.get('gene_symbol', '')
            subs_list = record.get('submissions', [])
            if subs_list:
                for s in subs_list:
                    sname = s.get('submitter', '')
                    if not sname:
                        continue
                    latest_submitter_total[sname] += 1
                    if gene:
                        latest_submitter_genes[sname][gene] += 1
                    scls = s.get('classification', '')
                    if scls:
                        latest_submitter_clinsig[sname][scls] += 1
                    srev = s.get('review_status', '')
                    if srev:
                        latest_submitter_review[sname][srev] += 1
            else:
                for sub in record.get('all_submitters', []):
                    if not sub:
                        continue
                    latest_submitter_total[sub] += 1
                    if gene:
                        latest_submitter_genes[sub][gene] += 1

        elapsed = time.time() - t0
        print(f"  [{snapshot_date}] {len(data):,} variants, "
              f"{len(active_submitters)} submitters ({elapsed:.1f}s)", flush=True)

        del data
        gc.collect()

    print(f"\nLoaded {n_checkpoints_loaded} checkpoints "
          f"({n_xml_patched} patched from XML), "
          f"tracked {total_reclass_events:,} reclassification events", flush=True)

    # ═══════════════════════════════════════════════════════════════════
    # BUILD OUTPUT JSON FILES
    # ═══════════════════════════════════════════════════════════════════

    output_dir.mkdir(parents=True, exist_ok=True)

    # ── A. Overall rankings ───────────────────────────────────────────
    top_n = 50
    overall_list = []
    total_all = sum(latest_submitter_total.values())
    for rank, (name, count) in enumerate(latest_submitter_total.most_common(top_n), 1):
        sig_counts = dict(latest_submitter_clinsig.get(name, Counter()).most_common())
        rev_counts = dict(latest_submitter_review.get(name, Counter()).most_common())
        top_genes = [{"gene": g, "count": c}
                     for g, c in latest_submitter_genes.get(name, Counter()).most_common(10)]
        overall_list.append({
            "rank": rank,
            "submitter": name,
            "total_submissions": count,
            "pct_of_total": round(100 * count / total_all, 2) if total_all else 0,
            "top_10_genes": top_genes,
            "classification_breakdown": sig_counts,
            "review_status_breakdown": rev_counts,
        })

    with open(output_dir / "submitter_overall.json", "w") as f:
        json.dump(overall_list, f, indent=2)
    print(f"Wrote submitter_overall.json ({len(overall_list)} entries)", flush=True)

    # ── B. Temporal ───────────────────────────────────────────────────
    top20_names = [e["submitter"] for e in overall_list[:20]]
    all_dates = sorted(temporal_counts.keys())

    temporal_out = {"dates": all_dates, "submitters": {}}
    for name in top20_names:
        temporal_out["submitters"][name] = {
            d: temporal_counts[d].get(name, 0) for d in all_dates
        }
    other_by_date = {}
    for d in all_dates:
        top20_sum = sum(temporal_counts[d].get(n, 0) for n in top20_names)
        total_d = sum(temporal_counts[d].values())
        other_by_date[d] = total_d - top20_sum
    temporal_out["submitters"]["_Other"] = other_by_date
    temporal_out["total_per_snapshot"] = {
        d: sum(temporal_counts[d].values()) for d in all_dates
    }
    temporal_out["active_submitters_per_snapshot"] = {
        d: snapshot_active_submitters.get(d, 0) for d in all_dates
    }

    with open(output_dir / "submitter_temporal.json", "w") as f:
        json.dump(temporal_out, f, indent=2)
    print(f"Wrote submitter_temporal.json ({len(all_dates)} snapshots)", flush=True)

    # ── Aggregate stats ───────────────────────────────────────────────
    total_subs_latest = sum(latest_submitter_total.values())
    top10_sum = sum(c for _, c in latest_submitter_total.most_common(10))
    top20_sum_val = sum(c for _, c in latest_submitter_total.most_common(20))

    stats = {
        "latest_snapshot": latest_date,
        "total_submissions_latest": total_subs_latest,
        "total_unique_submitters_latest": len(latest_submitter_total),
        "total_unique_submitters_ever": len(submitter_lifespan),
        "top10_concentration_pct": round(100 * top10_sum / total_subs_latest, 2) if total_subs_latest else 0,
        "top20_concentration_pct": round(100 * top20_sum_val / total_subs_latest, 2) if total_subs_latest else 0,
        "snapshots_processed": n_checkpoints_loaded,
        "snapshots_xml_patched": n_xml_patched,
        "total_reclassification_events_tracked": total_reclass_events,
    }

    with open(output_dir / "submitter_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Wrote submitter_stats.json", flush=True)

    # ── C. Submitter x gene matrix ────────────────────────────────────
    top30_genes = Counter()
    for name in top20_names:
        top30_genes.update(latest_submitter_genes.get(name, {}))
    top30_gene_names = [g for g, _ in top30_genes.most_common(30)]

    gene_matrix = {"submitters": top20_names, "genes": top30_gene_names, "matrix": {}}
    for name in top20_names:
        gene_counts = latest_submitter_genes.get(name, {})
        gene_matrix["matrix"][name] = {g: gene_counts.get(g, 0) for g in top30_gene_names}

    with open(output_dir / "submitter_gene_matrix.json", "w") as f:
        json.dump(gene_matrix, f, indent=2)
    print(f"Wrote submitter_gene_matrix.json", flush=True)

    # ── D. Reclassification linkage ───────────────────────────────────
    reclass_data = {
        "total_events_tracked": total_reclass_events,
        "events_with_new_submitter": sum(reclass_by_new_submitter.values()),
        "events_with_reinterpretation": sum(reclass_by_existing_submitter.values()),
        "top_30_new_submitter_triggers": [
            {"submitter": s, "events_triggered": c}
            for s, c in reclass_by_new_submitter.most_common(30)
        ],
        "top_30_reinterpretation": [
            {"submitter": s, "events_involved": c}
            for s, c in reclass_by_existing_submitter.most_common(30)
        ],
    }

    with open(output_dir / "submitter_reclassification.json", "w") as f:
        json.dump(reclass_data, f, indent=2)
    print(f"Wrote submitter_reclassification.json "
          f"({reclass_data['events_with_new_submitter']:,} new-submitter, "
          f"{reclass_data['events_with_reinterpretation']:,} reinterpretation)", flush=True)

    # ── E. Entry / exit ───────────────────────────────────────────────
    entry_exit = {
        "submitters": [
            {
                "submitter": name,
                "first_seen": dates[0],
                "last_seen": dates[1],
                "active_months": sum(
                    1 for d in all_dates
                    if dates[0] <= d <= dates[1]
                    and temporal_counts.get(d, {}).get(name, 0) > 0
                ),
            }
            for name, dates in sorted(submitter_lifespan.items(),
                                      key=lambda x: x[1][0])
        ]
    }

    with open(output_dir / "submitter_entry_exit.json", "w") as f:
        json.dump(entry_exit, f, indent=2)
    print(f"Wrote submitter_entry_exit.json ({len(entry_exit['submitters'])} submitters)", flush=True)

    # ── F. Conflict pairs ─────────────────────────────────────────────
    conflict_data = {
        "top_50_conflict_pairs": [
            {"submitter_a": pair[0], "submitter_b": pair[1], "variant_count": cnt}
            for pair, cnt in conflict_pair_counts.most_common(50)
        ],
        "total_conflict_variant_instances": sum(conflict_pair_counts.values()),
        "unique_conflict_pairs": len(conflict_pair_counts),
    }

    with open(output_dir / "submitter_conflict_pairs.json", "w") as f:
        json.dump(conflict_data, f, indent=2)
    print(f"Wrote submitter_conflict_pairs.json ({len(conflict_pair_counts)} unique pairs)", flush=True)

    print("\n" + "=" * 70)
    print("Submitter landscape analysis complete!")
    print(f"  Snapshots: {n_checkpoints_loaded} ({n_xml_patched} patched from XML)")
    print(f"  Unique submitters (ever): {len(submitter_lifespan)}")
    print(f"  Unique submitters (latest): {len(latest_submitter_total)}")
    print(f"  Reclassification events linked: {total_reclass_events:,}")
    print(f"  Conflict pairs found: {len(conflict_pair_counts)}")
    print(f"  Output directory: {output_dir}")
    print("=" * 70, flush=True)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="ClinVar Submitter Landscape Analysis (Server-Side)")
    parser.add_argument("--checkpoint-dir", required=True,
                        help="Directory containing *.checkpoint.pkl.gz files")
    parser.add_argument("--snapshot-dir", default=None,
                        help="Directory containing FTP XML snapshots "
                             "(e.g., clinvar_cache/ftp_snapshots). "
                             "Required for VCV-era submitter extraction.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for output JSON files")
    args = parser.parse_args()

    checkpoint_dir = Path(args.checkpoint_dir).expanduser()
    snapshot_dir = Path(args.snapshot_dir).expanduser() if args.snapshot_dir else None
    output_dir = Path(args.output_dir).expanduser()

    if not checkpoint_dir.exists():
        print(f"ERROR: Checkpoint directory not found: {checkpoint_dir}")
        sys.exit(1)

    print("=" * 70)
    print("ClinVar Submitter Landscape Analysis (Server-Side)")
    print("=" * 70)
    print(f"Checkpoint directory: {checkpoint_dir}")
    print(f"Snapshot directory:   {snapshot_dir or '(not provided — VCV submitters will be empty)'}")
    print(f"Output directory:     {output_dir}")
    print(flush=True)

    start = time.time()
    analyse_submitter_landscape(checkpoint_dir, snapshot_dir, output_dir)
    elapsed = time.time() - start
    print(f"\nTotal runtime: {elapsed / 60:.1f} minutes", flush=True)


if __name__ == "__main__":
    main()
