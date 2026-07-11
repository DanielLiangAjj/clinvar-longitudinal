#!/usr/bin/env python3
"""
ClinVar XML Structure Analyzer — Run on server with XML snapshots
=================================================================
Analyzes structural evolution of ClinVar XML releases for Q1, Q2, Q6, Q7.
Produces xml_structure_analysis.json for the visualization script.

Usage (on server):
    python3 clinvar_xml_structure_analyzer.py \
        --snapshot-dir ./clinvar_cache/ftp_snapshots \
        --output xml_structure_analysis.json \
        --max-records 3000

Then copy result to local machine:
    scp ly2680@aqua.dbmi.columbia.edu:/public/variant_reclassification/xml_structure_analysis.json ~/Downloads/
"""

import xml.etree.ElementTree as ET
import gzip
import json
import re
import sys
import time
import argparse
from pathlib import Path
from collections import defaultdict, Counter


def analyze_snapshot(xml_path, max_records=3000):
    """Analyze XML structure of a single ClinVar snapshot by sampling records."""
    elements = set()
    attributes_per_element = defaultdict(set)
    element_paths = set()          # full ancestor paths e.g. "VariationArchive/ClassifiedRecord/..."
    max_depth = 0
    depth_sum = 0
    depth_count = 0

    # Evidence & linkage counters
    xref_databases = Counter()
    citation_sources = Counter()
    records_with_xref = 0
    records_with_citation = 0
    classification_strings = Counter()

    record_count = 0
    record_tag = None
    current_path = []              # stack for building full paths
    has_xref_this_record = False
    has_citation_this_record = False

    print(f"  Analyzing: {xml_path.name} ...", end=" ", flush=True)
    t0 = time.time()

    try:
        with gzip.open(xml_path, "rb") as f:
            context = ET.iterparse(f, events=("start", "end"))
            root = None

            for event, elem in context:
                if event == "start":
                    if root is None:
                        root = elem

                    tag = elem.tag
                    current_path.append(tag)
                    depth = len(current_path)

                    elements.add(tag)
                    for attr_name in elem.attrib:
                        attributes_per_element[tag].add(attr_name)

                    # Track full XPath-like path (limit depth to avoid explosion)
                    if depth <= 8:
                        element_paths.add("/".join(current_path[:8]))

                    if depth > max_depth:
                        max_depth = depth
                    depth_sum += depth
                    depth_count += 1

                    # Detect record boundary
                    if record_tag is None and tag in ("ClinVarSet", "VariationArchive"):
                        record_tag = tag

                    # Track XRef databases
                    if tag == "XRef":
                        db = elem.get("DB", "")
                        if db:
                            xref_databases[db] += 1
                            has_xref_this_record = True

                    # Track citations
                    if tag == "ID":
                        src = elem.get("Source", "")
                        if src:
                            citation_sources[src] += 1
                            has_citation_this_record = True

                    # Track classification strings
                    if tag == "Description" and len(current_path) >= 2:
                        parent = current_path[-2] if len(current_path) >= 2 else ""
                        if parent in (
                            "GermlineClassification",
                            "SomaticClinicalImpact",
                            "OncogenicityClassification",
                            "ClinicalSignificance",
                            "Classification",
                        ):
                            if elem.text and elem.text.strip():
                                classification_strings[elem.text.strip()] += 1

                elif event == "end":
                    if current_path:
                        current_path.pop()

                    if elem.tag == record_tag:
                        record_count += 1
                        if has_xref_this_record:
                            records_with_xref += 1
                        if has_citation_this_record:
                            records_with_citation += 1
                        has_xref_this_record = False
                        has_citation_this_record = False

                        if root is not None:
                            root.clear()
                        if record_count >= max_records:
                            break

    except Exception as e:
        print(f"ERROR: {e}")
        return None

    elapsed = time.time() - t0
    print(f"{record_count:,} records in {elapsed:.1f}s")

    # Extract date from filename
    match = re.search(r"(\d{4}-\d{2})", xml_path.name)
    snapshot_date = match.group(1) if match else xml_path.stem

    fmt = "rcv" if record_tag == "ClinVarSet" else "vcv"

    return {
        "snapshot_name": xml_path.name,
        "snapshot_date": snapshot_date,
        "format": fmt,
        "records_sampled": record_count,
        "num_unique_elements": len(elements),
        "elements": sorted(elements),
        "total_unique_attributes": sum(len(v) for v in attributes_per_element.values()),
        "attributes_per_element": {k: sorted(v) for k, v in attributes_per_element.items()},
        "num_element_paths": len(element_paths),
        "max_depth": max_depth,
        "avg_depth": round(depth_sum / max(depth_count, 1), 2),
        "xref_databases": dict(xref_databases.most_common()),
        "citation_sources": dict(citation_sources.most_common()),
        "records_with_xref": records_with_xref,
        "records_with_citation": records_with_citation,
        "pct_with_xref": round(records_with_xref / max(record_count, 1) * 100, 1),
        "pct_with_citation": round(records_with_citation / max(record_count, 1) * 100, 1),
        "classification_strings": dict(classification_strings.most_common(50)),
        "num_distinct_classifications": len(classification_strings),
    }


def compute_lifecycle(snapshots_data):
    """Compute element lifecycle (new/retained/deprecated) across releases."""
    lifecycle = []
    prev_elements = None

    for snap in snapshots_data:
        current_elements = set(snap["elements"])
        entry = {
            "snapshot_date": snap["snapshot_date"],
            "total_elements": len(current_elements),
        }
        if prev_elements is not None:
            new = current_elements - prev_elements
            deprecated = prev_elements - current_elements
            retained = current_elements & prev_elements
            entry["new_elements"] = len(new)
            entry["deprecated_elements"] = len(deprecated)
            entry["retained_elements"] = len(retained)
            entry["new_element_names"] = sorted(new)
            entry["deprecated_element_names"] = sorted(deprecated)
            entry["pct_new"] = round(len(new) / max(len(current_elements), 1) * 100, 1)
            entry["pct_deprecated"] = round(
                len(deprecated) / max(len(prev_elements), 1) * 100, 1
            )
            entry["retention_rate"] = round(
                len(retained) / max(len(prev_elements), 1) * 100, 1
            )
        else:
            entry["new_elements"] = len(current_elements)
            entry["deprecated_elements"] = 0
            entry["retained_elements"] = 0
            entry["pct_new"] = 100.0
            entry["pct_deprecated"] = 0.0
            entry["retention_rate"] = 0.0

        lifecycle.append(entry)
        prev_elements = current_elements

    return lifecycle


def compute_field_availability(snapshots_data):
    """Track which key data-carrying elements are present in each release."""
    key_fields = [
        "VariationArchive", "ClinVarSet", "SimpleAllele", "Gene",
        "GermlineClassification", "ClinicalSignificance", "Classification",
        "SomaticClinicalImpact", "OncogenicityClassification",
        "ReviewStatus", "Description", "HGVSlist", "HGVS",
        "NucleotideExpression", "ProteinExpression", "NucleotideChange", "ProteinChange",
        "SequenceLocation", "XRef", "Citation", "MeasureSet", "Measure",
        "MeasureRelationship", "ClinVarSubmissionID",
        "ClinicalAssertion", "ClinVarAssertion",
        "ObservedIn", "TraitSet", "Trait",
        "AttributeSet", "Comment", "FunctionalConsequence",
        "Haplotype", "Genotype", "GenotypeSet",
        "Classifications", "ClassifiedRecord",
    ]

    availability = {}
    for snap in snapshots_data:
        elem_set = set(snap["elements"])
        availability[snap["snapshot_date"]] = {
            f: (f in elem_set) for f in key_fields
        }

    return {"fields": key_fields, "availability": availability}


def main():
    parser = argparse.ArgumentParser(
        description="ClinVar XML Structure Analyzer for Q1/Q2/Q6/Q7"
    )
    parser.add_argument(
        "--snapshot-dir",
        default="./clinvar_cache/ftp_snapshots",
        help="Directory containing .xml.gz snapshots",
    )
    parser.add_argument(
        "--output", "-o",
        default="./xml_structure_analysis.json",
        help="Output JSON file path",
    )
    parser.add_argument(
        "--max-records", "-n",
        type=int, default=3000,
        help="Max records to sample per snapshot (default: 3000)",
    )
    args = parser.parse_args()

    snapshot_dir = Path(args.snapshot_dir)
    if not snapshot_dir.exists():
        print(f"ERROR: Snapshot directory not found: {snapshot_dir}")
        sys.exit(1)

    # Find all XML snapshots
    xml_files = sorted(
        [f for f in snapshot_dir.iterdir()
         if f.name.endswith(".xml.gz") and f.stat().st_size > 1000],
        key=lambda f: f.name,
    )

    if not xml_files:
        print(f"ERROR: No .xml.gz files found in {snapshot_dir}")
        sys.exit(1)

    print(f"Found {len(xml_files)} XML snapshots in {snapshot_dir}")
    print(f"Sampling {args.max_records} records per snapshot\n")

    # Analyze each snapshot
    snapshots_data = []
    for xml_path in xml_files:
        result = analyze_snapshot(xml_path, max_records=args.max_records)
        if result:
            snapshots_data.append(result)

    if not snapshots_data:
        print("ERROR: No snapshots could be analyzed")
        sys.exit(1)

    # Sort by date
    snapshots_data.sort(key=lambda x: x["snapshot_date"])

    # Compute derived analyses
    lifecycle = compute_lifecycle(snapshots_data)
    field_availability = compute_field_availability(snapshots_data)

    # Build output
    output = {
        "analysis_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "num_snapshots": len(snapshots_data),
        "snapshots": snapshots_data,
        "element_lifecycle": lifecycle,
        "field_availability": field_availability,
    }

    output_path = Path(args.output)
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nAnalysis complete: {len(snapshots_data)} snapshots")
    print(f"Output saved to: {output_path}")
    print(f"\nCopy to local machine:")
    print(f"  scp <user>@<server>:{output_path.resolve()} ~/Downloads/")


if __name__ == "__main__":
    main()
