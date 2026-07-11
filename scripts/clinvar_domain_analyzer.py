#!/usr/bin/env python3
"""
clinvar_domain_analyzer.py
--------------------------
Disease-domain subgroup analysis on ClinVar variant data.

Reads a summary CSV produced by the main ClinVar pipeline, classifies each
variant's gene into a disease domain (loaded from disease_domains.json), and
computes per-domain statistics including reclassification rates, oscillation
counts, pathogenic-to-benign / benign-to-pathogenic direction counts, and
evidence-level metrics.

Outputs:
    domain_summary.json  -- aggregate statistics per domain and cancer vs
                            noncancer binary split
    domain_per_gene.csv  -- per-gene row for the top 200 genes across all
                            domains

Usage:
    python clinvar_domain_analyzer.py \
        --summary-csv '/public/variant_reclassification/clinvar_output/summary_*.csv' \
        --domain-map  disease_domains.json \
        --output-dir  /public/variant_reclassification/
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from typing import Any, Dict, List, Set, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROGRESS_INTERVAL = 500_000

# Change-type patterns that indicate pathogenic-to-benign reclassification.
# We check whether any of these substrings appear inside the Change_Type field
# (which may contain multi-hop patterns like "P->VUS->LB").
P_TO_B_PATTERNS = ("P->B", "P->LB", "LP->B", "LP->LB",
                    "P/LP->B", "P/LP->LB", "P/LP->B/LB")

# Benign-to-pathogenic direction.
B_TO_P_PATTERNS = ("B->P", "B->LP", "LB->P", "LB->LP",
                    "B/LB->P", "B/LB->LP", "B/LB->P/LP")


# ---------------------------------------------------------------------------
# Domain-mapping helpers
# ---------------------------------------------------------------------------

def load_domain_map(path: str) -> Tuple[Dict[str, str], Set[str], Dict[str, List[str]], Set[str]]:
    """Load disease_domains.json and build derived lookups.

    Returns
    -------
    gene_to_domain : dict[str, str]
        Reverse mapping from gene symbol to its domain name.
    cancer_genes : set[str]
        Set of gene symbols classified as cancer genes.
    domain_gene_lists : dict[str, list[str]]
        Original domain -> gene-list mapping (for JSON output).
    common_disease_genes : set[str]
        Set of gene symbols primarily associated with common (non-rare/Mendelian)
        diseases, pharmacogenomics, or common carrier screening conditions.
    """
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    cancer_genes: Set[str] = set(data.get("cancer_genes", []))
    domain_defs: Dict[str, List[str]] = data.get("disease_domains", {})
    common_disease_genes: Set[str] = set(data.get("common_disease_genes", []))

    gene_to_domain: Dict[str, str] = {}
    for domain_name, gene_list in domain_defs.items():
        for gene in gene_list:
            # First domain wins if a gene appears in multiple domains.
            if gene not in gene_to_domain:
                gene_to_domain[gene] = domain_name

    return gene_to_domain, cancer_genes, domain_defs, common_disease_genes


def resolve_domain(gene: str, gene_to_domain: Dict[str, str]) -> str:
    """Return the domain for *gene*, defaulting to 'Other/Unclassified'."""
    if not gene or gene.strip() == "":
        return "Other/Unclassified"
    return gene_to_domain.get(gene.strip(), "Other/Unclassified")


def resolve_cancer_label(gene: str, cancer_genes: Set[str]) -> str:
    """Return 'cancer' or 'noncancer' for the binary split."""
    if not gene or gene.strip() == "":
        return "noncancer"
    return "cancer" if gene.strip() in cancer_genes else "noncancer"


def resolve_rarity(gene: str, common_genes: Set[str]) -> str:
    """Return 'common' if *gene* is primarily associated with common diseases,
    pharmacogenomics, or common carrier screening; otherwise 'rare'.

    In ClinVar, ~95% of variants are for rare/Mendelian conditions, so the
    default is 'rare' unless the gene appears in the common_disease_genes list.
    """
    if not gene or gene.strip() == "":
        return "rare"
    return "common" if gene.strip() in common_genes else "rare"


# ---------------------------------------------------------------------------
# Accumulator class
# ---------------------------------------------------------------------------

class DomainAccumulator:
    """Collects running totals for a single domain (or binary bucket).

    All numeric sums are kept as integers or floats so that memory stays
    constant regardless of CSV size -- we never store individual rows.
    """

    def __init__(self) -> None:
        self.total: int = 0
        self.changed: int = 0
        self.oscillating: int = 0
        self.p_to_b: int = 0
        self.b_to_p: int = 0
        self.sum_submissions: int = 0
        self.sum_pmids: int = 0
        self.sum_acmg: int = 0
        self.classification_dist: Counter = Counter()
        self.change_type_dist: Counter = Counter()
        self.gene_counts: Counter = Counter()
        self.review_status_dist: Counter = Counter()

    def ingest(self, row: Dict[str, str]) -> None:
        """Update running totals with one CSV row."""
        self.total += 1

        # Classification changed
        if row.get("Classification_Changed", "").strip() == "True":
            self.changed += 1

        # Oscillation
        if row.get("Had_Oscillation", "").strip() == "True":
            self.oscillating += 1

        # Change-type direction detection
        change_type = row.get("Change_Type", "").strip()
        if change_type:
            self.change_type_dist[change_type] += 1
            if any(pat in change_type for pat in P_TO_B_PATTERNS):
                self.p_to_b += 1
            if any(pat in change_type for pat in B_TO_P_PATTERNS):
                self.b_to_p += 1

        # Classification distribution
        cat = row.get("Current_Category", "").strip()
        if cat:
            self.classification_dist[cat] += 1

        # Numeric evidence columns (tolerant of missing / non-numeric)
        self.sum_submissions += _safe_int(row.get("Num_Submissions"))
        self.sum_pmids += _safe_int(row.get("Num_PMIDs"))
        self.sum_acmg += _safe_int(row.get("Num_ACMG_Codes"))

        # Gene counts
        gene = (row.get("Gene") or "").strip()
        if gene:
            self.gene_counts[gene] += 1

        # Review status
        rs = (row.get("Review_Status") or "").strip()
        if rs:
            self.review_status_dist[rs] += 1

    # ------------------------------------------------------------------
    # Serialisation helpers
    # ------------------------------------------------------------------

    def top_change_types(self, n: int = 10) -> List[Dict[str, Any]]:
        return [{"type": ct, "count": c}
                for ct, c in self.change_type_dist.most_common(n)]

    def top_genes(self, n: int = 20) -> List[Dict[str, Any]]:
        return [{"gene": g, "count": c}
                for g, c in self.gene_counts.most_common(n)]

    def as_dict(self, include_genes: bool = True) -> Dict[str, Any]:
        """Return a JSON-serialisable dictionary of aggregated stats."""
        d: Dict[str, Any] = {
            "total": self.total,
            "changed": self.changed,
            "change_rate_pct": round(self.changed / self.total * 100, 4) if self.total else 0.0,
            "oscillating": self.oscillating,
            "classification_dist": dict(self.classification_dist),
            "top_change_types": self.top_change_types(10),
            "p_to_b": self.p_to_b,
            "b_to_p": self.b_to_p,
            "mean_submissions": round(self.sum_submissions / self.total, 4) if self.total else 0.0,
            "mean_pmids": round(self.sum_pmids / self.total, 4) if self.total else 0.0,
            "mean_acmg": round(self.sum_acmg / self.total, 4) if self.total else 0.0,
            "review_status_dist": dict(self.review_status_dist),
        }
        if include_genes:
            d["top_genes"] = self.top_genes(20)
        return d


# ---------------------------------------------------------------------------
# Per-gene accumulator (lightweight)
# ---------------------------------------------------------------------------

class GeneAccumulator:
    """Minimal per-gene running totals for domain_per_gene.csv output."""

    __slots__ = (
        "total", "changed", "sum_submissions", "sum_pmids",
        "cat_P", "cat_LP", "cat_VUS", "cat_LB", "cat_B", "cat_Other",
    )

    def __init__(self) -> None:
        self.total: int = 0
        self.changed: int = 0
        self.sum_submissions: int = 0
        self.sum_pmids: int = 0
        self.cat_P: int = 0
        self.cat_LP: int = 0
        self.cat_VUS: int = 0
        self.cat_LB: int = 0
        self.cat_B: int = 0
        self.cat_Other: int = 0

    def ingest(self, row: Dict[str, str]) -> None:
        self.total += 1
        if row.get("Classification_Changed", "").strip() == "True":
            self.changed += 1
        self.sum_submissions += _safe_int(row.get("Num_Submissions"))
        self.sum_pmids += _safe_int(row.get("Num_PMIDs"))

        cat = row.get("Current_Category", "").strip()
        if cat == "P":
            self.cat_P += 1
        elif cat == "LP":
            self.cat_LP += 1
        elif cat in ("VUS", "Conflicting"):
            self.cat_VUS += 1
        elif cat == "LB":
            self.cat_LB += 1
        elif cat == "B":
            self.cat_B += 1
        else:
            self.cat_Other += 1


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _safe_int(val: str | None) -> int:
    """Convert a string to int, returning 0 on failure."""
    if val is None:
        return 0
    try:
        return int(val.strip())
    except (ValueError, AttributeError):
        return 0


def resolve_csv_path(pattern: str) -> str:
    """Resolve a path that may contain a glob pattern.

    Returns the first match sorted alphabetically so that the most recent
    timestamped file (summary_YYYYMMDD_HHMMSS.csv) is deterministic.
    """
    matches = sorted(glob.glob(pattern))
    if not matches:
        print(f"ERROR: No files matched the pattern: {pattern}", file=sys.stderr)
        sys.exit(1)
    chosen = matches[0]
    if len(matches) > 1:
        print(f"INFO: {len(matches)} files matched; using first: {chosen}")
    return chosen


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_csv(
    csv_path: str,
    gene_to_domain: Dict[str, str],
    cancer_genes: Set[str],
    common_disease_genes: Set[str],
) -> Tuple[Dict[str, DomainAccumulator], Dict[str, DomainAccumulator], Dict[str, DomainAccumulator], Dict[str, GeneAccumulator], int]:
    """Stream through the summary CSV and populate accumulators.

    Returns
    -------
    domain_accs : dict[str, DomainAccumulator]
        Keyed by domain name.
    binary_accs : dict[str, DomainAccumulator]
        Keyed by 'cancer' or 'noncancer'.
    rarity_accs : dict[str, DomainAccumulator]
        Keyed by 'rare' or 'common' (rare/Mendelian vs common disease).
    gene_accs : dict[str, GeneAccumulator]
        Keyed by gene symbol.
    total_rows : int
        Total rows processed.
    """
    domain_accs: Dict[str, DomainAccumulator] = defaultdict(DomainAccumulator)
    binary_accs: Dict[str, DomainAccumulator] = defaultdict(DomainAccumulator)
    rarity_accs: Dict[str, DomainAccumulator] = defaultdict(DomainAccumulator)
    gene_accs: Dict[str, GeneAccumulator] = defaultdict(GeneAccumulator)

    total_rows = 0

    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)

        for row in reader:
            total_rows += 1

            gene = (row.get("Gene") or "").strip()
            domain = resolve_domain(gene, gene_to_domain)
            cancer_label = resolve_cancer_label(gene, cancer_genes)
            rarity_label = resolve_rarity(gene, common_disease_genes)

            # Feed domain accumulator
            domain_accs[domain].ingest(row)

            # Feed binary accumulator
            binary_accs[cancer_label].ingest(row)

            # Feed rarity accumulator
            rarity_accs[rarity_label].ingest(row)

            # Feed per-gene accumulator
            gene_key = gene if gene else "(no_gene)"
            gene_accs[gene_key].ingest(row)

            # Progress logging
            if total_rows % PROGRESS_INTERVAL == 0:
                print(f"  ... processed {total_rows:,} rows")

    return dict(domain_accs), dict(binary_accs), dict(rarity_accs), dict(gene_accs), total_rows


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_domain_summary_json(
    output_dir: str,
    csv_path: str,
    total_rows: int,
    domain_accs: Dict[str, DomainAccumulator],
    binary_accs: Dict[str, DomainAccumulator],
    rarity_accs: Dict[str, DomainAccumulator],
    gene_to_domain: Dict[str, str],
) -> str:
    """Write domain_summary.json and return the output path."""

    result: Dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_csv": os.path.basename(csv_path),
        "total_variants_analyzed": total_rows,
        "cancer_vs_noncancer": {
            label: acc.as_dict(include_genes=False)
            for label, acc in sorted(binary_accs.items())
        },
        "rare_vs_common": {
            label: acc.as_dict(include_genes=True)
            for label, acc in sorted(rarity_accs.items())
        },
        "by_domain": {
            domain: acc.as_dict(include_genes=True)
            for domain, acc in sorted(domain_accs.items())
        },
        "gene_domain_map_used": {
            gene: domain for gene, domain in sorted(gene_to_domain.items())
        },
    }

    out_path = os.path.join(output_dir, "domain_summary.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2, ensure_ascii=False)

    return out_path


def write_domain_per_gene_csv(
    output_dir: str,
    gene_accs: Dict[str, GeneAccumulator],
    gene_to_domain: Dict[str, str],
    top_n: int = 200,
) -> str:
    """Write domain_per_gene.csv for the top *top_n* genes by variant count.

    Returns the output path.
    """
    # Sort genes by total variant count descending, take top N.
    ranked = sorted(gene_accs.items(), key=lambda kv: kv[1].total, reverse=True)[:top_n]

    out_path = os.path.join(output_dir, "domain_per_gene.csv")
    fieldnames = [
        "Gene", "Domain", "Total_Variants", "Changed_Variants",
        "Change_Rate_Pct", "P", "LP", "VUS", "LB", "B", "Other",
        "Mean_Submissions", "Mean_PMIDs",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for gene, acc in ranked:
            domain = resolve_domain(gene, gene_to_domain)
            change_rate = round(acc.changed / acc.total * 100, 4) if acc.total else 0.0
            mean_subs = round(acc.sum_submissions / acc.total, 4) if acc.total else 0.0
            mean_pmids = round(acc.sum_pmids / acc.total, 4) if acc.total else 0.0

            writer.writerow({
                "Gene": gene,
                "Domain": domain,
                "Total_Variants": acc.total,
                "Changed_Variants": acc.changed,
                "Change_Rate_Pct": change_rate,
                "P": acc.cat_P,
                "LP": acc.cat_LP,
                "VUS": acc.cat_VUS,
                "LB": acc.cat_LB,
                "B": acc.cat_B,
                "Other": acc.cat_Other,
                "Mean_Submissions": mean_subs,
                "Mean_PMIDs": mean_pmids,
            })

    return out_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Disease-domain subgroup analysis on ClinVar variant data. "
            "Classifies variants by gene into disease domains and computes "
            "per-domain statistics."
        ),
    )
    parser.add_argument(
        "--summary-csv",
        required=True,
        help=(
            "Path to the summary CSV file produced by the ClinVar pipeline. "
            "Supports glob patterns (e.g. 'summary_*.csv'); the first match "
            "is used."
        ),
    )
    parser.add_argument(
        "--domain-map",
        required=True,
        help="Path to disease_domains.json mapping genes to disease domains.",
    )
    parser.add_argument(
        "--output-dir",
        default="/public/variant_reclassification/",
        help=(
            "Directory where output files (domain_summary.json, "
            "domain_per_gene.csv) will be written. Default: "
            "/public/variant_reclassification/"
        ),
    )
    return parser.parse_args(argv)


def main(argv: List[str] | None = None) -> None:
    args = parse_args(argv)

    # --- Resolve inputs ---------------------------------------------------
    csv_path = resolve_csv_path(args.summary_csv)
    domain_map_path = args.domain_map
    output_dir = args.output_dir

    print(f"Summary CSV : {csv_path}")
    print(f"Domain map  : {domain_map_path}")
    print(f"Output dir  : {output_dir}")

    if not os.path.isfile(domain_map_path):
        print(f"ERROR: Domain map file not found: {domain_map_path}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # --- Load domain mapping ----------------------------------------------
    print("Loading disease domain mapping ...")
    gene_to_domain, cancer_genes, _, common_disease_genes = load_domain_map(domain_map_path)
    print(f"  {len(gene_to_domain):,} genes mapped to domains")
    print(f"  {len(cancer_genes):,} cancer genes defined")
    print(f"  {len(common_disease_genes):,} common disease genes defined")

    # --- Stream through CSV -----------------------------------------------
    print(f"Processing {csv_path} ...")
    domain_accs, binary_accs, rarity_accs, gene_accs, total_rows = process_csv(
        csv_path, gene_to_domain, cancer_genes, common_disease_genes,
    )
    print(f"Finished processing {total_rows:,} rows.")
    print(f"  Domains encountered: {sorted(domain_accs.keys())}")

    # --- Write outputs ----------------------------------------------------
    json_path = write_domain_summary_json(
        output_dir, csv_path, total_rows, domain_accs, binary_accs, rarity_accs,
        gene_to_domain,
    )
    print(f"Wrote domain summary  : {json_path}")

    csv_out_path = write_domain_per_gene_csv(
        output_dir, gene_accs, gene_to_domain, top_n=200,
    )
    print(f"Wrote per-gene CSV    : {csv_out_path}")

    # --- Quick sanity summary ---------------------------------------------
    print("\n--- Quick Summary ---")
    for label in ("cancer", "noncancer"):
        acc = binary_accs.get(label)
        if acc:
            print(
                f"  {label:>10s}: {acc.total:>10,} variants, "
                f"change rate {acc.changed / acc.total * 100:.2f}%, "
                f"oscillating {acc.oscillating:,}"
            )
    print()
    print("  Rare vs Common disease breakdown:")
    for label in ("rare", "common"):
        acc = rarity_accs.get(label)
        if acc:
            print(
                f"    {label:>7s}: {acc.total:>10,} variants, "
                f"change rate {acc.changed / acc.total * 100:.2f}%, "
                f"oscillating {acc.oscillating:,}"
            )
    print()
    for domain in sorted(domain_accs):
        acc = domain_accs[domain]
        print(
            f"  {domain:<25s}: {acc.total:>10,} variants, "
            f"change rate {acc.changed / acc.total * 100:.2f}%"
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
