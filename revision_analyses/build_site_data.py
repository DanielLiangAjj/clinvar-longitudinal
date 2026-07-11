#!/usr/bin/env python
"""Assemble slim, safe JSON datasets for the interactive site into docs/data/.
Sources everything from the ALREADY-ANONYMIZED repo files plus the aggregate figure cache
(counts only, no names). Verifies no commercial-lab name leaks into docs/data."""
import json, csv, shutil
from pathlib import Path

REPO = Path("/Users/danielliang/projects/clinvar-longitudinal")
SRC = Path("/Users/danielliang/Downloads/ClinVar Research Project")
DD = REPO / "docs" / "data"
DD.mkdir(parents=True, exist_ok=True)

def csv_rows(p):
    with open(p) as f:
        return list(csv.DictReader(f))

def dump(name, obj):
    json.dump(obj, open(DD / name, "w"))

# --- aggregate figure cache (growth, classification, transition matrix) : counts only
for c in ["summary_agg.json", "history_agg.json"]:
    src = SRC / "figures" / ".cache" / c
    if src.exists():
        shutil.copy(src, DD / c)

# --- copy safe JSONs straight from the anonymized repo
for name in ["tiered_reclassification.json", "survival_overall.json", "stats_univariate.json",
             "its_results.json", "citation_filtering.json", "citation_distribution_by_domain.json"]:
    src = REPO / "revision_analyses" / "outputs" / name
    if src.exists():
        shutil.copy(src, DD / name)
for name in ["domain_summary.json", "reclassification_drivers.json", "batch_event_decomposition.json",
             "concordance_timeseries_v2.json", "submitter_overall.json"]:
    src = REPO / "data" / name
    if src.exists():
        shutil.copy(src, DD / name)
shutil.copy(REPO / "data" / "submitter_analysis" / "submitter_conflict_pairs.json", DD / "submitter_conflict_pairs.json")

# --- convert CSVs -> JSON
km = csv_rows(REPO / "revision_analyses" / "outputs" / "km_curve_overall.csv")
dump("km_curve.json", [{"yr": float(r["age_years"]), "ci": float(r["cuminc"]),
                        "lo": float(r["cuminc_lo"]), "hi": float(r["cuminc_hi"])} for r in km])
mv = csv_rows(REPO / "revision_analyses" / "outputs" / "monthly_timeseries.csv")
dump("monthly_timeseries.json", [{"m": r["month"], "vus": float(r["pct_VUS"]),
     "conc": float(r["pct_concordant_of_multi"]), "conf": float(r["pct_conflicting_of_all"]),
     "n": int(r["n_variants"]), "multi": int(r["n_multi"])} for r in mv])
reg = csv_rows(REPO / "revision_analyses" / "outputs" / "stats_multivariable_baseline.csv")
dump("stats_multivariable.json", [{"model": r["model"], "term": r["term"], "OR": float(r["OR"]),
     "lo": float(r["OR_lo"]), "hi": float(r["OR_hi"]), "p": r["p"]} for r in reg
     if r["term"] != "Intercept"])
sev = csv_rows(REPO / "data" / "conflict_severity_stratification.csv")
dump("conflict_severity.json", [{"sev": r["severity"], "events": int(r["count"]),
     "variants": int(r["unique_variants"])} for r in sev])

# --- slim schema evolution (per-snapshot metrics only; xref already whitelisted)
xs = json.load(open(REPO / "data" / "xml_structure_analysis.json"))
dump("schema_evolution.json", [{"date": s["snapshot_date"], "format": s.get("format"),
     "elements": s.get("num_unique_elements"), "attributes": s.get("total_unique_attributes"),
     "max_depth": s.get("max_depth"), "n_classifications": s.get("num_distinct_classifications"),
     "xref": s.get("xref_databases", {})} for s in xs["snapshots"]])

# --- verify no commercial-lab names in docs/data
FRAGS = ["invitae", "ambry", "genedx", "fulgent", "counsyl", "centogene", "natera"]
bad = []
for f in DD.glob("*.json"):
    low = f.read_text(errors="ignore").lower()
    for fr in FRAGS:
        if fr in low:
            bad.append((f.name, fr))
print("site data files:", len(list(DD.glob('*.json'))))
print("LEAK CHECK:", "FAIL " + str(bad) if bad else "clean ✓")
print("docs/data size:", sum(f.stat().st_size for f in DD.glob('*.json')) // 1024, "KB")
