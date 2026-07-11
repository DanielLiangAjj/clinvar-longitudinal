#!/usr/bin/env python
"""
Assemble the public GitHub repo for the ClinVar longitudinal analysis.
- Anonymizes ALL submitter identities to 'Laboratory A/B/C...' by submission rank (to match the
  manuscript), verifying that no real lab name leaks.
- Copies analysis scripts, small result JSONs, revision-analysis scripts+outputs, and figures.
- EXCLUDES the multi-GB raw CSVs, checkpoints, the raw NCBI organization_summary.txt, and the
  manuscript/response documents.
- Stages the interactive site's data (docs/data), converting a few CSVs to JSON.
"""
import json, csv, re, shutil, string
from pathlib import Path

SRC = Path("/Users/danielliang/Downloads/ClinVar Research Project")
REV = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses")
REPO = Path("/Users/danielliang/projects/clinvar-longitudinal")
if REPO.exists():
    shutil.rmtree(REPO)
for sub in ["scripts", "data", "data/submitter_analysis", "revision_analyses",
            "revision_analyses/outputs", "figures", "docs", "docs/data"]:
    (REPO / sub).mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- anonymization
REAL_NAME_FRAGMENTS = ["invitae", "ambry", "labcorp", "genedx", "illumina", "myriad", "counsyl",
                       "blueprint", "fulgent", "quest", "mayo", "prevention genetics", "centogene",
                       "baylor", "omim", "genereviews", "clingen", "emory", "partners", "arup",
                       "mgz", "3billion", "natera", "color", "sema4", "revvity", "greenwood"]

def letter_code(i):  # 0->A ... 25->Z, 26->AA ...
    s = ""
    i += 1
    while i > 0:
        i, r = divmod(i - 1, 26)
        s = string.ascii_uppercase[r] + s
    return f"Laboratory {s}"

def collect_names():
    """Union of every real submitter name across the submitter files, ranked by submissions."""
    counts = {}
    for f in ["data/submitter_overall.json", "data/submitter_analysis/submitter_overall.json"]:
        for row in json.load(open(SRC / f)):
            counts[row["submitter"]] = max(counts.get(row["submitter"], 0), row.get("total_submissions", 0))
    # sweep the other files for any additional names
    def add(n):
        if n and n not in counts:
            counts[n] = counts.get(n, 0)
    cp = json.load(open(SRC / "data/submitter_analysis/submitter_conflict_pairs.json"))
    for p in cp["top_50_conflict_pairs"]:
        add(p["submitter_a"]); add(p["submitter_b"])
    gm = json.load(open(SRC / "data/submitter_analysis/submitter_gene_matrix.json"))
    for n in gm["submitters"]:
        add(n)
    for f in ["data/submitter_analysis/submitter_temporal.json", "data/submitter_temporal.json"]:
        for n in json.load(open(SRC / f))["submitters"].keys():
            add(n)
    for row in json.load(open(SRC / "data/submitter_analysis/submitter_entry_exit.json"))["submitters"]:
        add(row["submitter"])
    rc = json.load(open(SRC / "data/submitter_analysis/submitter_reclassification.json"))
    for k in ["top_30_new_submitter_triggers", "top_30_reinterpretation"]:
        for row in rc.get(k, []):
            add(row["submitter"])
    return counts

counts = collect_names()
# rank by submissions desc, then name; assign codes
ranked = sorted(counts, key=lambda n: (-counts[n], n))
NAME2CODE = {n: letter_code(i) for i, n in enumerate(ranked)}
# unify obvious aliases of the top lab (Labcorp/Invitae) to one code
labcorp_code = next((c for n, c in NAME2CODE.items() if "invitae" in n.lower() or "labcorp" in n.lower()), None)
for n in list(NAME2CODE):
    low = n.lower()
    if labcorp_code and ("invitae" in low or low == "labcorp"):
        NAME2CODE[n] = labcorp_code
# longest names first for safe substring replacement
NAMES_SORTED = sorted(NAME2CODE, key=len, reverse=True)

def anonymize(obj):
    if isinstance(obj, dict):
        return {anonymize_str(k): anonymize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [anonymize(x) for x in obj]
    if isinstance(obj, str):
        return anonymize_str(obj)
    return obj

def anonymize_str(s):
    for name in NAMES_SORTED:
        if name in s:
            s = s.replace(name, NAME2CODE[name])
    return s

# ---------------------------------------------------------------- copy data (JSON)
SUBMITTER_FILES = {
    "data/submitter_overall.json", "data/submitter_temporal.json", "data/submitter_stats.json",
    "data/submitter_analysis/submitter_overall.json", "data/submitter_analysis/submitter_conflict_pairs.json",
    "data/submitter_analysis/submitter_gene_matrix.json", "data/submitter_analysis/submitter_reclassification.json",
    "data/submitter_analysis/submitter_temporal.json", "data/submitter_analysis/submitter_stats.json",
    "data/submitter_analysis/submitter_entry_exit.json",
}
EXCLUDE = {"data/organization_summary.txt", "data/submission_summary.txt.gz",
           "data/acmg_criteria_detailed.csv", "data/reclassification_events_detail.csv"}

# public cross-reference databases to keep in xml_structure_analysis.json's xref_databases;
# lab-registered XRef sources (everything else) are bucketed as 'Other' so no lab name leaks.
XREF_WHITELIST = {"omim", "dbsnp", "dbvar", "gene", "ncbi", "ncbi curation", "orphanet",
    "sequence ontology", "genereviews", "snomed ct", "medgen", "hgnc", "uniprotkb", "uniprot",
    "clingen", "clinvar", "mondo", "human phenotype ontology", "hpo", "pharmgkb", "cosmic",
    "office of rare diseases", "genetics home reference", "genbank", "ensembl", "refseq",
    "omim allelic variant", "gtr", "genetic testing registry", "clingen allele registry"}

def scrub_xref(data):
    for sn in data.get("snapshots", []):
        xr = sn.get("xref_databases")
        if isinstance(xr, dict):
            kept, other = {}, 0
            for k, v in xr.items():
                if k.lower() in XREF_WHITELIST:
                    kept[k] = v
                else:
                    other += v if isinstance(v, (int, float)) else 0
            if other:
                kept["Other sources"] = other
            sn["xref_databases"] = kept
    return data

copied = []
for f in SRC.glob("data/**/*"):
    rel = f.relative_to(SRC).as_posix()
    if f.is_dir() or rel in EXCLUDE:
        continue
    dst = REPO / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if rel == "data/xml_structure_analysis.json":
        data = scrub_xref(anonymize(json.load(open(f))))
        json.dump(data, open(dst, "w"))
    elif rel in SUBMITTER_FILES:
        data = json.load(open(f))
        json.dump(anonymize(data), open(dst, "w"))
    elif f.suffix == ".json":
        shutil.copy(f, dst)
    elif f.suffix == ".csv":
        # small CSVs only (<5MB); anonymize if they contain names
        if f.stat().st_size < 5_000_000:
            txt = f.read_text(errors="ignore")
            txt = anonymize_str(txt)
            dst.write_text(txt)
    copied.append(rel)

# ---------------------------------------------------------------- scripts, revision, figures
for f in SRC.glob("scripts/*.py"):
    shutil.copy(f, REPO / "scripts" / f.name)
for f in REV.glob("*.py"):
    shutil.copy(f, REPO / "revision_analyses" / f.name)
for f in REV.glob("outputs/*"):
    if f.is_file() and f.stat().st_size < 5_000_000:  # skip variant_pmids_latest.csv (94MB)
        shutil.copy(f, REPO / "revision_analyses/outputs" / f.name)
for f in list(SRC.glob("figures/*.png")) + list(REV.glob("figures/*.png")):
    shutil.copy(f, REPO / "figures" / f.name)

# ---------------------------------------------------------------- VERIFY no real names leak
# Allowed benign matches: SEMA4* is a gene (not the Sema4 lab); OMIM/ClinGen/GeneReviews are
# public resources openly named as gene-list / method sources (as in the manuscript).
# SEMA4* is a gene; OMIM/ClinGen/GeneReviews are public resources the manuscript names openly
# (they are anonymized where they appear as ClinVar SUBMITTERS, but are fine as resource refs).
PUBLIC_RESOURCES = {"omim", "clingen", "genereviews"}
DETECT = [f for f in REAL_NAME_FRAGMENTS if f not in ({"sema4"} | PUBLIC_RESOURCES)]
ALLOW_IN = {}
leaked = []
for f in REPO.rglob("*"):
    if f.is_file() and f.suffix in (".json", ".csv", ".txt"):
        rel = f.relative_to(REPO).as_posix()
        low = f.read_text(errors="ignore").lower()
        for frag in DETECT:
            if frag in low and frag not in ALLOW_IN.get(rel, set()):
                leaked.append((rel, frag))
print(f"code map: {len(NAME2CODE)} submitter names -> Laboratory codes (top: {ranked[0][:40]!r} -> {NAME2CODE[ranked[0]]})")
print(f"files copied to repo data/: {len([c for c in copied])}")
if leaked:
    print("!!! LEAK DETECTED (real name fragment found in published file):")
    for f, frag in leaked[:20]:
        print(f"    {f}: {frag}")
else:
    print("VERIFY: no real lab-name fragments found in any published data file ✓")
print(f"repo staged at: {REPO}")
