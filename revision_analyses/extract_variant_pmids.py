#!/usr/bin/env python3
"""
Run this ON THE SERVER (aqua) to export per-variant PubMed IDs from the latest checkpoint.

It reads the most recent checkpoint in clinvar_output/checkpoints/ (a gzip pickle of
{'snapshot_date', 'data': {variation_id -> record}} where each record has 'all_pmids'),
and writes a small CSV of Variation_ID, VCV_ID, Gene, PMIDs for variants with >=1 PMID.

Usage (from /public/variant_reclassification/clinvar_output):
    python3 extract_variant_pmids.py \
        --checkpoints ./checkpoints \
        --out variant_pmids_latest.csv

Then copy variant_pmids_latest.csv into Daniel's Dropbox/ClinVar_GenomMedicine/ folder.
No third-party packages needed (stdlib only).
"""
import argparse, csv, glob, gzip, os, pickle, re, sys

def pick_latest(checkpoint_dir):
    files = glob.glob(os.path.join(checkpoint_dir, "*.checkpoint.pkl.gz"))
    if not files:
        sys.exit(f"No *.checkpoint.pkl.gz found in {checkpoint_dir}")
    def ymkey(p):
        m = re.search(r"(\d{4})-(\d{2})", os.path.basename(p))
        return (int(m.group(1)), int(m.group(2))) if m else (0, 0)
    files.sort(key=ymkey)
    return files[-1]                      # newest by YYYY-MM in filename

def gene_of(rec):
    for k in ("gene_symbol", "gene", "Gene", "GeneSymbol"):
        if rec.get(k):
            return rec[k]
    return ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoints", default="./checkpoints")
    ap.add_argument("--out", default="variant_pmids_latest.csv")
    ap.add_argument("--checkpoint-file", default=None,
                    help="explicit checkpoint path (overrides auto-latest)")
    a = ap.parse_args()

    path = a.checkpoint_file or pick_latest(a.checkpoints)
    print(f"Loading checkpoint: {path}", flush=True)
    with gzip.open(path, "rb") as f:
        cp = pickle.load(f)
    date = cp.get("snapshot_date", "?")
    data = cp["data"]
    print(f"  snapshot_date={date}  variants={len(data):,}", flush=True)

    n_written = 0
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Variation_ID", "VCV_ID", "Gene", "PMIDs"])
        for vid, rec in data.items():
            pmids = rec.get("all_pmids") or []
            if not pmids:
                continue
            w.writerow([rec.get("variation_id", vid), rec.get("vcv_id", ""),
                        gene_of(rec), ";".join(str(p) for p in pmids)])
            n_written += 1
    print(f"Wrote {n_written:,} variants with >=1 PMID to {a.out}  (snapshot {date})")

if __name__ == "__main__":
    main()
