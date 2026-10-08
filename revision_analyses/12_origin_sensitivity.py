#!/usr/bin/env python
"""
Revision analysis #12 — Snapshot-vintage audit (T1), origin extraction (T2), and
origin-stratified sensitivity analysis of the cancer-vs-non-cancer reclassification
comparison (T3).

Co-author's objection: the headline contrast (cancer genes reclassify 16.6% vs 8.2%)
uses the COSMIC Cancer Gene Census, a largely SOMATIC-driver list. If somatic-origin
submissions drive the excess reclassification, the germline interpretation of the
finding is confounded. Here we (a) audit whether the '2025-07' snapshot really is
July-2025 ClinVar content, (b) extract per-SCV allele origin, and (c) re-estimate the
contrast excluding somatic-origin variants.

Methodology is mirrored from 02_statistics.py:
  - two_by_two()  -> 02_statistics.py:31-49 (Katz RR CI, Woolf OR CI, chi2 no correction)
  - changed flag  -> 02_statistics.py:57
  - is_cancer     -> 02_statistics.py:62

Outputs:
  outputs/variant_origin_map.csv.gz      (T2: variation_id -> origins, has_somatic, has_origin)
  outputs/snapshot_vintage_audit.json    (T1 evidence)
  outputs/origin_sensitivity.json        (T3 results)
"""
import gc, gzip, json, pickle, time, warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

DROPBOX = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/outputs")
CKPT = DROPBOX / "ClinVarVariationRelease_2025-07.checkpoint.pkl.gz"
SUMMARY = DROPBOX / "summary_20260219_024805.csv"
ORIGIN_MAP = OUT / "variant_origin_map.csv.gz"

# origin strings that denote a germline/heritable allele origin (used only for the
# supplementary STRICT definition; the primary definition follows the pre-registered
# spec: germline_only = has_origin AND NOT has_somatic)
GERMLINE_ISH = {"germline", "de novo", "maternal", "paternal", "inherited",
                "biparental", "uniparental"}
UNINFORMATIVE = {"unknown", "not provided", "not applicable", ""}


# ---------------------------------------------------------------- 02_statistics.py:31-49
def two_by_two(x1, n1, x2, n2, label1, label2):
    """x=changed, n=total in each group. Returns RR, OR, chi2 p with CIs.
    Verbatim port of 02_statistics.py:31-49 so results are directly comparable."""
    a, b = x1, n1 - x1
    c, d = x2, n2 - x2
    p1, p2 = a / n1, c / n2
    rr = p1 / p2
    se_lnrr = np.sqrt(b / (a * n1) + d / (c * n2))
    rr_lo, rr_hi = rr * np.exp(-1.96 * se_lnrr), rr * np.exp(1.96 * se_lnrr)
    orr = (a * d) / (b * c)
    se_lnor = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    or_lo, or_hi = orr * np.exp(-1.96 * se_lnor), orr * np.exp(1.96 * se_lnor)
    chi2, p, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=False)
    return {
        "group1": label1, "n1": int(n1), "changed1": int(a), "rate1_pct": round(100 * p1, 3),
        "group2": label2, "n2": int(n2), "changed2": int(c), "rate2_pct": round(100 * p2, 3),
        "risk_ratio": round(rr, 3), "rr_95ci": [round(rr_lo, 3), round(rr_hi, 3)],
        "odds_ratio": round(orr, 3), "or_95ci": [round(or_lo, 3), round(or_hi, 3)],
        "chi2": round(chi2, 1), "p_value": (f"{p:.2e}" if p > 0 else "<1e-300"),
    }


# ---------------------------------------------------------------- T1 + T2
def build_origin_map():
    """Load the checkpoint once: emit T1 vintage evidence and the T2 origin map."""
    t0 = time.time()
    print(f"loading checkpoint {CKPT.name} ...", flush=True)
    with gzip.open(CKPT, "rb") as fh:
        ck = pickle.load(fh)
    print(f"  loaded in {time.time()-t0:.1f}s", flush=True)

    wrapper = {k: ck[k] for k in ck if k != "data"}
    data = ck["data"]
    print(f"  wrapper: {wrapper}")
    print(f"  variants in checkpoint: {len(data):,}", flush=True)

    # ---- T1: content-vintage evidence over ALL records ----
    rec_created, rec_updated = Counter(), Counter()
    sub_submitted, sub_updated = Counter(), Counter()
    classification_vc, review_vc = Counter(), Counter()
    rec_keys, sub_keys = Counter(), Counter()

    # ---- T2: origin extraction ----
    vids, joined, has_som, has_org, strict_germ = [], [], [], [], []
    origin_scv_vc = Counter()      # per-SCV origin strings
    origin_variant_vc = Counter()  # per-variant distinct origin strings
    n_scv = 0
    n_scv_with_origin = 0

    for i, (vid, rec) in enumerate(data.items()):
        rec_keys.update(rec.keys())
        dc = rec.get("date_created") or ""
        du = rec.get("date_updated") or ""
        if dc: rec_created[dc[:7]] += 1
        if du: rec_updated[du[:7]] += 1
        classification_vc[str(rec.get("classification", ""))] += 1
        review_vc[str(rec.get("review_status", ""))] += 1

        origins = set()
        for s in rec.get("submissions", []) or []:
            n_scv += 1
            sub_keys.update(s.keys())
            ds = s.get("date_submitted") or ""
            dsu = s.get("date_updated") or ""
            if ds: sub_submitted[ds[:7]] += 1
            if dsu: sub_updated[dsu[:7]] += 1
            o = s.get("origin")
            o = "" if o is None else str(o).strip()
            origin_scv_vc[o] += 1
            if o:
                n_scv_with_origin += 1
                origins.add(o)

        for o in origins:
            origin_variant_vc[o] += 1

        # populated-origin set excludes the empty string only (spec: "populated origin")
        ho = len(origins) > 0
        hs = any("somatic" in o.lower() for o in origins)
        sg = (not hs) and any(o.lower() in GERMLINE_ISH for o in origins)

        vids.append(vid)
        joined.append("|".join(sorted(origins)))
        has_som.append(hs)
        has_org.append(ho)
        strict_germ.append(sg)

        if (i + 1) % 500000 == 0:
            print(f"    {i+1:,} variants scanned ({time.time()-t0:.0f}s)", flush=True)

    n_var = len(vids)
    om = pd.DataFrame({
        "variation_id": vids,
        "origins_joined": joined,
        "has_somatic": has_som,
        "has_origin": has_org,
        "strict_germline": strict_germ,
    })
    om.to_csv(ORIGIN_MAP, index=False, compression="gzip")
    print(f"  wrote {ORIGIN_MAP} ({len(om):,} rows)", flush=True)

    def top_months(c, n=14):
        return {k: v for k, v in sorted(c.items(), reverse=True)[:n]}

    audit = {
        "checkpoint_file": str(CKPT),
        "wrapper_metadata": wrapper,
        "n_variants": n_var,
        "n_scv": n_scv,
        "n_scv_with_populated_origin": n_scv_with_origin,
        "record_keys": dict(rec_keys),
        "submission_keys": dict(sub_keys),
        "date_created_max": max(rec_created) if rec_created else None,
        "date_updated_max": max(rec_updated) if rec_updated else None,
        "sub_date_submitted_max": max(sub_submitted) if sub_submitted else None,
        "sub_date_updated_max": max(sub_updated) if sub_updated else None,
        "date_created_min": min(rec_created) if rec_created else None,
        "date_updated_min": min(rec_updated) if rec_updated else None,
        "date_updated_top_months_desc": top_months(rec_updated),
        "sub_date_updated_top_months_desc": top_months(sub_updated),
        "date_created_top_months_desc": top_months(rec_created),
        "classification_value_counts_top": dict(classification_vc.most_common(25)),
        "review_status_value_counts_top": dict(review_vc.most_common(15)),
        "origin_value_counts_per_scv": dict(origin_scv_vc.most_common()),
        "origin_value_counts_per_variant": dict(origin_variant_vc.most_common()),
        "schema_split_fields_present": {
            "any_germline_classification_field": any(
                "germline" in k.lower() for k in rec_keys),
            "any_somatic_clinical_impact_field": any(
                "somatic" in k.lower() and "clinical" in k.lower() for k in rec_keys),
            "any_oncogenicity_field": any("oncogenic" in k.lower() for k in rec_keys),
            "per_scv_classification_field": any(
                k in sub_keys for k in ("classification", "clinical_significance",
                                        "interpretation", "significance")),
        },
    }
    json.dump(audit, open(OUT / "snapshot_vintage_audit.json", "w"), indent=2, default=str)
    print(f"  wrote {OUT/'snapshot_vintage_audit.json'}", flush=True)

    del ck, data
    gc.collect()
    return audit


def load_origin_map():
    if ORIGIN_MAP.exists():
        print(f"using cached {ORIGIN_MAP.name}", flush=True)
        return pd.read_csv(ORIGIN_MAP, dtype={"variation_id": str},
                           keep_default_na=False)
    build_origin_map()
    gc.collect()
    return pd.read_csv(ORIGIN_MAP, dtype={"variation_id": str}, keep_default_na=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    om = load_origin_map()
    for c in ("has_somatic", "has_origin", "strict_germline"):
        if om[c].dtype == object:
            om[c] = om[c].map({"True": True, "False": False}).astype(bool)
    print(f"origin map: {len(om):,} variants", flush=True)

    # ---- summary CSV, mirroring 02_statistics.py:53-64 ----
    print("loading summary ...", flush=True)
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "Classification_Changed"],
                     dtype={"Variation_ID": str, "Gene": str})
    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    N = len(df)
    print(f"  {N:,} variants, {df['changed'].sum():,} changed", flush=True)

    DOMAINS = json.load(open(PROJ / "disease_domains.json"))
    cancer = set(DOMAINS["cancer_genes"])
    df["is_cancer"] = df["Gene"].isin(cancer)

    # ---- SANITY GATE: reproduce the published baseline exactly ----
    def grp(mask):
        sub = df[mask]
        return int(sub["changed"].sum()), int(len(sub))

    xc, nc = grp(df["is_cancer"])
    xnc, nnc = grp(~df["is_cancer"])
    baseline = two_by_two(xc, nc, xnc, nnc, "Cancer (COSMIC CGC)", "Non-cancer")
    published = json.load(open(OUT / "stats_univariate.json"))["cancer_vs_noncancer"]
    gate = {k: (baseline[k], published[k]) for k in
            ("n1", "rate1_pct", "n2", "rate2_pct", "risk_ratio", "odds_ratio")}
    ok = all(baseline[k] == published[k] for k in gate)
    print("\n=== SANITY GATE (reproduce published baseline) ===")
    for k, (got, exp) in gate.items():
        flag = "OK " if got == exp else "MISMATCH"
        print(f"  {flag} {k:12s} got={got} published={exp}")
    if not ok:
        raise SystemExit("BASELINE NOT REPRODUCED — stopping per protocol.")
    print("  -> baseline reproduced exactly.")

    # ---- coverage ----
    df = df.merge(om, left_on="Variation_ID", right_on="variation_id", how="left")
    df["in_checkpoint"] = df["variation_id"].notna()
    df["has_somatic"] = df["has_somatic"].fillna(False).astype(bool)
    df["has_origin"] = df["has_origin"].fillna(False).astype(bool)
    df["strict_germline"] = df["strict_germline"].fillna(False).astype(bool)

    cov = {
        "summary_rows": int(N),
        "origin_map_rows": int(len(om)),
        "matched_in_checkpoint": int(df["in_checkpoint"].sum()),
        "matched_pct": round(100 * df["in_checkpoint"].mean(), 3),
        "not_in_checkpoint": int((~df["in_checkpoint"]).sum()),
        "with_populated_origin": int(df["has_origin"].sum()),
        "with_populated_origin_pct": round(100 * df["has_origin"].mean(), 3),
        "no_origin_data": int((~df["has_origin"]).sum()),
        "has_somatic": int(df["has_somatic"].sum()),
        "has_somatic_pct": round(100 * df["has_somatic"].mean(), 3),
        "germline_only_spec": int((df["has_origin"] & ~df["has_somatic"]).sum()),
        "strict_germline": int(df["strict_germline"].sum()),
        "reclass_rate_in_checkpoint_pct": round(
            100 * df.loc[df["in_checkpoint"], "changed"].mean(), 3),
        "reclass_rate_not_in_checkpoint_pct": round(
            100 * df.loc[~df["in_checkpoint"], "changed"].mean(), 3) if (~df["in_checkpoint"]).any() else None,
    }
    print("\n=== T3 COVERAGE ===")
    for k, v in cov.items():
        print(f"  {k:38s} {v}")

    # somatic prevalence by gene set
    cov["has_somatic_pct_among_cancer"] = round(
        100 * df.loc[df["is_cancer"], "has_somatic"].mean(), 3)
    cov["has_somatic_pct_among_noncancer"] = round(
        100 * df.loc[~df["is_cancer"], "has_somatic"].mean(), 3)
    print(f"  has_somatic among CGC-cancer genes:     {cov['has_somatic_pct_among_cancer']}%")
    print(f"  has_somatic among non-cancer genes:     {cov['has_somatic_pct_among_noncancer']}%")

    # ---- T3 scenarios ----
    res = {"baseline_published": published, "sanity_gate_passed": ok, "coverage": cov}
    scen = {}

    scen["a_all_variants"] = baseline

    m = ~df["has_somatic"]
    x1, n1 = grp(df["is_cancer"] & m); x2, n2 = grp(~df["is_cancer"] & m)
    scen["b_exclude_has_somatic"] = two_by_two(x1, n1, x2, n2,
                                              "Cancer (CGC), somatic-origin excluded",
                                              "Non-cancer, somatic-origin excluded")

    m = df["has_origin"] & ~df["has_somatic"]
    x1, n1 = grp(df["is_cancer"] & m); x2, n2 = grp(~df["is_cancer"] & m)
    scen["c_germline_only"] = two_by_two(x1, n1, x2, n2,
                                        "Cancer (CGC), germline_only",
                                        "Non-cancer, germline_only")

    # supplementary: restrict to checkpoint-covered variants only (removes the
    # coverage artefact that (b) silently keeps origin-less variants)
    m = df["in_checkpoint"]
    x1, n1 = grp(df["is_cancer"] & m); x2, n2 = grp(~df["is_cancer"] & m)
    scen["d_in_checkpoint_only"] = two_by_two(x1, n1, x2, n2,
                                             "Cancer (CGC), in 2025-07 checkpoint",
                                             "Non-cancer, in 2025-07 checkpoint")

    m = df["in_checkpoint"] & ~df["has_somatic"]
    x1, n1 = grp(df["is_cancer"] & m); x2, n2 = grp(~df["is_cancer"] & m)
    scen["e_in_checkpoint_exclude_somatic"] = two_by_two(x1, n1, x2, n2,
                                                        "Cancer (CGC), in ckpt, somatic excluded",
                                                        "Non-cancer, in ckpt, somatic excluded")

    m = df["strict_germline"]
    x1, n1 = grp(df["is_cancer"] & m); x2, n2 = grp(~df["is_cancer"] & m)
    scen["f_strict_germline"] = two_by_two(x1, n1, x2, n2,
                                          "Cancer (CGC), strict germline origin",
                                          "Non-cancer, strict germline origin")

    res["scenarios"] = scen
    json.dump(res, open(OUT / "origin_sensitivity.json", "w"), indent=2)

    print("\n=== T3 ORIGIN SENSITIVITY (cancer vs non-cancer) ===")
    print(f"{'scenario':38s} {'n_cancer':>9} {'rate1':>7} {'n_nonc':>10} {'rate2':>7} "
          f"{'RR':>6} {'RR 95% CI':>16} {'OR':>6} {'OR 95% CI':>16}")
    for k, v in scen.items():
        print(f"{k:38s} {v['n1']:>9,} {v['rate1_pct']:>6.2f}% {v['n2']:>10,} "
              f"{v['rate2_pct']:>6.2f}% {v['risk_ratio']:>6.3f} "
              f"[{v['rr_95ci'][0]:.3f},{v['rr_95ci'][1]:.3f}]".rjust(0)
              + f" {v['odds_ratio']:>6.3f} [{v['or_95ci'][0]:.3f},{v['or_95ci'][1]:.3f}]")
    print(f"\nwritten: {OUT/'origin_sensitivity.json'}")


if __name__ == "__main__":
    main()
