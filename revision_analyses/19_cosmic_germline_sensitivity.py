#!/usr/bin/env python
"""
Revision analysis #19 — COSMIC-own-germline sensitivity analysis.

Co-author objection: the "cancer genes reclassify more" result uses the 304 COSMIC
Cancer Gene Census (CGC) genes, which are MAJORITY SOMATIC drivers, and the germline
robustness check already run (script 13) uses a hand-curated 106-gene "Hereditary
Cancer" subset with no documented provenance.

Fix (this script): re-run the cancer-vs-non-cancer reclassification comparison using an
EXTERNALLY-PROVENANCED germline cancer-gene set — the CGC genes that COSMIC's own
Cancer Gene Census flags as GERMLINE (GERMLINE == 'y' in the v101 CGC TSV).

Methodology mirrored EXACTLY from (READ, not modified):
  - two_by_two()        -> 02_statistics.py:31-49  (Katz RR CI, Woolf OR CI, chi2 no correction)
  - changed flag        -> 02_statistics.py:57 / 12_origin_sensitivity.py:216
  - fit_adjusted() / M1 -> 13_cgc_sensitivity.py:72-94  (== 05_multivariable_baseline.py:82-84)
  - baseline covariates -> 13_cgc_sensitivity.py:200-227 (== 05_multivariable_baseline.py:58-68)
  - origin map / has_somatic -> 12_origin_sensitivity.py (variant_origin_map.csv.gz)

Baseline sanity gates reproduced before trusting anything:
  crude CGC(304): rate1 16.619%, rate2 8.191%, RR 2.029, OR 2.234  (stats_univariate.json)
  adjusted CGC(304): OR 1.484 [1.470, 1.499]                       (origin_cgc_sensitivity_summary.json)

Outputs:
  outputs/cosmic_germline_sensitivity.json
  outputs/cosmic_germline_sensitivity.md
"""
import gc, json, re, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

DROPBOX = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/outputs")
SUMMARY = DROPBOX / "summary_20260219_024805.csv"
ORIGIN_MAP = OUT / "variant_origin_map.csv.gz"
CGC_TSV = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/inputs/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/"
               "Cosmic_CancerGeneCensus_v101_GRCh37.tsv")
CGC_SUMMARY = OUT / "origin_cgc_sensitivity_summary.json"


# ---------------------------------------------------------------- 02_statistics.py:31-49
def two_by_two(x1, n1, x2, n2, label1, label2):
    """Verbatim port of 02_statistics.py:31-49 so results are directly comparable."""
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


def vclass(s):
    """05_multivariable_baseline.py:28-32"""
    s = str(s).lower()
    for k in ["missense", "nonsense", "frameshift", "splice", "synonymous", "inframe"]:
        if k in s:
            return k
    return "other"


def fit_adjusted(d, exposure, tag):
    """Re-fit M1 from 05_multivariable_baseline.py:82-84 (== 13_cgc_sensitivity.py:72-94)
    with `exposure` swapped in for is_cancer. Returns the exposure adjusted OR + 95% CI."""
    f = (f"changed ~ {exposure} + log_base_subs + follow_up_yr + log_gene_span "
         "+ C(era, Treatment('RCV_2015_2019')) + C(vclass, Treatment('missense')) "
         "+ C(base_review_tier, Treatment('single'))")
    m = smf.logit(f, data=d).fit(disp=0, maxiter=200)
    p, ci, pv = m.params, m.conf_int(), m.pvalues
    out = {
        "model": tag, "term": exposure, "n": int(m.nobs),
        "n_exposed": int(d[exposure].sum()),
        "pseudo_r2": round(float(m.prsquared), 4),
        "OR": round(float(np.exp(p[exposure])), 3),
        "OR_lo": round(float(np.exp(ci.loc[exposure, 0])), 3),
        "OR_hi": round(float(np.exp(ci.loc[exposure, 1])), 3),
        "p": f"{pv[exposure]:.2e}",
    }
    print(f"  {tag:44s} n={out['n']:>9,} exposed={out['n_exposed']:>8,} "
          f"OR={out['OR']:.3f} [{out['OR_lo']:.3f}, {out['OR_hi']:.3f}] p={out['p']}", flush=True)
    del m
    gc.collect()
    return out


# a synonym token counts as a gene symbol only if it looks like an HGNC-style symbol
# (letters/digits/hyphen, starts with a letter), NOT an Ensembl/UniProt/Entrez id.
_ENSEMBL = re.compile(r"^ENSG\d+", re.I)
_UNIPROT = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")


def symbol_like(tok):
    tok = tok.strip()
    if not tok or tok.isdigit():           # Entrez id / empty
        return False
    if _ENSEMBL.match(tok) or _UNIPROT.match(tok):
        return False
    if not re.match(r"^[A-Za-z]", tok):     # must start with a letter
        return False
    return True


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOMAINS = json.load(open(PROJ / "disease_domains.json"))
    cgc304 = set(DOMAINS["cancer_genes"])                      # paper's 304-gene CGC set
    hc106 = set(DOMAINS["disease_domains"]["Hereditary Cancer"])

    # ---------------- T1: build the germline-CGC gene set ----------------
    print("=== T1: build germline-CGC set from COSMIC's own GERMLINE flag ===", flush=True)
    tsv = pd.read_csv(CGC_TSV, sep="\t", dtype=str)
    print(f"  CGC TSV rows (genes): {len(tsv):,}; GERMLINE vc: "
          f"{tsv['GERMLINE'].value_counts().to_dict()}")
    gdf = tsv[tsv["GERMLINE"].str.lower() == "y"].copy()
    germline_symbols = sorted(gdf["GENE_SYMBOL"].tolist())
    germline_set = set(germline_symbols)
    print(f"  germline-CGC genes (GERMLINE=='y'): {len(germline_set)}")

    # ---- gene-symbol matching to the summary CSV Gene column ----
    print("\nloading summary Gene column for matching + crude analysis ...", flush=True)
    crude = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "Classification_Changed"],
                        dtype={"Variation_ID": str, "Gene": str})
    crude["changed"] = crude["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    N = len(crude)
    dataset_genes = set(crude["Gene"].dropna().unique())
    print(f"  summary rows: {N:,}; distinct Gene symbols in dataset: {len(dataset_genes):,}")

    # synonym lookup for germline genes not directly present
    syn_map = {}   # germline GENE_SYMBOL -> synonym actually found in dataset
    for _, r in gdf.iterrows():
        toks = [t.strip() for t in str(r["SYNONYMS"]).split(",")]
        syn_map[r["GENE_SYMBOL"]] = [t for t in toks
                                     if symbol_like(t) and t != r["GENE_SYMBOL"]]

    matched_direct, matched_via_syn, unmatched = [], [], []
    resolved_symbols = set()          # symbols (as they appear in dataset) to flag as germline-CGC
    syn_resolution = {}               # germline gene -> synonym used
    for g in germline_symbols:
        if g in dataset_genes:
            matched_direct.append(g)
            resolved_symbols.add(g)
        else:
            hits = [s for s in syn_map[g] if s in dataset_genes]
            if hits:
                matched_via_syn.append(g)
                syn_resolution[g] = hits
                resolved_symbols.update(hits)
            else:
                unmatched.append(g)

    print(f"\n  matched directly (symbol present in dataset): {len(matched_direct)}")
    print(f"  matched only via COSMIC SYNONYMS:              {len(matched_via_syn)}  {syn_resolution}")
    print(f"  UNMATCHED (no symbol/synonym in dataset):      {len(unmatched)}  {unmatched}")

    overlap = {
        "germline_cgc_n": len(germline_set),
        "germline_cgc_genes": germline_symbols,
        "matched_direct_n": len(matched_direct),
        "matched_via_synonym_n": len(matched_via_syn),
        "matched_via_synonym": syn_resolution,
        "unmatched_n": len(unmatched),
        "unmatched_genes": unmatched,
        "n_effective_flag_symbols": len(resolved_symbols),
        "overlap_with_paper_304_cgc": len(germline_set & cgc304),
        "germline_cgc_not_in_paper_304": sorted(germline_set - cgc304),
        "overlap_with_106_hereditary_cancer": len(germline_set & hc106),
        "germline_cgc_not_in_106_hc": sorted(germline_set - hc106),
        "n_106_hc_not_germline_cgc": len(hc106 - germline_set),
    }
    print(f"\n  overlap with paper's 304 CGC: {overlap['overlap_with_paper_304_cgc']} / {len(germline_set)}"
          f"  (germline-CGC genes absent from paper 304: {len(overlap['germline_cgc_not_in_paper_304'])})")
    print(f"  overlap with 106 Hereditary Cancer: {overlap['overlap_with_106_hereditary_cancer']} / {len(germline_set)}"
          f"  (106-HC genes not germline-CGC: {overlap['n_106_hc_not_germline_cgc']})")

    # ---------------- T2: crude contrasts ----------------
    print("\n=== T2: crude germline-CGC vs all-other-genes ===", flush=True)
    om = pd.read_csv(ORIGIN_MAP, dtype={"variation_id": str}, keep_default_na=False)
    for c in ("has_somatic", "has_origin", "strict_germline"):
        if om[c].dtype == object:
            om[c] = om[c].map({"True": True, "False": False}).astype(bool)
    crude = crude.merge(om[["variation_id", "has_somatic"]], left_on="Variation_ID",
                        right_on="variation_id", how="left")
    crude["has_somatic"] = crude["has_somatic"].fillna(False).astype(bool)
    crude["is_gcgc"] = crude["Gene"].isin(resolved_symbols)
    crude["is_cancer"] = crude["Gene"].isin(cgc304)          # paper's 304-gene CGC set

    def grp(d, mask):
        s = d[mask]
        return int(s["changed"].sum()), int(len(s))

    crude_res = {}
    x1, n1 = grp(crude, crude["is_gcgc"]); x2, n2 = grp(crude, ~crude["is_gcgc"])
    crude_res["all_variants"] = two_by_two(
        x1, n1, x2, n2, "Germline-CGC (COSMIC GERMLINE=='y')", "All other genes")
    m = ~crude["has_somatic"]
    x1, n1 = grp(crude, crude["is_gcgc"] & m); x2, n2 = grp(crude, ~crude["is_gcgc"] & m)
    crude_res["exclude_somatic"] = two_by_two(
        x1, n1, x2, n2, "Germline-CGC, somatic-origin excluded",
        "All other genes, somatic-origin excluded")
    for k, v in crude_res.items():
        print(f"  {k:16s} exposed {v['rate1_pct']:.3f}% (n={v['n1']:,})  vs others "
              f"{v['rate2_pct']:.3f}% (n={v['n2']:,})  RR={v['risk_ratio']} {v['rr_95ci']}  "
              f"OR={v['odds_ratio']} {v['or_95ci']}")

    # ---- T4: germline-CGC vs TRUE non-cancer (comparator drops ALL Census/cancer genes:
    #         not in paper's 304 CGC AND not in germline-CGC; mirrors 13's HC-vs-true-non-cancer) ----
    print("\n=== T4: crude germline-CGC vs TRUE non-cancer (all CGC/Census genes dropped) ===",
          flush=True)
    comp = ~crude["is_cancer"] & ~crude["is_gcgc"]      # true non-cancer comparator
    tn_crude = {}
    x1, n1 = grp(crude, crude["is_gcgc"]); x2, n2 = grp(crude, comp)
    tn_crude["all_variants"] = two_by_two(
        x1, n1, x2, n2, "Germline-CGC (114)", "True non-cancer (no CGC/Census gene)")
    m = ~crude["has_somatic"]
    x1, n1 = grp(crude, crude["is_gcgc"] & m); x2, n2 = grp(crude, comp & m)
    tn_crude["exclude_somatic"] = two_by_two(
        x1, n1, x2, n2, "Germline-CGC (114), somatic excl.",
        "True non-cancer, somatic excl.")
    for k, v in tn_crude.items():
        print(f"  {k:16s} exposed {v['rate1_pct']:.3f}% (n={v['n1']:,})  vs true-non-cancer "
              f"{v['rate2_pct']:.3f}% (n={v['n2']:,})  RR={v['risk_ratio']} {v['rr_95ci']}  "
              f"OR={v['odds_ratio']} {v['or_95ci']}")

    del crude
    gc.collect()

    # ---------------- adjusted models ----------------
    print("\n=== T2 adjusted (mirrors 05/13 M1; baseline covariates) ===", flush=True)
    base = pd.read_csv(OUT / "baseline_covariates.csv",
                       dtype={"Variation_ID": str, "base_num_subs": "Int64"})
    print(f"  baseline rows: {len(base):,}")
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "First_Date",
                                       "Classification_Changed", "Molecular_Consequence",
                                       "Start", "Stop"],
                     dtype={"Variation_ID": str, "Gene": str, "First_Date": str})
    df = df.merge(base[["Variation_ID", "base_num_subs", "base_review_tier"]],
                  on="Variation_ID", how="inner")
    print(f"  merged rows: {len(df):,}")
    del base
    gc.collect()

    # covariates verbatim from 13_cgc_sensitivity.py:200-212 (== 05:58-68)
    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true").astype(int)
    df["log_base_subs"] = np.log10(df["base_num_subs"].fillna(1).clip(lower=1).astype(float))
    yr = df["First_Date"].str.slice(0, 4).astype(float)
    mo = df["First_Date"].str.slice(5, 7).astype(float)
    df["follow_up_yr"] = (2025 + 7 / 12) - (yr + mo / 12)
    df["era"] = pd.cut(yr, [0, 2014.99, 2019.99, 3000],
                       labels=["pre2015", "RCV_2015_2019", "VCV_2020plus"])
    gene_span = (df.groupby("Gene").apply(
        lambda g: g["Stop"].astype(float).max() - g["Start"].astype(float).min(),
        include_groups=False))
    df["log_gene_span"] = np.log10(df["Gene"].map(gene_span).clip(lower=1).fillna(1) + 1)
    df["vclass"] = df["Molecular_Consequence"].map(vclass)

    df["is_cancer"] = df["Gene"].isin(cgc304).astype(int)          # 304-CGC baseline exposure
    df["is_germline_cgc"] = df["Gene"].isin(resolved_symbols).astype(int)
    # genes in a Census/cancer set but NOT germline-flagged: dropped for the "vs true
    # non-cancer" contrast so the comparator holds only non-CGC, non-germline genes
    # (mirrors 13_cgc_sensitivity.py's is_cgc_not_hc, line 216/246)
    df["is_cgc_not_gcgc"] = ((df["is_cancer"] == 1) & (df["is_germline_cgc"] == 0)).astype(int)

    df = df.merge(om[["variation_id", "has_somatic"]], left_on="Variation_ID",
                  right_on="variation_id", how="left")
    df["has_somatic"] = df["has_somatic"].fillna(False).astype(bool)

    keep = ["changed", "log_base_subs", "follow_up_yr", "log_gene_span", "vclass",
            "base_review_tier", "is_cancer", "is_germline_cgc", "is_cgc_not_gcgc",
            "era", "has_somatic"]
    df = df[keep].dropna()
    print(f"  modeling rows: {len(df):,}", flush=True)

    # --- SANITY GATE: reproduce the published adjusted 304-CGC OR = 1.484 [1.470, 1.499] ---
    print("\n  [gate] reproduce published adjusted 304-CGC baseline:", flush=True)
    gate_row = fit_adjusted(df, "is_cancer", "gate_baseline_CGC304_all_variants")
    gate_ok = (abs(gate_row["OR"] - 1.484) < 0.002 and abs(gate_row["OR_lo"] - 1.470) < 0.002
               and abs(gate_row["OR_hi"] - 1.499) < 0.002)
    print(f"    -> adjusted baseline gate: got OR={gate_row['OR']} "
          f"[{gate_row['OR_lo']}, {gate_row['OR_hi']}], published 1.484 [1.470, 1.499] "
          f"({'MATCH' if gate_ok else 'MISMATCH'})")
    if not gate_ok:
        raise SystemExit("ADJUSTED BASELINE NOT REPRODUCED — stopping per protocol.")

    print("\n  germline-CGC adjusted models (vs all other genes):", flush=True)
    adj = {}
    adj["all_variants"] = fit_adjusted(df, "is_germline_cgc", "gcgc_all_variants")
    adj["exclude_somatic"] = fit_adjusted(df[~df["has_somatic"]], "is_germline_cgc",
                                          "gcgc_exclude_somatic")

    print("\n  germline-CGC adjusted models (vs TRUE non-cancer; all CGC/Census genes dropped):",
          flush=True)
    adj_tn = {}
    adj_tn["all_variants"] = fit_adjusted(
        df[df["is_cgc_not_gcgc"] == 0], "is_germline_cgc", "gcgc_vs_true_noncancer_all")
    adj_tn["exclude_somatic"] = fit_adjusted(
        df[(df["is_cgc_not_gcgc"] == 0) & ~df["has_somatic"]], "is_germline_cgc",
        "gcgc_vs_true_noncancer_excl_somatic")

    # ---------------- T3: side-by-side (pull 304 / 106 from prior run) ----------------
    print("\n=== T3: side-by-side table ===", flush=True)
    prior = json.load(open(CGC_SUMMARY))
    grid = prior["grid"]

    def find(gene_set, origin_filter):
        for r in grid:
            if r["gene_set"] == gene_set and r["origin_filter"] == origin_filter:
                return r
        return None

    def row_from_prior(label, gene_set, of):
        r = find(gene_set, of)
        return {
            "gene_set": label, "origin_filter": of, "n_exposed": r["n_exposed"],
            "crude_rate_pct": r["rate_exposed_pct"], "comparator_rate_pct": r["rate_comparator_pct"],
            "RR": r["RR"], "RR_95CI": r["RR_95CI"], "OR_crude": r["OR_crude"],
            "OR_crude_95CI": r["OR_crude_95CI"], "OR_adjusted": r["OR_adjusted"],
            "OR_adjusted_95CI": r["OR_adjusted_95CI"], "adjusted_model_n": r["adjusted_model_n"],
        }

    def row_from_new(of_key, of_label):
        cr = crude_res[of_key]; aj = adj[of_key]
        return {
            "gene_set": f"Germline-CGC ({len(germline_set)})", "origin_filter": of_label,
            "n_exposed": cr["n1"], "crude_rate_pct": cr["rate1_pct"],
            "comparator_rate_pct": cr["rate2_pct"], "RR": cr["risk_ratio"],
            "RR_95CI": cr["rr_95ci"], "OR_crude": cr["odds_ratio"],
            "OR_crude_95CI": cr["or_95ci"], "OR_adjusted": aj["OR"],
            "OR_adjusted_95CI": [aj["OR_lo"], aj["OR_hi"]], "adjusted_model_n": aj["n"],
        }

    table = [
        row_from_prior("COSMIC CGC (304) [baseline]", "COSMIC CGC (304)", "none (baseline)"),
        row_from_prior("COSMIC CGC (304) [baseline]", "COSMIC CGC (304)", "exclude has_somatic"),
        row_from_prior("Hereditary Cancer (106)", "Hereditary Cancer (106)", "none"),
        row_from_prior("Hereditary Cancer (106)", "Hereditary Cancer (106)", "exclude has_somatic"),
        row_from_new("all_variants", "none"),
        row_from_new("exclude_somatic", "exclude has_somatic"),
    ]

    # ---------------- T4 table: germline-CGC vs true non-cancer ----------------
    def row_tn(of_key, of_label):
        cr = tn_crude[of_key]; aj = adj_tn[of_key]
        return {
            "gene_set": f"Germline-CGC ({len(germline_set)})",
            "comparator": "True non-cancer (no CGC/Census gene)",
            "origin_filter": of_label,
            "n_exposed": cr["n1"], "n_comparator": cr["n2"],
            "crude_rate_pct": cr["rate1_pct"], "comparator_rate_pct": cr["rate2_pct"],
            "RR": cr["risk_ratio"], "RR_95CI": cr["rr_95ci"],
            "OR_crude": cr["odds_ratio"], "OR_crude_95CI": cr["or_95ci"],
            "OR_adjusted": aj["OR"], "OR_adjusted_95CI": [aj["OR_lo"], aj["OR_hi"]],
            "adjusted_model_n": aj["n"],
        }
    table_tn = [row_tn("all_variants", "none"),
                row_tn("exclude_somatic", "exclude has_somatic")]

    # ---------------- write outputs ----------------
    res = {
        "generated_for": "ClinVar revision — CGC's OWN germline annotation as the cancer exposure",
        "cgc_tsv": str(CGC_TSV),
        "cgc_tsv_version": "Cosmic_CancerGeneCensus_v101_GRCh37",
        "T1_germline_set_and_overlap": overlap,
        "T2_crude": crude_res,
        "T2_adjusted": adj,
        "sanity_gates": {
            "crude_published_304_rate1_pct": 16.619,
            "crude_published_304_rate2_pct": 8.191,
            "crude_published_304_RR": 2.029, "crude_published_304_OR": 2.234,
            "adjusted_published_304_OR": 1.484,
            "adjusted_published_304_OR_95CI": [1.470, 1.499],
            "adjusted_304_reproduced_OR": gate_row["OR"],
            "adjusted_304_reproduced_OR_95CI": [gate_row["OR_lo"], gate_row["OR_hi"]],
            "adjusted_baseline_match": bool(gate_ok),
        },
        "T3_side_by_side": table,
        "T4_vs_true_noncancer": {
            "definition": "exposure = germline-CGC (114, GERMLINE=='y'); comparator = genes "
                          "NOT in paper's 304 cancer_genes AND not in germline-CGC (all "
                          "Census/cancer genes dropped). Adjusted fit drops the "
                          "in-304-but-not-germline genes (is_cgc_not_gcgc==1); mirrors the "
                          "106-HC 'vs true non-cancer' contrast in 13_cgc_sensitivity.py.",
            "crude": tn_crude,
            "adjusted": adj_tn,
            "table": table_tn,
        },
        "notes": {
            "germline_set_definition": "CGC v101 genes with GERMLINE=='y' (COSMIC's own flag)",
            "comparator": "all genes NOT in the germline-CGC set (straight substitution, "
                          "matching the 'vs all others' contrast used for HC in 13_cgc_sensitivity.py)",
            "exposure_flag_symbols": "GENE_SYMBOL if present in dataset, else COSMIC SYNONYM present in dataset",
            "adjusted_model_formula": "05_multivariable_baseline.py:82-84 (M1), exposure swapped",
            "baseline_covariates_source": str(OUT / "baseline_covariates.csv"),
            "origin_source": str(ORIGIN_MAP),
        },
    }
    json.dump(res, open(OUT / "cosmic_germline_sensitivity.json", "w"), indent=2)
    print(f"\nwritten: {OUT/'cosmic_germline_sensitivity.json'}")

    # markdown
    def ci(x):
        return f"[{x[0]}, {x[1]}]" if x and x[0] is not None else "—"

    def adjci(r):
        return f"{r['OR_adjusted']} {ci(r['OR_adjusted_95CI'])}" if r["OR_adjusted"] is not None else "—"

    lines = []
    lines.append("# COSMIC-own-germline sensitivity — cancer vs non-cancer reclassification\n")
    lines.append(f"Germline-CGC set = {len(germline_set)} CGC v101 genes with GERMLINE=='y' "
                 f"(COSMIC's own annotation).\n")
    lines.append(f"- Matched to dataset: {len(matched_direct)} direct + {len(matched_via_syn)} via "
                 f"synonym = {len(matched_direct)+len(matched_via_syn)}/{len(germline_set)}; "
                 f"unmatched {len(unmatched)} {unmatched or ''}.")
    lines.append(f"- Overlap with paper's 304 CGC: {overlap['overlap_with_paper_304_cgc']}/"
                 f"{len(germline_set)}. Overlap with 106 Hereditary-Cancer: "
                 f"{overlap['overlap_with_106_hereditary_cancer']}/{len(germline_set)}.\n")
    lines.append("| Gene set | Origin filter | n (exposed) | Crude rate % | Comparator rate % | "
                 "RR (95% CI) | Crude OR (95% CI) | Adjusted OR (95% CI) |")
    lines.append("|---|---|---:|---:|---:|---|---|---|")
    for r in table:
        lines.append(
            f"| {r['gene_set']} | {r['origin_filter']} | {r['n_exposed']:,} | "
            f"{r['crude_rate_pct']} | {r['comparator_rate_pct']} | "
            f"{r['RR']} {ci(r['RR_95CI'])} | {r['OR_crude']} {ci(r['OR_crude_95CI'])} | "
            f"{adjci(r)} |")
    lines.append("")
    lines.append("## Germline-CGC (114) vs TRUE non-cancer (all CGC/Census genes dropped from comparator)\n")
    lines.append("| Gene set | Origin filter | n (exposed) | Crude rate % | Comparator rate % | "
                 "RR (95% CI) | Crude OR (95% CI) | Adjusted OR (95% CI) |")
    lines.append("|---|---|---:|---:|---:|---|---|---|")
    for r in table_tn:
        lines.append(
            f"| {r['gene_set']} vs true non-cancer | {r['origin_filter']} | {r['n_exposed']:,} | "
            f"{r['crude_rate_pct']} | {r['comparator_rate_pct']} | "
            f"{r['RR']} {ci(r['RR_95CI'])} | {r['OR_crude']} {ci(r['OR_crude_95CI'])} | "
            f"{adjci(r)} |")
    lines.append("")
    lines.append(f"Sanity gate: adjusted 304-CGC baseline reproduced OR="
                 f"{gate_row['OR']} [{gate_row['OR_lo']}, {gate_row['OR_hi']}] "
                 f"(published 1.484 [1.470, 1.499]) — "
                 f"{'MATCH' if gate_ok else 'MISMATCH'}.")
    md = "\n".join(lines) + "\n"
    open(OUT / "cosmic_germline_sensitivity.md", "w").write(md)
    print(f"written: {OUT/'cosmic_germline_sensitivity.md'}")
    print("\n" + md)


if __name__ == "__main__":
    main()
