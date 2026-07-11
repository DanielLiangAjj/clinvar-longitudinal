#!/usr/bin/env python
"""
Revision analysis #1 — Censored time-to-event (survival) analysis of variant reclassification.

Addresses Reviewer 2, Major point #2 (editor-flagged essential):
  "The time-to-reclassification analysis requires proper survival/time-to-event methods...
   Because variants enter ClinVar at different times and have unequal follow-up, this should
   be handled as a censored time-to-event analysis, not only among reclassified variants."

Approach (right-censored Kaplan-Meier on the VARIANT-AGE time scale):
  - Every variant is observed from its own first ClinVar appearance (t=0). No left truncation.
  - Entry date  = First_Date (summary CSV), month resolution.
  - Study end   = last snapshot = 2025-07.
  - Event (=1)  = first change in normalized category; event time = earliest Curr_Date in the
                  event-detail file minus First_Date.
  - Censor (=0) = variant never reclassified; censored at 2025-07 minus First_Date.
  - KM survivor S(t)=P(not yet reclassified by age t); cumulative incidence = 1-S(t).
  - Greenwood variance; log-log 95% CI (bounded in [0,1]).

Outputs:
  outputs/survival_overall.json          headline numbers (overall + by stratum)
  outputs/km_curve_overall.csv           full KM curve for plotting
  figures/survival_km.png/.pdf           KM cumulative-incidence figure
"""
import json, csv, sys
from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
FIG  = DATA / "revision_analyses" / "figures"
SUMMARY = DATA / "summary_20260219_024805.csv"
EVENTS  = PROJ / "reclassification_events_detail.csv"
STUDY_END = 2025 * 12 + 7          # 2025-07 as month index
DOMAINS = json.load(open(PROJ / "disease_domains.json"))

def month_idx(s):
    """'2012-11' -> 24155 (year*12+month). Returns NaN on bad input."""
    try:
        y, m = s.split("-")[:2]
        return int(y) * 12 + int(m)
    except Exception:
        return np.nan

def km_estimate(times, events):
    """Right-censored Kaplan-Meier. times in months (int), events 0/1.
    Returns dict with per-event-time arrays and Greenwood SE of S(t)."""
    t = np.asarray(times, float); e = np.asarray(events, int)
    order = np.argsort(t, kind="mergesort"); t, e = t[order], e[order]
    n = len(t)
    uniq = np.unique(t[e == 1])                      # distinct event times
    S = 1.0; cum_var_term = 0.0
    rows = []
    for ti in uniq:
        n_i = np.sum(t >= ti)                         # at risk
        d_i = np.sum((t == ti) & (e == 1))            # events at ti
        if n_i == 0:
            continue
        S *= (1.0 - d_i / n_i)
        cum_var_term += d_i / (n_i * (n_i - d_i)) if n_i > d_i else 0.0
        se_S = S * np.sqrt(cum_var_term)              # Greenwood
        rows.append((ti, n_i, d_i, S, se_S))
    return rows, n

def loglog_ci(S, seS, z=1.96):
    """log-log transform CI for a survival probability (bounded [0,1])."""
    S = np.clip(S, 1e-12, 1 - 1e-12)
    theta = np.log(-np.log(S))
    se_theta = seS / (S * np.abs(np.log(S)))
    lo = np.exp(-np.exp(theta + z * se_theta))
    hi = np.exp(-np.exp(theta - z * se_theta))
    return lo, hi

def summarize(rows, n, label):
    """From KM rows, extract median (if reached) and cumulative incidence at fixed ages."""
    if not rows:
        return {"label": label, "n": int(n), "n_events": 0}
    arr = np.array([(r[0], r[3], r[4]) for r in rows])   # ti, S, seS
    ti, S, seS = arr[:, 0], arr[:, 1], arr[:, 2]
    n_events = int(sum(r[2] for r in rows))
    def ci_at(age_years):
        m = age_years * 12.0
        idx = np.searchsorted(ti, m, side="right") - 1
        if idx < 0:
            return 0.0, 0.0, 0.0
        s = S[idx]; lo, hi = loglog_ci(s, seS[idx])
        # cumulative incidence = 1 - S; CI flips
        return float(1 - s), float(1 - hi), float(1 - lo)
    # median time-to-event: smallest age where S <= 0.5
    below = np.where(S <= 0.5)[0]
    median_years = float(ti[below[0]] / 12.0) if len(below) else None
    ci = {f"cuminc_{y}yr": dict(zip(["est", "lo", "hi"], ci_at(y))) for y in (1, 2, 5, 10)}
    return {"label": label, "n": int(n), "n_events": n_events,
            "median_time_to_reclass_years": median_years,
            "median_reached": median_years is not None,
            **ci}

def main():
    print("Loading first-event dates from event-detail file ...", flush=True)
    # earliest Curr_Date (the month the new category first appears) per variant
    first_ev = {}
    with open(EVENTS) as fh:
        r = csv.DictReader(fh)
        for row in r:
            vid = row["Variation_ID"]; mi = month_idx(row["Curr_Date"])
            if vid not in first_ev or mi < first_ev[vid]:
                first_ev[vid] = mi
    print(f"  variants with >=1 event: {len(first_ev):,}", flush=True)

    print("Loading summary (3.68M variants) ...", flush=True)
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "Gene", "First_Date",
                                       "Classification_Changed", "Num_Submissions"],
                     dtype={"Variation_ID": str, "Gene": str, "First_Date": str})
    print(f"  rows: {len(df):,}", flush=True)

    df["entry"] = df["First_Date"].map(month_idx)
    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true")
    df["ev_month"] = df["Variation_ID"].map(first_ev)

    # event time (changed): first event month - entry ; censor time (unchanged): STUDY_END - entry
    df["event"] = df["changed"].astype(int)
    df.loc[df["changed"], "time_m"] = df.loc[df["changed"], "ev_month"] - df.loc[df["changed"], "entry"]
    df.loc[~df["changed"], "time_m"] = STUDY_END - df.loc[~df["changed"], "entry"]

    # data hygiene
    before = len(df)
    df = df.dropna(subset=["entry", "time_m"])
    df = df[df["time_m"] >= 0]
    # a changed variant with no event month (shouldn't happen) -> drop
    df = df[~(df["changed"] & df["ev_month"].isna())]
    print(f"  usable after hygiene: {len(df):,} (dropped {before-len(df):,})", flush=True)
    print(f"  events: {int(df['event'].sum()):,}  censored: {int((1-df['event']).sum()):,}", flush=True)

    results = {}
    # ---- overall
    rows, n = km_estimate(df["time_m"].values, df["event"].values)
    results["overall"] = summarize(rows, n, "All variants")
    # save full curve
    with open(OUT / "km_curve_overall.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["age_months", "age_years", "n_at_risk", "n_events", "S", "cuminc", "cuminc_lo", "cuminc_hi"])
        for ti, n_i, d_i, S, seS in rows:
            lo, hi = loglog_ci(S, seS)
            w.writerow([int(ti), round(ti/12, 4), int(n_i), int(d_i), round(S, 6),
                        round(1-S, 6), round(1-hi, 6), round(1-lo, 6)])

    # ---- cancer vs non-cancer
    cancer_genes = set(DOMAINS["cancer_genes"])
    df["cancer"] = df["Gene"].isin(cancer_genes)
    for lab, mask in [("Cancer genes (COSMIC CGC)", df["cancer"]),
                      ("Non-cancer genes", ~df["cancer"])]:
        rows, n = km_estimate(df.loc[mask, "time_m"].values, df.loc[mask, "event"].values)
        results[f"cancer::{lab}"] = summarize(rows, n, lab)

    # ---- by disease domain
    gene2dom = {}
    for dom, genes in DOMAINS["disease_domains"].items():
        for g in genes:
            gene2dom[g] = dom
    df["domain"] = df["Gene"].map(gene2dom).fillna("Other/Unclassified")
    dom_curves = {}
    for dom in list(DOMAINS["disease_domains"].keys()):
        mask = df["domain"] == dom
        rows, n = km_estimate(df.loc[mask, "time_m"].values, df.loc[mask, "event"].values)
        results[f"domain::{dom}"] = summarize(rows, n, dom)
        dom_curves[dom] = [(ti/12.0, 1-S) for ti, _, _, S, _ in rows]

    json.dump(results, open(OUT / "survival_overall.json", "w"), indent=2)

    # ---- figure: cumulative incidence, overall + by domain
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    # panel a: overall
    cur = pd.read_csv(OUT / "km_curve_overall.csv")
    ax[0].step(cur["age_years"], cur["cuminc"]*100, where="post", color="#2b3a67", lw=1.8)
    ax[0].fill_between(cur["age_years"], cur["cuminc_lo"]*100, cur["cuminc_hi"]*100,
                       step="post", alpha=0.18, color="#2b3a67")
    ax[0].set_xlabel("Years since first ClinVar appearance")
    ax[0].set_ylabel("Cumulative incidence of first\nreclassification (%)")
    ax[0].set_title("a  Overall (censored KM, n=3.68M)")
    ax[0].set_xlim(0, 13); ax[0].grid(alpha=0.25)
    # panel b: by domain
    palette = ["#b2182b", "#ef8a62", "#2166ac", "#67a9cf", "#1b7837"]
    for (dom, curve), c in zip(dom_curves.items(), palette):
        if curve:
            xs, ys = zip(*curve)
            ax[1].step(xs, [y*100 for y in ys], where="post", label=dom, color=c, lw=1.6)
    ax[1].set_xlabel("Years since first ClinVar appearance")
    ax[1].set_ylabel("Cumulative incidence (%)")
    ax[1].set_title("b  By disease domain")
    ax[1].set_xlim(0, 13); ax[1].grid(alpha=0.25); ax[1].legend(fontsize=8, frameon=False)
    plt.tight_layout()
    fig.savefig(FIG / "survival_km.png", dpi=300)
    fig.savefig(FIG / "survival_km.pdf")
    print("\n=== DONE. Key numbers: ===")
    o = results["overall"]
    print(f"Overall: n={o['n']:,}, events={o['n_events']:,}, median reached={o['median_reached']}")
    for y in (1, 2, 5, 10):
        c = o[f"cuminc_{y}yr"]
        print(f"  cumulative incidence @ {y:>2}yr: {c['est']*100:5.2f}%  (95% CI {c['lo']*100:.2f}-{c['hi']*100:.2f})")

if __name__ == "__main__":
    main()
