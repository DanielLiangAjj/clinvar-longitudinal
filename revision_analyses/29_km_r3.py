#!/usr/bin/env python
"""
Revision-3 analysis #29 — competing-risks / sensitivity extensions of the KM analysis.

Addresses Reviewer 2, ROUND-2 Major comment #2:
  (i)   The naive 2.58-yr median among reclassified variants (and the per-transition naive
        medians 1.92/1.92/3.17/0.67 yr in the same paragraph and in Supplementary Figure 1)
        is biased by censoring: recently entered variants can only contribute SHORT observed
        times. All of these are REMOVED from the manuscript. In their place, this script
        estimates CAUSE-SPECIFIC cumulative incidence of the FIRST reclassification by
        transition type (Aalen-Johansen estimator, other first-events treated as competing),
        which is the proper analogue and becomes the new Supplementary Figure 1.
  (ii)  The Methods must define event time, censoring, first-snapshot handling, and
        disappeared/merged/withdrawn variants. This script quantifies the last two:
          - sensitivity A: exclude the variants already present in the first archived
            snapshot (2012-11), whose pre-archive history is unobservable;
          - sensitivity B: exclude the variants absent from the final (2025-07) release
            (merged/withdrawn/removed), whose follow-up is overstated by censoring at
            study end. IDs come from the 2025-07 checkpoint.
  (iii) Conditional statement check: F(5)/F(10) — "of variants estimated to reclassify
        within ten years, approximately half do so within five".

Event/censoring definitions are identical to 01_survival_analysis.py (verified against
outputs/survival_overall.json before anything else runs).

Outputs:
  outputs/km_r3.json                                  all numbers for text + letter
  figures/r3/SuppFig1_transition_cif.(png|pdf)        new Supplementary Figure 1
"""
import csv, gzip, json, pickle
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
REV  = DATA / "revision_analyses"
OUT  = REV / "outputs"
FIGR3 = REV / "figures" / "r3"
FIGR3.mkdir(parents=True, exist_ok=True)
SUMMARY = DATA / "summary_20260219_024805.csv"
EVENTS  = Path("/Users/danielliang/Downloads/ClinVar_Research_Project/data/reclassification_events_detail.csv")
CKPT    = DATA / "ClinVarVariationRelease_2025-07.checkpoint.pkl.gz"
STUDY_END = 2025 * 12 + 7
FIRST_SNAPSHOT = "2012-11"
RNG = np.random.default_rng(20260821)
N_BOOT = 200

P_SIDE = {"P", "LP", "P/LP"}
B_SIDE = {"B", "LB", "B/LB"}
CAUSES = ["VUS->P/LP", "VUS->B/LB", "VUS->Conflicting",
          "P/LP->B/LB (drastic)", "B/LB->P/LP (drastic)", "Other first transition"]


def month_idx(s):
    try:
        y, m = s.split("-")[:2]
        return int(y) * 12 + int(m)
    except Exception:
        return np.nan


def cause_of(prev, curr):
    if prev == "VUS" and curr in P_SIDE: return "VUS->P/LP"
    if prev == "VUS" and curr in B_SIDE: return "VUS->B/LB"
    if prev == "VUS" and curr == "Conflicting": return "VUS->Conflicting"
    if prev in P_SIDE and curr in B_SIDE: return "P/LP->B/LB (drastic)"
    if prev in B_SIDE and curr in P_SIDE: return "B/LB->P/LP (drastic)"
    return "Other first transition"


def km_fixed_ages(times, events, ages=(1, 2, 5, 10)):
    """Overall KM cumulative incidence at fixed ages (per 01_; Greenwood + log-log CI)."""
    t = np.asarray(times, float); e = np.asarray(events, int)
    order = np.argsort(t, kind="mergesort"); t, e = t[order], e[order]
    te = t[e == 1]
    uniq = np.unique(te)
    n = len(t)
    S = 1.0; cum = 0.0
    ti_l, S_l, se_l = [], [], []
    # counts per unique event time via searchsorted on the sorted arrays
    for ti in uniq:
        n_i = n - np.searchsorted(t, ti, side="left")
        d_i = np.searchsorted(te, ti, side="right") - np.searchsorted(te, ti, side="left")
        if n_i <= 0: continue
        S *= (1 - d_i / n_i)
        cum += d_i / (n_i * (n_i - d_i)) if n_i > d_i else 0.0
        ti_l.append(ti); S_l.append(S); se_l.append(S * np.sqrt(cum))
    ti_a, S_a, se_a = np.array(ti_l), np.array(S_l), np.array(se_l)
    res = {}
    for y in ages:
        idx = np.searchsorted(ti_a, y * 12, side="right") - 1
        if idx < 0:
            res[f"cuminc_{y}yr"] = {"est": 0.0, "lo": 0.0, "hi": 0.0}
            continue
        s, se = S_a[idx], se_a[idx]
        sc = np.clip(s, 1e-12, 1 - 1e-12)
        th = np.log(-np.log(sc)); se_th = se / (sc * abs(np.log(sc)))
        lo_s = np.exp(-np.exp(th + 1.96 * se_th)); hi_s = np.exp(-np.exp(th - 1.96 * se_th))
        res[f"cuminc_{y}yr"] = {"est": float(1 - s), "lo": float(1 - hi_s), "hi": float(1 - lo_s)}
    return res, int(e.sum()), n


def aj_from_cells(cells, months):
    """Aalen-Johansen CIFs from aggregated counts.
    cells: dict month -> np.array([d_cause1..d_cause6, n_censored_at_month])
    Returns dict cause -> array of CIF over the sorted month grid, plus the grid."""
    grid = np.array(sorted(months))
    K = len(CAUSES)
    d = np.zeros((len(grid), K)); c = np.zeros(len(grid))
    for i, m in enumerate(grid):
        v = cells.get(m)
        if v is not None:
            d[i, :] = v[:K]; c[i] = v[K]
    total = d.sum() + c.sum()
    # at-risk just before each time: total minus everything that exited earlier
    exited_before = np.concatenate([[0.0], np.cumsum(d.sum(axis=1) + c)[:-1]])
    n_at = total - exited_before
    with np.errstate(divide="ignore", invalid="ignore"):
        h_any = np.where(n_at > 0, d.sum(axis=1) / n_at, 0.0)
    S_left = np.concatenate([[1.0], np.cumprod(1 - h_any)[:-1]])   # S(t-)
    cif = {}
    for k, cause in enumerate(CAUSES):
        with np.errstate(divide="ignore", invalid="ignore"):
            h_k = np.where(n_at > 0, d[:, k] / n_at, 0.0)
        cif[cause] = np.cumsum(S_left * h_k)
    return grid, cif


def main():
    print("Loading first events with causes ...", flush=True)
    first_ev = {}   # vid -> (month, cause)
    with open(EVENTS) as fh:
        for row in csv.DictReader(fh):
            vid = row["Variation_ID"]; mi = month_idx(row["Curr_Date"])
            if vid not in first_ev or mi < first_ev[vid][0]:
                first_ev[vid] = (mi, cause_of(row["Prev_Category"], row["Curr_Category"]))
    print(f"  variants with >=1 event: {len(first_ev):,}")

    print("Loading summary ...", flush=True)
    df = pd.read_csv(SUMMARY, usecols=["Variation_ID", "First_Date", "Classification_Changed"],
                     dtype={"Variation_ID": str, "First_Date": str})
    df["entry"] = df["First_Date"].map(month_idx)
    df["changed"] = df["Classification_Changed"].astype(str).str.lower().eq("true")
    ev = df["Variation_ID"].map(lambda v: first_ev.get(v))
    df["ev_month"] = ev.map(lambda x: x[0] if x else np.nan)
    df["cause"] = ev.map(lambda x: x[1] if x else None)
    df["event"] = df["changed"].astype(int)
    df["time_m"] = np.where(df["changed"], df["ev_month"] - df["entry"], STUDY_END - df["entry"])
    df = df.dropna(subset=["entry", "time_m"])
    df = df[df["time_m"] >= 0]
    df = df[~(df["changed"] & df["ev_month"].isna())]
    print(f"  usable: {len(df):,}  events: {int(df['event'].sum()):,}")

    res = {"n": int(len(df)), "n_events": int(df["event"].sum())}

    # ---- 0. verification against the round-2 headline numbers ------------------
    main_km, ne, n = km_fixed_ages(df["time_m"].values, df["event"].values)
    res["overall"] = main_km
    print("overall cuminc (must match survival_overall.json 2.544/15.528/30.85):")
    for y in (1, 2, 5, 10):
        c = main_km[f"cuminc_{y}yr"]
        print(f"  {y:>2}yr {c['est']*100:6.3f}  [{c['lo']*100:.3f}-{c['hi']*100:.3f}]")
    res["conditional_F5_over_F10"] = round(main_km["cuminc_5yr"]["est"] / main_km["cuminc_10yr"]["est"], 4)

    # ---- 1. cause-specific CIF (Aalen-Johansen) --------------------------------
    print("Aalen-Johansen cause-specific CIF ...", flush=True)
    cells = {}
    ev_rows = df[df["event"] == 1]
    for m, cause in zip(ev_rows["time_m"].astype(int), ev_rows["cause"]):
        v = cells.setdefault(m, np.zeros(len(CAUSES) + 1))
        v[CAUSES.index(cause)] += 1
    for m, cnt in df.loc[df["event"] == 0, "time_m"].astype(int).value_counts().items():
        v = cells.setdefault(m, np.zeros(len(CAUSES) + 1))
        v[len(CAUSES)] += cnt
    grid, cif = aj_from_cells(cells, cells.keys())
    res["first_event_cause_counts"] = {c: int(ev_rows["cause"].eq(c).sum()) for c in CAUSES}

    def at_age(arr, yr):
        idx = np.searchsorted(grid, yr * 12, side="right") - 1
        return float(arr[idx]) if idx >= 0 else 0.0

    res["cif_at_horizons"] = {c: {f"{y}yr": round(100 * at_age(cif[c], y), 3) for y in (1, 2, 5, 10)}
                              for c in CAUSES}

    # multinomial bootstrap over aggregated cells for horizon CIs
    print(f"bootstrap ({N_BOOT} reps) ...", flush=True)
    months_sorted = sorted(cells.keys())
    flat = np.array([cells[m] for m in months_sorted])          # (M, K+1)
    p = flat.flatten() / flat.sum()
    N = int(flat.sum())
    boot = {c: {y: [] for y in (5, 10)} for c in CAUSES}
    for _ in range(N_BOOT):
        draw = RNG.multinomial(N, p).reshape(flat.shape)
        cells_b = {m: draw[i] for i, m in enumerate(months_sorted)}
        g_b, cif_b = aj_from_cells(cells_b, months_sorted)
        for c in CAUSES:
            for y in (5, 10):
                idx = np.searchsorted(g_b, y * 12, side="right") - 1
                boot[c][y].append(100 * cif_b[c][idx])
    res["cif_ci"] = {c: {f"{y}yr": [round(float(np.percentile(boot[c][y], 2.5)), 3),
                                    round(float(np.percentile(boot[c][y], 97.5)), 3)]
                         for y in (5, 10)} for c in CAUSES}

    # ---- 2. sensitivity A: exclude first-snapshot variants ---------------------
    maskA = df["First_Date"] != FIRST_SNAPSHOT
    kmA, neA, nA = km_fixed_ages(df.loc[maskA, "time_m"].values, df.loc[maskA, "event"].values)
    res["sensitivity_excl_first_snapshot"] = {"n_excluded": int((~maskA).sum()), "n": nA,
                                              "n_events": neA, **kmA}

    # ---- 3. sensitivity B: exclude variants absent from the final release ------
    try:
        print("Loading 2025-07 checkpoint for final-release membership ...", flush=True)
        with gzip.open(CKPT, "rb") as fh:
            obj = pickle.load(fh)
        # checkpoint layout: {"snapshot_date": "2025-07", ..., "num_variants": 3594583,
        #                     "data": {variation_id: {...}, ...}}
        if isinstance(obj, dict) and isinstance(obj.get("data"), dict):
            final_ids = set(map(str, obj["data"].keys()))
        elif isinstance(obj, dict) and isinstance(obj.get("variants"), dict):
            final_ids = set(map(str, obj["variants"].keys()))
        else:
            raise TypeError(f"unexpected checkpoint layout: {type(obj)} keys={list(obj)[:6] if isinstance(obj, dict) else None}")
        del obj
        maskB = df["Variation_ID"].isin(final_ids)
        kmB, neB, nB = km_fixed_ages(df.loc[maskB, "time_m"].values, df.loc[maskB, "event"].values)
        res["sensitivity_excl_absent_from_final"] = {
            "n_final_release": len(final_ids),
            "n_absent": int((~maskB).sum()), "n": nB, "n_events": neB, **kmB}
        print(f"  final release holds {len(final_ids):,}; absent from it: {int((~maskB).sum()):,}")
    except Exception as e:
        print("  checkpoint unavailable:", e)
        res["sensitivity_excl_absent_from_final"] = {"error": str(e)}

    json.dump(res, open(OUT / "km_r3.json", "w"), indent=2)

    # ---- 4. new Supplementary Figure 1 -----------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COLORS = {"VUS->P/LP": "#b2182b", "VUS->B/LB": "#2166ac", "VUS->Conflicting": "#7b3294",
              "P/LP->B/LB (drastic)": "#e08214", "B/LB->P/LP (drastic)": "#1b7837",
              "Other first transition": "#999999"}
    PRETTY = {"VUS->P/LP": "VUS → P/LP", "VUS->B/LB": "VUS → B/LB",
              "VUS->Conflicting": "VUS → Conflicting",
              "P/LP->B/LB (drastic)": "P/LP → B/LB (drastic)",
              "B/LB->P/LP (drastic)": "B/LB → P/LP (drastic)",
              "Other first transition": "Other first transition"}
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 4.0))
    yrs = grid / 12.0
    for c in CAUSES:
        ax[0].step(yrs, 100 * cif[c], where="post", lw=1.6, color=COLORS[c], label=PRETTY[c])
    ax[0].set_xlim(0, 13); ax[0].set_xlabel("Years since first ClinVar appearance")
    ax[0].set_ylabel("Cause-specific cumulative incidence of\nfirst reclassification (%)")
    ax[0].set_title("a  First reclassification by transition type\n(Aalen–Johansen, competing risks)",
                    loc="left", fontsize=9)
    ax[0].grid(alpha=0.25); ax[0].legend(fontsize=6.5, frameon=False, loc="upper left")
    x = np.arange(len(CAUSES)); bw = 0.36
    v5 = [res["cif_at_horizons"][c]["5yr"] for c in CAUSES]
    v10 = [res["cif_at_horizons"][c]["10yr"] for c in CAUSES]
    e5 = np.array([[v5[i] - res["cif_ci"][c]["5yr"][0], res["cif_ci"][c]["5yr"][1] - v5[i]]
                   for i, c in enumerate(CAUSES)]).T
    e10 = np.array([[v10[i] - res["cif_ci"][c]["10yr"][0], res["cif_ci"][c]["10yr"][1] - v10[i]]
                    for i, c in enumerate(CAUSES)]).T
    ax[1].bar(x - bw/2, v5, bw, yerr=np.abs(e5), capsize=2, color="#9ecae1",
              edgecolor="white", label="5 years")
    ax[1].bar(x + bw/2, v10, bw, yerr=np.abs(e10), capsize=2, color="#2166ac",
              edgecolor="white", label="10 years")
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([PRETTY[c].replace(" (drastic)", "\n(drastic)").replace("Other first ", "Other first\n")
                           for c in CAUSES], fontsize=6.2)
    ax[1].set_ylabel("Cumulative incidence (%)")
    ax[1].set_title("b  Fixed-horizon estimates (95% bootstrap CI)", loc="left", fontsize=9)
    ax[1].grid(axis="y", alpha=0.25); ax[1].legend(fontsize=7, frameon=False)
    plt.tight_layout()
    fig.savefig(FIGR3 / "SuppFig1_transition_cif.png", dpi=400)
    fig.savefig(FIGR3 / "SuppFig1_transition_cif.pdf")
    print("\nwrote", FIGR3 / "SuppFig1_transition_cif.png")
    print(json.dumps(res["cif_at_horizons"], indent=2))


if __name__ == "__main__":
    main()
