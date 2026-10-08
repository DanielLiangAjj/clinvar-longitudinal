#!/usr/bin/env python
"""
Revision-3 analysis #31 — flow-based ITS outcomes + rebuilt Figure 10.

Addresses Reviewer 2, ROUND-2 Major comment #5:
  "The overall proportion of VUS is a cumulative database measure ... A more informative
   primary outcome would be the monthly proportion of newly appearing variants whose initial
   classification was VUS, complemented by monthly VUS-resolution rates among existing
   variants. The manuscript should also clarify that the current analysis demonstrates
   temporal associations rather than a causal effect ..."

Two new PRIMARY outcomes (both flow measures, insensitive to database growth/stock effects):
  (A) monthly share of NEWLY APPEARING variants whose initial normalized category is VUS
      (summary CSV First_Date x First_Classification, normalized with the pipeline's own
      categorize_classification)
  (B) monthly VUS RESOLUTION rate = reclassification events VUS -> P/LP or B/LB in month t
      divided by the VUS stock at month t-1 (events detail + monthly_timeseries n_VUS);
      VUS -> Conflicting departures are computed alongside for context.

Both get the same segmented (interrupted time-series) regression used in 06_its.py
(level+slope changes at 2015-06 ACMG/AMP and 2020-01 VCV; HAC/Newey-West SEs), and the
round-2 outcome (database-wide VUS proportion) is RETAINED as a descriptive panel.

Figure 10 (r3): a  VCV-era concordance (unchanged data)
                b  ITS of the new-variant initial-VUS share      (new primary)
                c  ITS of the monthly VUS resolution rate        (new primary)
                d  database-wide VUS proportion                  (retained, descriptive)

Outputs: outputs/its_r3.json, figures/r3/Figure_10_concordance_ITS_r3.(png|pdf)
"""
import csv, json, sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
REV  = DATA / "revision_analyses"
EXT  = REV / "external_dependencies" / "ClinVar Research Project"
OUT  = REV / "outputs"
FIGR3 = REV / "figures" / "r3"
FIGR3.mkdir(parents=True, exist_ok=True)
SUMMARY = DATA / "summary_20260219_024805.csv"
EVENTS  = Path("/Users/danielliang/Downloads/ClinVar_Research_Project/data/reclassification_events_detail.csv")

sys.path.insert(0, str(EXT / "scripts"))
from clinvar_visualizations import categorize_classification  # noqa: E402

P_SIDE = {"P", "LP", "P/LP"}
B_SIDE = {"B", "LB", "B/LB"}
B2015 = 2015 * 12 + 6
B2020 = 2020 * 12 + 1


def midx(s):
    y, m = s.split("-")[:2]
    return int(y) * 12 + int(m)


def segmented(t, y, breaks, label, hac_lags=6):
    """Same spec as 06_its.py: y ~ time + sum_k(level_k + slope_k), HAC(6)."""
    t = np.asarray(t, float); y = np.asarray(y, float)
    X = {"time": t - t[0]}
    for i, b in enumerate(breaks):
        X[f"level{i}"] = (t >= b).astype(float)
        X[f"slope{i}"] = np.clip(t - b, 0, None)
    Xd = sm.add_constant(pd.DataFrame(X))
    m = sm.OLS(y, Xd).fit(cov_type="HAC", cov_kwds={"maxlags": hac_lags})
    out = {"label": label, "n": int(len(y)), "r2": round(float(m.rsquared), 3), "terms": {}}
    for name in Xd.columns:
        out["terms"][name] = {"coef": float(m.params[name]), "p": float(m.pvalues[name]),
                              "ci": [float(m.conf_int().loc[name, 0]), float(m.conf_int().loc[name, 1])]}
    return out, m


def main():
    ts = pd.read_csv(OUT / "monthly_timeseries.csv")
    ts["ti"] = ts["month"].map(midx)
    ts = ts.sort_values("ti").reset_index(drop=True)

    # ---- (A) monthly share of new variants entering as VUS ---------------------
    print("Streaming summary for initial classifications ...", flush=True)
    new_tot, new_vus = Counter(), Counter()
    with open(SUMMARY) as fh:
        for row in csv.DictReader(fh):
            fd = row.get("First_Date", "")
            if not fd:
                continue
            new_tot[fd] += 1
            if categorize_classification(row.get("First_Classification", "")) == "VUS":
                new_vus[fd] += 1
    ts["n_new"] = ts["month"].map(lambda m: new_tot.get(m, 0))
    ts["n_new_vus"] = ts["month"].map(lambda m: new_vus.get(m, 0))
    ts["pct_new_vus"] = 100 * ts["n_new_vus"] / ts["n_new"].replace(0, np.nan)

    # ---- (B) monthly VUS resolution rate ---------------------------------------
    print("Counting VUS departures per month ...", flush=True)
    res_ev, conf_ev = Counter(), Counter()
    with open(EVENTS) as fh:
        for row in csv.DictReader(fh):
            if row["Prev_Category"] != "VUS":
                continue
            cur = row["Curr_Category"]
            if cur in P_SIDE or cur in B_SIDE:
                res_ev[row["Curr_Date"]] += 1
            elif cur == "Conflicting":
                conf_ev[row["Curr_Date"]] += 1
    ts["n_vus_resolved"] = ts["month"].map(lambda m: res_ev.get(m, 0))
    ts["n_vus_to_conf"] = ts["month"].map(lambda m: conf_ev.get(m, 0))
    ts["vus_stock_prev"] = ts["n_VUS"].shift(1)
    ts["res_rate_pct"] = 100 * ts["n_vus_resolved"] / ts["vus_stock_prev"]
    ts["conf_rate_pct"] = 100 * ts["n_vus_to_conf"] / ts["vus_stock_prev"]

    results = {"definitions": {
        "pct_new_vus": "share of variants first appearing in month t whose initial normalized category is VUS",
        "res_rate_pct": "events VUS->P/LP/B/LB in month t / VUS stock at month t-1 x 100",
        "conf_rate_pct": "events VUS->Conflicting in month t / VUS stock at month t-1 x 100"}}

    # descriptive means by era for the text
    def era_means(col):
        pre = ts.loc[ts["ti"] < B2015, col].mean()
        mid = ts.loc[(ts["ti"] >= B2015) & (ts["ti"] < B2020), col].mean()
        vcv = ts.loc[ts["ti"] >= B2020, col].mean()
        return {"pre2015": round(float(pre), 3), "2015_2019": round(float(mid), 3),
                "vcv_2020plus": round(float(vcv), 3)}
    results["era_means"] = {c: era_means(c) for c in ["pct_new_vus", "res_rate_pct", "conf_rate_pct"]}

    # ITS fits (drop months with missing outcome; note small-month exposure)
    a = ts.dropna(subset=["pct_new_vus"])
    rA, mA = segmented(a["ti"], a["pct_new_vus"], [B2015, B2020], "Share of new variants entering as VUS (%)")
    results["new_vus_its"] = rA
    b = ts.dropna(subset=["res_rate_pct"])
    rB, mB = segmented(b["ti"], b["res_rate_pct"], [B2015, B2020], "Monthly VUS resolution rate (%)")
    results["vus_resolution_its"] = rB
    # retained round-2 outcome, refit for the descriptive panel overlay
    rD, mD = segmented(ts["ti"], ts["pct_VUS"], [B2015, B2020], "Database-wide VUS proportion (%)")
    results["vus_stock_its_retained"] = rD

    # robustness: exclude months with < 500 new variants (early sparse archive)
    a2 = a[a["n_new"] >= 500]
    rA2, _ = segmented(a2["ti"], a2["pct_new_vus"], [B2015, B2020], "new-VUS share, months with >=500 new variants")
    results["new_vus_its_min500"] = rA2

    json.dump(results, open(OUT / "its_r3.json", "w"), indent=2)
    for tag in ["new_vus_its", "vus_resolution_its", "vus_stock_its_retained"]:
        r = results[tag]
        print(f"\n=== {r['label']} (R2={r['r2']}, n={r['n']}) ===")
        for name, v in r["terms"].items():
            star = "***" if v["p"] < 0.001 else ("**" if v["p"] < 0.01 else ("*" if v["p"] < 0.05 else ""))
            print(f"  {name:8s} coef={v['coef']:+.5f}  p={v['p']:.2e} {star}")
    print("\nera means:", json.dumps(results["era_means"], indent=1))

    # ---- Figure 10 (r3) --------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 1, figsize=(7.6, 11.6))
    fig.subplots_adjust(hspace=0.52)
    yrs = ts["ti"] / 12.0

    def vlines(ax):
        ax.axvline(2015 + 6/12, color="#b2182b", ls="--", lw=1.0)
        ax.axvline(2020 + 1/12, color="#2166ac", ls="--", lw=1.0)

    def fitline(ax, m, tvals, color):
        pred = m.predict(sm.add_constant(pd.DataFrame({
            "time": tvals - tvals.iloc[0],
            "level0": (tvals >= B2015).astype(float),
            "slope0": np.clip(tvals - B2015, 0, None),
            "level1": (tvals >= B2020).astype(float),
            "slope1": np.clip(tvals - B2020, 0, None)})))
        for lo, hi in [(tvals.min(), B2015 - 1), (B2015, B2020 - 1), (B2020, tvals.max())]:
            seg = (tvals >= lo) & (tvals <= hi)
            ax.plot(tvals[seg] / 12.0, np.asarray(pred)[seg], color=color, lw=1.8, zorder=4)

    # (a) VCV-era concordance — unchanged from round 2
    ax = axes[0]
    vcv = ts[ts["ti"] >= B2020]
    ax.plot(vcv["ti"]/12, vcv["pct_concordant_of_multi"], color="#1a7a3a", lw=1.6,
            label="Concordant (% of multi-submitter)")
    ax.plot(vcv["ti"]/12, 100 * vcv["n_conflicting_reviewstatus"] / vcv["n_multi"],
            color="#b2182b", lw=1.6, label="Discordant / conflicting (% of multi-submitter)")
    ax.plot(vcv["ti"]/12, 100 * vcv["n_multi"] / vcv["n_variants"], color="#2166ac",
            ls="--", lw=1.3, label="Multi-submitter (% of all variants)")
    ax.set_title("a  Inter-submitter concordance among multi-submitter variants (VCV era)",
                 loc="left", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("Percentage (%)"); ax.set_ylim(0, 80)
    ax.grid(alpha=0.25); ax.legend(fontsize=6.5, frameon=False, loc="center left")

    # (b) new-variant initial-VUS share — new primary
    ax = axes[1]
    ax.plot(a["ti"]/12, a["pct_new_vus"], ".", ms=3.2, color="#666666", alpha=0.7)
    fitline(ax, mA, a["ti"], "#7b3294"); vlines(ax)
    ax.set_title("b  Interrupted time-series: share of newly appearing variants entering as VUS",
                 loc="left", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("New variants entering\nas VUS (%)")
    ax.text(2015.55, ax.get_ylim()[1]*0.97, "ACMG/AMP 2015", color="#b2182b", fontsize=6.5, va="top")
    ax.text(2020.15, ax.get_ylim()[1]*0.97, "VCV 2020", color="#2166ac", fontsize=6.5, va="top")
    ax.grid(alpha=0.25)

    # (c) VUS resolution rate — new primary
    ax = axes[2]
    ax.plot(b["ti"]/12, b["res_rate_pct"], ".", ms=3.2, color="#666666", alpha=0.7)
    fitline(ax, mB, b["ti"], "#e08214"); vlines(ax)
    ax.set_title("c  Interrupted time-series: monthly VUS resolution rate (VUS → P/LP or B/LB)",
                 loc="left", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("VUS resolved per month\n(% of VUS stock)")
    ax.grid(alpha=0.25)

    # (d) database-wide VUS proportion — retained, descriptive
    ax = axes[3]
    ax.plot(yrs, ts["pct_VUS"], ".", ms=3.2, color="#888888", alpha=0.6)
    fitline(ax, mD, ts["ti"], "#2b3a67"); vlines(ax)
    ax.set_title("d  Database-wide VUS proportion (cumulative stock measure; descriptive)",
                 loc="left", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("VUS proportion (%)"); ax.set_xlabel("Year")
    ax.grid(alpha=0.25)

    for ax in axes:
        ax.set_xlim(2012.5, 2026.0)

    fig.savefig(FIGR3 / "Figure_10_concordance_ITS_r3.png", dpi=400, bbox_inches="tight")
    fig.savefig(FIGR3 / "Figure_10_concordance_ITS_r3.pdf", bbox_inches="tight")
    print("\nwrote", FIGR3 / "Figure_10_concordance_ITS_r3.png")


if __name__ == "__main__":
    main()
