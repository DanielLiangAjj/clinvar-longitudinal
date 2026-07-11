#!/usr/bin/env python
"""
Revision analysis #6 — Interrupted time-series (segmented regression).

Addresses Reviewer 2, Major point #5:
  "...since the study spans pre- and post-ACMG/AMP guideline periods, the authors should
   evaluate whether ... reclassification rates changed after 2015. An interrupted time-series
   analysis would substantially strengthen the manuscript."

We fit segmented regressions of two monthly outcomes available across the FULL period:
  (1) monthly reclassification rate = events that month / variants present that month
  (2) monthly VUS proportion
with level+slope change terms at two interventions:
  - 2015-06  ACMG/AMP guidelines (Richards et al., May 2015)
  - 2020-01  RCV->VCV schema transition (modeled to separate the schema break from real trends)
HAC (Newey-West) standard errors account for autocorrelation.

NOTE on concordance/conflict metrics: ClinVar's RCV-era review vocabulary did not distinguish
concordant vs conflicting multi-submitter variants, so those series are only interpretable in
the VCV era (2020+). We therefore restrict a concordance ITS to the VCV era where meaningful.

Outputs: outputs/its_results.json, figures/its.png/.pdf
"""
import csv, json, collections
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

DATA = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine")
PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project/data")
OUT  = DATA / "revision_analyses" / "outputs"
FIG  = DATA / "revision_analyses" / "figures"

def midx(s):
    y, m = s.split("-")[:2]; return int(y)*12 + int(m)

def segmented(t, y, breaks, label, hac_lags=6):
    """y_t ~ time + sum_k [ level_k * post_k + slope_k * (t - break_k)_+ ]."""
    t = np.asarray(t, float); y = np.asarray(y, float)
    X = {"time": t - t[0]}
    for i, b in enumerate(breaks):
        bt = b - t[0]
        X[f"level{i}"] = (t >= b).astype(float)
        X[f"slope{i}"] = np.clip(t - b, 0, None)
    Xd = pd.DataFrame(X); Xd = sm.add_constant(Xd)
    m = sm.OLS(y, Xd).fit(cov_type="HAC", cov_kwds={"maxlags": hac_lags})
    out = {"label": label, "n": int(len(y)), "r2": round(float(m.rsquared), 3), "terms": {}}
    for name in Xd.columns:
        out["terms"][name] = {"coef": float(m.params[name]),
                              "p": float(m.pvalues[name]),
                              "ci": [float(m.conf_int().loc[name, 0]), float(m.conf_int().loc[name, 1])]}
    return out, m

def main():
    ts = pd.read_csv(OUT / "monthly_timeseries.csv")
    ts["ti"] = ts["month"].map(midx)
    ts = ts.sort_values("ti").reset_index(drop=True)

    # monthly reclassification events (from event-detail file)
    ev = collections.Counter()
    with open(PROJ / "reclassification_events_detail.csv") as fh:
        for row in csv.DictReader(fh):
            ev[row["Curr_Date"]] += 1
    ts["n_events"] = ts["month"].map(lambda m: ev.get(m, 0))
    ts["reclass_rate_pct"] = 100 * ts["n_events"] / ts["n_variants"]

    B2015 = midx("2015-06"); B2020 = midx("2020-01")
    results = {}

    # (1) reclassification rate — full period, two breaks
    r1, m1 = segmented(ts["ti"], ts["reclass_rate_pct"], [B2015, B2020], "Monthly reclassification rate (%)")
    results["reclass_rate"] = r1
    # (2) VUS proportion — full period, two breaks
    r2, m2 = segmented(ts["ti"], ts["pct_VUS"], [B2015, B2020], "Monthly VUS proportion (%)")
    results["vus_pct"] = r2
    # (3) concordance among multi-submitter variants — VCV era only, one trend
    vcv = ts[ts["ti"] >= B2020].copy()
    r3, m3 = segmented(vcv["ti"], vcv["pct_concordant_of_multi"], [], "VCV-era concordance among multi-submitter (%)")
    results["concordance_vcv"] = r3

    json.dump(results, open(OUT / "its_results.json", "w"), indent=2)

    def describe(tag, r):
        print(f"\n=== {tag} (R2={r['r2']}, n={r['n']}) ===")
        for name, v in r["terms"].items():
            star = "***" if v["p"] < 0.001 else ("**" if v["p"] < 0.01 else ("*" if v["p"] < 0.05 else ""))
            print(f"  {name:8s} coef={v['coef']:+.5f}  p={v['p']:.2e} {star}")
    describe("Reclassification rate: level0/slope0=2015 break, level1/slope1=2020 break", results["reclass_rate"])
    describe("VUS proportion", results["vus_pct"])
    describe("VCV-era concordance trend", results["concordance_vcv"])

    # figure
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    yrs = ts["ti"]/12.0
    for a, col, ttl, mdl in [(ax[0], "reclass_rate_pct", "a  Monthly reclassification rate", m1),
                             (ax[1], "pct_VUS", "b  VUS proportion", m2)]:
        a.plot(yrs, ts[col], ".", ms=3, color="#555", alpha=0.6)
        a.axvline(2015+5/12, color="#b2182b", ls="--", lw=1, label="ACMG/AMP 2015")
        a.axvline(2020+0/12, color="#2166ac", ls="--", lw=1, label="VCV 2020")
        a.set_title(ttl); a.set_xlabel("Year"); a.grid(alpha=0.2); a.legend(fontsize=7, frameon=False)
    vyrs = vcv["ti"]/12.0
    ax[2].plot(vyrs, vcv["pct_concordant_of_multi"], ".", ms=4, color="#1b7837")
    ax[2].set_title("c  VCV-era concordance (multi-submitter)"); ax[2].set_xlabel("Year"); ax[2].grid(alpha=0.2)
    plt.tight_layout(); fig.savefig(FIG/"its.png", dpi=300); fig.savefig(FIG/"its.pdf")
    print("\nfigure written: figures/its.png")

if __name__ == "__main__":
    main()
