#!/usr/bin/env python
"""
Fix for submitted Fig 10 (inter-submitter concordance) — Reviewer 2 minor:
  "Avoid dual y-axes in Figure 10 ... two y-axes ... usually confusing."

Redraws it on a SINGLE percentage axis (concordant %, discordant %, and % of variants that are
multi-submitter), VCV era only (concordance vocabulary not meaningful pre-2020). Data from the
14 GB-history monthly pass.
"""
from pathlib import Path
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/outputs")
FIG = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses/figures")
ts = pd.read_csv(OUT / "monthly_timeseries.csv")
ts["ti"] = ts["month"].map(lambda s: int(s[:4]) + (int(s[5:7]) - 1) / 12)
v = ts[ts["ti"] >= 2020].copy()
v["disc_of_multi"] = 100 * v["n_conflicting_reviewstatus"] / v["n_multi"]
v["multi_of_all"] = 100 * v["n_multi"] / v["n_variants"]

fig, ax = plt.subplots(figsize=(8.4, 4.6))
ax.plot(v["ti"], v["pct_concordant_of_multi"], "-", color="#1b7837", lw=2, label="Concordant (% of multi-submitter)")
ax.plot(v["ti"], v["disc_of_multi"], "-", color="#b2182b", lw=2, label="Discordant / conflicting (% of multi-submitter)")
ax.plot(v["ti"], v["multi_of_all"], "--", color="#2166ac", lw=1.6, label="Multi-submitter (% of all variants)")
ax.set_xlabel("Year"); ax.set_ylabel("Percentage (%)")
ax.set_title("Fig 10 (revised): inter-submitter concordance, VCV era — single axis")
ax.set_ylim(0, 80); ax.grid(alpha=0.25); ax.legend(frameon=False, fontsize=9)
plt.tight_layout()
fig.savefig(FIG / "fig10_concordance_single_axis.png", dpi=300)
fig.savefig(FIG / "fig10_concordance_single_axis.pdf")
print("wrote figures/fig10_concordance_single_axis.png")
print(f"July-2025 check: concordant={v['pct_concordant_of_multi'].iloc[-1]:.1f}%  "
      f"discordant={v['disc_of_multi'].iloc[-1]:.1f}%  multi={v['multi_of_all'].iloc[-1]:.1f}%")
