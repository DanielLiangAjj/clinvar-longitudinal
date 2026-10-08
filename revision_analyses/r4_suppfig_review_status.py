#!/usr/bin/env python
"""R4: standalone Additional file 1 Fig. S4 (review status evolution) at full resolution,
legend outside the plotting area (the earlier version was a low-resolution crop of a
multi-panel figure with the legend over the data)."""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
REV = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses")
sys.path.insert(0, str(REV / "external_dependencies" / "ClinVar Research Project" / "scripts"))
import clinvar_visualizations as V
hist = json.load(open("/Users/danielliang/Downloads/ClinVar_Research_Project/figures/.cache/history_agg.json"))
dates = sorted(hist.keys())
x = np.array([V.date_to_float(d) for d in dates])
fig, ax = plt.subplots(figsize=(7.0, 3.6))
bottom = np.zeros_like(x)
for cat in V.REVIEW_ORDER:
    pct = np.array([hist[d]["review_statuses"].get(cat, 0) / max(hist[d]["total"], 1) * 100 for d in dates])
    ax.fill_between(x, bottom, bottom + pct, color=V.REVIEW_COLORS[cat], linewidth=0.15,
                    edgecolor="white", alpha=0.9)
    bottom = bottom + pct
ax.set_ylim(0, 100); ax.set_xlim(x.min(), x.max())
ax.set_ylabel("Variants (%)", fontsize=10); ax.set_xlabel("Year", fontsize=10)
ax.tick_params(labelsize=9)
ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
ax.set_title("Review status evolution", fontsize=11)
handles = [mpatches.Patch(color=V.REVIEW_COLORS[c], label=c) for c in V.REVIEW_ORDER]
ax.legend(handles=handles, fontsize=8.5, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
fig.tight_layout()
out = REV / "figures" / "r4" / "FigS4_review_status_evolution.png"
fig.savefig(out, dpi=600)
print("wrote", out)
