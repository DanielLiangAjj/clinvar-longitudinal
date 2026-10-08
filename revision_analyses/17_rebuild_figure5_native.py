#!/usr/bin/env python
"""
Rebuild manuscript Figure 5 (a-f) as ONE native matplotlib figure — replacing the PIL image-composite
(compose_revised_figures.py) that pasted a pre-baked a-d block next to shrunken e/f panels.

Fixes Kai Wang's review comment on Figure 5:
  "This figure needs to be re-made. Panel c has overlapping words. The 'e' and 'f' are huge compared
   to others, yet the font in e and f are tiny to read."

Root cause of both problems: it was never a single figure. Panels a-d were a 7.08x6.0in block rendered
with 4-5.5pt fonts; e/f were separate PNGs with 150px letters pasted on and their content scaled down
to fit two-across. Panel c's title/label/annotation overlap lived inside the a-d block.

This script renders all six panels in one 3x2 GridSpec with a single, uniform typographic scale and
one panel_label() helper, so every panel letter is the same size and every interior font is legible.
Panel c is re-laid-out so the title, the "c" label and the oscillating callout no longer collide.

Data (no re-run of the 3.68M-variant pass; all from saved outputs):
  a  summary_agg.transition_matrix            (figures/.cache/summary_agg.json)
  b  statistics.change_type_distribution      (data/statistics_20260219_024805.json)
  c  summary_agg stable/changed/oscillating   (figures/.cache/summary_agg.json)
  d  summary_agg.p_to_b / b_to_p              (figures/.cache/summary_agg.json)
  e  tiered_reclassification.json             (revision_analyses/outputs/)
  f  survival_overall.json + km_curve_overall.csv

Out: figures/manuscript_numbered/Figure_5_reclassification_patterns.(png|pdf|tiff)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches

PROJ = Path("/Users/danielliang/Downloads/ClinVar Research Project")
REV  = Path("/Users/danielliang/Library/CloudStorage/Dropbox/ClinVar_GenomMedicine/revision_analyses")
OUT  = REV / "outputs"
# R3 (round-2 review): outputs now go to figures/r3 so the round-2 file in
# manuscript_numbered/ is preserved untouched. Panel f: y-axis extended to 100%
# and the annotation rephrased per Reviewer 2 (round 2, Major #2).
FIGN = REV / "figures" / "r4"
TOTAL = 3678303

CATEGORY_ORDER = ["P", "LP", "P/LP", "VUS", "Conflicting", "B/LB", "LB", "B", "Other"]
NAT = {"blue": "#0064B5", "dark_red": "#B30024", "grey": "#7A7A7A",
       "yellow": "#FFB300", "lt_grey": "#BDBDBD"}

# ---- one typographic scale for the whole figure (readable at 100% zoom) ----
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "axes.labelsize": 9.5, "xtick.labelsize": 8.2, "ytick.labelsize": 8.2,
    "legend.fontsize": 8, "figure.dpi": 300, "savefig.dpi": 400,
})

def fmt(v):
    v = float(v)
    if v >= 1e6: return f"{v/1e6:.1f}M"
    if v >= 1e3: return f"{v/1e3:.1f}K"
    return f"{int(v)}"

def panel_label(ax, letter, x=-0.16, y=1.14):
    # place the letter in figure space to the upper-left of the axes so it never
    # sits under a (centered) panel title, regardless of title length
    ax.annotate(letter, xy=(x, y), xycoords="axes fraction",
                fontsize=13, fontweight="bold", va="top", ha="left")

# ---------------------------------------------------------------- data
summary = json.load(open(PROJ / "figures/.cache/summary_agg.json"))
stats   = json.load(open(PROJ / "data/statistics_20260219_024805.json"))
tj      = json.load(open(OUT / "tiered_reclassification.json"))
sv      = json.load(open(OUT / "survival_overall.json"))["overall"]
svroot  = json.load(open(OUT / "survival_overall.json"))
km      = pd.read_csv(OUT / "km_curve_overall.csv")

# =============================================================== figure
fig = plt.figure(figsize=(7.5, 10.6))
gs = GridSpec(3, 2, figure=fig, hspace=0.62, wspace=0.62,
              left=0.115, right=0.955, top=0.955, bottom=0.085)

# ---- (a) transition matrix -------------------------------------------
cats = CATEGORY_ORDER
ax = fig.add_subplot(gs[0, 0])
M = np.zeros((len(cats), len(cats)))
for key, count in summary["transition_matrix"].items():
    parts = key.split("|")
    if len(parts) == 2 and parts[0] in cats and parts[1] in cats:
        M[cats.index(parts[0]), cats.index(parts[1])] = count
masked = np.ma.masked_where(M == 0, M)
im = ax.imshow(masked, cmap="YlOrRd", norm=LogNorm(vmin=1, vmax=max(M.max(), 1)),
               aspect="auto", interpolation="nearest")
ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7)
ax.set_xlabel("Current category"); ax.set_ylabel("First category", labelpad=1)
ax.set_title("Classification transitions", pad=8)
for i in range(len(cats)):
    for j in range(len(cats)):
        v = M[i, j]
        if v > 0:
            ax.text(j, i, (f"{v/1e6:.1f}M" if v >= 1e6 else f"{v/1e3:.0f}K" if v >= 1e5 else fmt(v)),
                    ha="center", va="center", fontsize=5.0,
                    color="white" if v > M.max() * 0.15 else "black")
cbar = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.03)
cbar.ax.tick_params(labelsize=6); cbar.set_label("Variants (log scale)", fontsize=7, labelpad=1)
panel_label(ax, "a")

# ---- (b) top reclassification patterns -------------------------------
ax = fig.add_subplot(gs[0, 1])
cd = stats["change_type_distribution"]
top = sorted(cd.items(), key=lambda x: -x[1])[:18]
labels_b = [k for k, _ in reversed(top)]
values_b = [v for _, v in reversed(top)]
p_side, b_side = {"P", "LP", "P/LP"}, {"B", "LB", "B/LB"}
cols = []
for lbl in labels_b:
    st = lbl.split("->"); a, z = st[0], st[-1]
    cols.append(NAT["blue"] if a in p_side and z in b_side
                else NAT["dark_red"] if a in b_side and z in p_side else NAT["grey"])
y = np.arange(len(labels_b))
ax.barh(y, values_b, color=cols, height=0.72, edgecolor="white", linewidth=0.3)
ax.set_yticks(y); ax.set_yticklabels(labels_b, fontsize=6.5)
ax.set_xlabel("Variants"); ax.set_title("Top reclassification patterns", pad=8)
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
present = set(cols)
defs = [(NAT["dark_red"], "B-side → P-side"), (NAT["blue"], "P-side → B-side"),
        (NAT["grey"], "Lateral / mixed")]
patches = [mpatches.Patch(color=c, label=l) for c, l in defs if c in present]
if patches:
    ax.legend(handles=patches, loc="lower right", frameon=True, facecolor="white",
              edgecolor="none", framealpha=0.85, fontsize=6.8)
panel_label(ax, "b")

# ---- (c) reclassification overview (donut) — RE-LAID-OUT --------------
ax = fig.add_subplot(gs[1, 0])
vals = [summary["stable"], summary["changed"] - summary["oscillating"], summary["oscillating"]]
ov_labels = ["Stable", "Changed (non-oscillating)", "Changed (oscillating)"]
ov_colors = [NAT["grey"], NAT["yellow"], NAT["dark_red"]]
wedges, _, autos = ax.pie(
    vals, labels=None, autopct=lambda p: f"{p:.1f}%", colors=ov_colors,
    startangle=90, pctdistance=0.80, explode=(0, 0, 0.12),
    wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.2),
    textprops=dict(fontsize=8))
for at in autos:
    at.set_fontsize(8); at.set_fontweight("bold"); at.set_color("white")
autos[2].set_text("")                                   # oscillating wedge: external callout
osc_pct = vals[2] / sum(vals) * 100
ang = (wedges[2].theta1 + wedges[2].theta2) / 2
xp, yp = np.cos(np.radians(ang)), np.sin(np.radians(ang))
ax.annotate(f"{osc_pct:.1f}% oscillating", xy=(0.80 * xp, 0.80 * yp),
            xytext=(1.15 * xp + 0.15, 1.42 * yp), fontsize=7.5, fontweight="bold",
            color=NAT["dark_red"], ha="center", va="center",
            arrowprops=dict(arrowstyle="-|>", color=NAT["dark_red"], lw=1.0))
ax.text(0, 0, f"{summary['total']:,}\ntotal", ha="center", va="center",
        fontsize=9.5, fontweight="bold", color="#333333")
ax.set_ylim(-1.25, 1.55)                                # headroom so title/callout don't collide
ax.set_title("Reclassification overview", pad=14)
ax.legend(wedges, ov_labels, loc="lower center", bbox_to_anchor=(0.5, -0.17),
          ncol=1, fontsize=7, frameon=False)
panel_label(ax, "c", x=-0.16, y=1.18)

# ---- (d) P/LP <-> B/LB reversals -------------------------------------
ax = fig.add_subplot(gs[1, 1])
p2b = sorted(summary.get("p_to_b", {}).items(), key=lambda x: -x[1])
b2p = sorted(summary.get("b_to_p", {}).items(), key=lambda x: -x[1])
items = ([(l, v, NAT["blue"]) for l, v in p2b[:6]] +
         [(l, v, NAT["dark_red"]) for l, v in b2p[:6]])
items.sort(key=lambda x: -x[1])
bl = [i[0] for i in reversed(items)]; bv = [i[1] for i in reversed(items)]; bc = [i[2] for i in reversed(items)]
y = np.arange(len(bl))
ax.barh(y, bv, color=bc, height=0.66, edgecolor="white", linewidth=0.3)
ax.set_yticks(y); ax.set_yticklabels(bl, fontsize=7)
ax.set_xlabel("Transition events"); ax.set_title("Pathogenic ↔ Benign reversals", pad=8)
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
ax.text(0.98, 0.06, f"P-side → B-side: {fmt(sum(dict(p2b).values()))}   "
        f"B-side → P-side: {fmt(sum(dict(b2p).values()))}",
        transform=ax.transAxes, fontsize=6.8, va="bottom", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=NAT["lt_grey"], linewidth=0.4))
panel_label(ax, "d")

# ---- (e) tiered clinical-significance --------------------------------
ax = fig.add_subplot(gs[2, 0])
part = tj["partition_by_most_severe_change"]
def find(sub): return next(v for k, v in part.items() if sub in k)
nvar = tj["reclassified_any_change"]["n"]
med, comp = find("Medically significant"), find("composite")
nonmed, adj = find("Non-medically"), find("Within-tier")
# the first bar is the TOTAL; the other four partition it (sum = the total). Separate them so
# they don't read as five peer categories (Reviewer comment on Figure 5e).
labels_e = ["Any change (total)", "Within-tier",
            "Non-medically signif.", "Composite/uncertain", "Medically signif."]
ns = [nvar, adj["n"], nonmed["n"], comp["n"], med["n"]]
vals_e = [100 * n / TOTAL for n in ns]
colors_e = ["#6b7785", "#d9d9d9", "#c6dbef", "#fdae6b", "#b2182b"]
xpos = [0, 1.55, 2.75, 3.95, 5.15]
ax.axvspan(0.98, 5.72, color="#f2f2f2", zorder=0)            # decomposition region
bars = ax.bar(xpos, vals_e, width=0.82, color=colors_e, edgecolor="white", zorder=3)
for b, n in zip(bars, ns):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.10,
            f"{n:,}\n({100*n/TOTAL:.2f}%)", ha="center", va="bottom", fontsize=6.0, zorder=4)
ax.set_xticks(xpos); ax.set_xticklabels(labels_e, fontsize=6.4, rotation=30, ha="right", rotation_mode="anchor")
ax.set_ylabel("% of all 3.68M ClinVar variants"); ax.set_ylim(0, 11.4); ax.set_xlim(-0.72, 5.85)
ax.grid(axis="y", alpha=0.25)
ax.text(0.77, 4.7, "=", ha="center", va="center", fontsize=13, color="#666", fontweight="bold", zorder=5)
ax.text(3.35, 11.1, "four parts sum to the total (9.30%)", ha="center", va="top",
        fontsize=6.0, style="italic", color="#555")
ax.set_title("Reclassification by clinical significance", pad=8, fontsize=8.8)
panel_label(ax, "e", x=-0.24)

# ---- (f) censored Kaplan-Meier ---------------------------------------
ax = fig.add_subplot(gs[2, 1])
cancer10 = svroot["cancer::Cancer genes (COSMIC CGC)"]["cuminc_10yr"]["est"]
noncanc10 = svroot["cancer::Non-cancer genes"]["cuminc_10yr"]["est"]
ax.step(km["age_years"], km["cuminc"] * 100, where="post", color="#2b3a67", lw=2.0)
ax.fill_between(km["age_years"], km["cuminc_lo"] * 100, km["cuminc_hi"] * 100,
                step="post", alpha=0.18, color="#2b3a67", linewidth=0)
for yr, key in [(1, "cuminc_1yr"), (5, "cuminc_5yr"), (10, "cuminc_10yr")]:
    est = sv[key]["est"] * 100
    ax.plot(yr, est, "o", color="#b2182b", ms=6, zorder=5)
    ax.annotate(f"{est:.1f}%", (yr, est), textcoords="offset points",
                xytext=(4, 10) if yr == 1 else (7, -11),
                fontsize=8, color="#b2182b", fontweight="bold")
ax.set_xlabel("Years since first ClinVar appearance")
ax.set_ylabel("Cumulative incidence of\nfirst reclassification (%)")
ax.set_title("Time to first reclassification (censored)", pad=8, fontsize=8.8)
ax.set_xlim(0, 13); ax.set_ylim(0, 100); ax.grid(alpha=0.25)
ax.axhline(50, color="#999999", lw=0.9, ls="--", zorder=1)
ax.text(12.8, 51.5, "50%", ha="right", va="bottom", fontsize=6.8, color="#777777")
ax.text(0.03, 0.97, "Cumulative incidence did not reach\n50% during follow-up; median\nnot estimable\n"
        f"10-yr: cancer genes {cancer10*100:.1f}%\nvs non-cancer {noncanc10*100:.1f}%",
        transform=ax.transAxes, va="top", ha="left", fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="#bbbbbb"))
panel_label(ax, "f", x=-0.24)

# =============================================================== save
FIGN.mkdir(parents=True, exist_ok=True)
base = FIGN / "Figure_5_reclassification_patterns"
fig.savefig(base.with_suffix(".png"))
fig.savefig(base.with_suffix(".pdf"))
try:
    fig.savefig(base.with_suffix(".tiff"), pil_kwargs={"compression": "tiff_lzw"})
except Exception as e:
    print("tiff skipped:", e)
plt.close(fig)
print("wrote:", base.with_suffix(".png"))
print("figure size:", fig.get_size_inches(), "in @ 400 dpi")
