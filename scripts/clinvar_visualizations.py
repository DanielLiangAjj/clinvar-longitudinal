#!/usr/bin/env python3
"""
ClinVar Longitudinal Evolution — Publication-Quality Visualizations
===================================================================
Generates Nature-journal-quality figures addressing research questions Q1-Q7
on the structural and content evolution of ClinVar XML releases.

Requirements:
    pip install matplotlib numpy

Data inputs (from clinvar_ftp_only_v5.py output):
    - summary_*.csv
    - classification_history_*.csv
    - statistics_*.json

Optional (from clinvar_xml_structure_analyzer.py):
    - xml_structure_analysis.json

Usage:
    python3 clinvar_visualizations.py \\
        --data-dir ~/Downloads \\
        --output-dir ~/Downloads/clinvar_figures

First run streams through CSVs (~25-40 min for 143M-row history file)
and caches aggregated data. Subsequent runs load from cache instantly.
"""

import os
import sys
import json
import csv
import time
import argparse
from pathlib import Path
from collections import Counter, defaultdict, OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm, Normalize
from matplotlib.gridspec import GridSpec


# ═══════════════════════════════════════════════════════════════════════
# STYLE & COLOUR CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

def set_nature_style():
    """Configure matplotlib for Nature-journal aesthetics."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica",
                            "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "axes.titleweight": "bold",
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6,
        "axes.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
        "legend.frameon": False,
        "legend.borderpad": 0.3,
        "axes.grid": False,
    })


# ── Nature Reviews colour palette (CMYK → hex) ──────────────────
NAT = {
    "blue":      "#0064B5",   # True Blue  (C100 M38 Y0  K8)
    "dark_red":  "#B30024",   # Dark Red   (C0   M100 Y80 K30)
    "pink":      "#D94B8C",   # Pink       (C12  M83  Y10 K0)
    "salmon":    "#FF5A50",   # Salmon     (C0   M73  Y50 K0)
    "purple":    "#5C0099",   # Purple     (C80  M100 Y0  K0)
    "violet":    "#7C2DB8",   # Violet     (C53  M88  Y0  K0)
    "aqua":      "#007A6E",   # Aqua       (C85  M0   Y50 K10)
    "lime":      "#4D9900",   # Lime Green (C60  M10  Y100 K10)
    "yellow":    "#FFB300",   # Yellow     (C0   M29  Y100 K0)
    "grey":      "#7A7A7A",   # Grey       (C5   M5   Y5   K50)
    "cyan":      "#005C99",   # Cyan       (C80  M20  Y0   K20)
    "lt_grey":   "#BDBDBD",   # Light grey (neutral)
}

# Classification categories — warm tones for pathogenic, cool for benign
CATEGORY_ORDER = ["P", "LP", "P/LP", "VUS", "Conflicting", "B/LB", "LB", "B", "Other"]

CATEGORY_COLORS = {
    "P":           NAT["dark_red"],
    "LP":          NAT["salmon"],
    "P/LP":        NAT["pink"],
    "VUS":         NAT["yellow"],
    "Conflicting": NAT["violet"],
    "LB":          NAT["blue"],
    "B":           NAT["cyan"],
    "B/LB":        NAT["aqua"],
    "Other":       NAT["grey"],
}

CATEGORY_LABELS = {
    "P": "Pathogenic",
    "LP": "Likely pathogenic",
    "P/LP": "Path./Likely path.",
    "VUS": "VUS",
    "Conflicting": "Conflicting",
    "LB": "Likely benign",
    "B": "Benign",
    "B/LB": "Benign/Likely benign",
    "Other": "Other",
}

# Review-status grouping — darker = higher review authority
REVIEW_ORDER = [
    "Practice guideline",
    "Expert panel",
    "Multiple, concordant",
    "Conflicting",
    "Single submitter",
    "No criteria",
    "Not provided",
]

REVIEW_COLORS = {
    "Practice guideline":   NAT["cyan"],
    "Expert panel":         NAT["blue"],
    "Multiple, concordant": NAT["lime"],
    "Conflicting":          NAT["violet"],
    "Single submitter":     NAT["yellow"],
    "No criteria":          NAT["salmon"],
    "Not provided":         NAT["lt_grey"],
}


# ═══════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def categorize_classification(cls_str):
    """Mirror the categorisation from clinvar_ftp_only_v5.py."""
    if not cls_str:
        return "Other"
    c = cls_str.lower().strip()
    if "conflicting" in c:
        return "Conflicting"
    if "uncertain" in c or "vus" in c:
        return "VUS"
    if "pathogenic/likely" in c or "pathogenic, likely" in c:
        return "P/LP"
    if "benign/likely" in c or "benign, likely" in c:
        return "B/LB"
    if "likely pathogenic" in c:
        return "LP"
    if "pathogenic" in c:
        return "P"
    if "likely benign" in c:
        return "LB"
    if "benign" in c:
        return "B"
    return "Other"


def categorize_review_status(rs):
    """Group review-status strings into standard tiers."""
    if not rs:
        return "Not provided"
    r = rs.lower().strip()
    if "practice guideline" in r:
        return "Practice guideline"
    if "expert panel" in r or "professional society" in r:
        return "Expert panel"
    if "classified by multiple" in r:
        return "Multiple, concordant"
    if "multiple submitters" in r and "no conflicts" in r:
        return "Multiple, concordant"
    if "conflicting" in r:
        return "Conflicting"
    if "single submitter" in r or "classified by single" in r:
        return "Single submitter"
    if "criteria provided" in r:
        return "Single submitter"
    if any(k in r for k in ("no assertion", "no classification",
                             "not classified", "flagged")):
        return "No criteria"
    return "Not provided"


def date_to_float(date_str):
    """Convert '2012-11' to 2012.833 for accurate temporal spacing."""
    parts = date_str.split("-")
    year = int(parts[0])
    month = int(parts[1]) if len(parts) > 1 else 1
    return year + (month - 1) / 12.0


def fmt(n):
    """Format large numbers compactly: 1234567 -> '1.2M', 54321 -> '54.3K'."""
    if n >= 1_000_000:
        return f"{n / 1e6:.1f}M"
    if n >= 1_000:
        return f"{n / 1e3:.1f}K"
    return str(int(n))


def panel_label(ax, label, x=-0.14, y=1.08):
    """Add bold lowercase panel label (Nature convention)."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top", ha="left")


def find_file(directory, pattern):
    """Find the most recent file matching a glob pattern."""
    matches = sorted(directory.glob(pattern))
    return matches[-1] if matches else None


# ═══════════════════════════════════════════════════════════════════════
# DATA AGGREGATION — Streaming for large files, with disk cache
# ═══════════════════════════════════════════════════════════════════════

def aggregate_history(history_path, cache_path):
    """Stream classification_history CSV and aggregate per-snapshot statistics.

    For the 143 M-row file this takes ~25-40 min on first run;
    results are cached to JSON for instant reload.
    """
    if cache_path.exists():
        print(f"  [cache hit] Loading history aggregation from {cache_path.name}")
        with open(cache_path) as f:
            return json.load(f)

    print(f"  Streaming {history_path.name} ({history_path.stat().st_size / 1e9:.1f} GB) ...")
    t0 = time.time()
    snapshot_data = {}
    row_count = 0

    with open(history_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row_count += 1
            if row_count % 10_000_000 == 0:
                elapsed = time.time() - t0
                rate = row_count / elapsed
                eta = (143_000_000 - row_count) / rate if rate > 0 else 0
                print(f"    {row_count / 1e6:.0f}M rows  "
                      f"({rate / 1e6:.2f}M/s, ~{eta / 60:.0f} min remaining)")

            date = row["Snapshot_Date"]
            if date not in snapshot_data:
                snapshot_data[date] = {
                    "total": 0,
                    "categories": {},
                    "review_statuses": {},
                    "n_cls_strings": {},
                    "sub_sum": 0,
                    "multi_sub": 0,
                }
            sd = snapshot_data[date]
            sd["total"] += 1

            cat = row.get("Category", "Other") or "Other"
            sd["categories"][cat] = sd["categories"].get(cat, 0) + 1

            rs_raw = row.get("Review_Status", "")
            rs_group = categorize_review_status(rs_raw)
            sd["review_statuses"][rs_group] = sd["review_statuses"].get(rs_group, 0) + 1

            cls_str = row.get("Classification", "")
            if cls_str:
                sd["n_cls_strings"][cls_str] = sd["n_cls_strings"].get(cls_str, 0) + 1

            nsub = int(row.get("Num_Submissions", 0) or 0)
            sd["sub_sum"] += nsub
            if nsub > 1:
                sd["multi_sub"] += 1

    # Derive summary stats per snapshot
    for d, sd in snapshot_data.items():
        sd["sub_mean"] = round(sd["sub_sum"] / max(sd["total"], 1), 3)
        sd["pct_multi_sub"] = round(sd["multi_sub"] / max(sd["total"], 1) * 100, 2)
        sd["classification_diversity"] = len(sd["n_cls_strings"])
        # Drop bulky raw strings to keep cache small
        del sd["n_cls_strings"]

    elapsed = time.time() - t0
    print(f"    Done: {row_count:,} rows, {len(snapshot_data)} snapshots in {elapsed / 60:.1f} min")

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(snapshot_data, f)
    print(f"    Cached to {cache_path.name}")
    return snapshot_data


def aggregate_summary(summary_path, cache_path):
    """Stream summary CSV and aggregate statistics for figures 2 & 4."""
    if cache_path.exists():
        print(f"  [cache hit] Loading summary aggregation from {cache_path.name}")
        with open(cache_path) as f:
            return json.load(f)

    print(f"  Streaming {summary_path.name} ({summary_path.stat().st_size / 1e9:.1f} GB) ...")
    t0 = time.time()

    agg = {
        "total": 0,
        "changed": 0,
        "stable": 0,
        "oscillating": 0,
        "field_completeness": {
            "Gene": 0, "Variant_Type": 0, "Molecular_Consequence": 0,
            "HGVS_c": 0, "HGVS_p": 0, "dbSNP": 0,
            "Chromosome": 0, "Assembly": 0, "Review_Status": 0,
        },
        "transition_matrix": {},
        "pmid_bins": {"0": 0, "1": 0, "2-5": 0, "6-10": 0, "11-20": 0, ">20": 0},
        "acmg_bins": {"0": 0, "1-2": 0, "3-5": 0, "6-10": 0, ">10": 0},
        "first_date_counts": {},
        "dbsnp_by_first_date": {},
        "pmid_sum_by_first_date": {},
        "p_to_b": {},
        "b_to_p": {},
    }

    p_side = {"P", "LP", "P/LP"}
    b_side = {"B", "LB", "B/LB"}

    with open(summary_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            agg["total"] += 1
            if agg["total"] % 500_000 == 0:
                print(f"    {agg['total'] / 1e6:.1f}M rows ...")

            # Field completeness
            for field in agg["field_completeness"]:
                if row.get(field, "").strip():
                    agg["field_completeness"][field] += 1

            # Transition matrix
            cur_cat = row.get("Current_Category", "Other") or "Other"
            first_cat = categorize_classification(row.get("First_Classification", ""))
            key = f"{first_cat}|{cur_cat}"
            agg["transition_matrix"][key] = agg["transition_matrix"].get(key, 0) + 1

            # Changed / stable / oscillating
            changed = row.get("Classification_Changed", "") == "True"
            if changed:
                agg["changed"] += 1
            else:
                agg["stable"] += 1
            if row.get("Had_Oscillation", "") == "True":
                agg["oscillating"] += 1

            # PMID bins
            npm = int(row.get("Num_PMIDs", 0) or 0)
            if npm == 0:
                agg["pmid_bins"]["0"] += 1
            elif npm == 1:
                agg["pmid_bins"]["1"] += 1
            elif npm <= 5:
                agg["pmid_bins"]["2-5"] += 1
            elif npm <= 10:
                agg["pmid_bins"]["6-10"] += 1
            elif npm <= 20:
                agg["pmid_bins"]["11-20"] += 1
            else:
                agg["pmid_bins"][">20"] += 1

            # ACMG bins
            nacmg = int(row.get("Num_ACMG_Codes", 0) or 0)
            if nacmg == 0:
                agg["acmg_bins"]["0"] += 1
            elif nacmg <= 2:
                agg["acmg_bins"]["1-2"] += 1
            elif nacmg <= 5:
                agg["acmg_bins"]["3-5"] += 1
            elif nacmg <= 10:
                agg["acmg_bins"]["6-10"] += 1
            else:
                agg["acmg_bins"][">10"] += 1

            # First-date tracking
            fd = row.get("First_Date", "")
            if fd:
                agg["first_date_counts"][fd] = agg["first_date_counts"].get(fd, 0) + 1
                if row.get("dbSNP", "").strip():
                    agg["dbsnp_by_first_date"][fd] = agg["dbsnp_by_first_date"].get(fd, 0) + 1
                agg["pmid_sum_by_first_date"][fd] = (
                    agg["pmid_sum_by_first_date"].get(fd, 0) + npm
                )

            # P/LP <-> B/LB direct transitions
            if changed:
                ct = row.get("Change_Type", "")
                if ct:
                    steps = ct.split("->")
                    for i in range(len(steps) - 1):
                        a, b = steps[i], steps[i + 1]
                        if a in p_side and b in b_side:
                            tk = f"{a}->{b}"
                            agg["p_to_b"][tk] = agg["p_to_b"].get(tk, 0) + 1
                        elif a in b_side and b in p_side:
                            tk = f"{a}->{b}"
                            agg["b_to_p"][tk] = agg["b_to_p"].get(tk, 0) + 1

    # Convert field completeness to percentages
    for field in agg["field_completeness"]:
        agg["field_completeness"][field] = round(
            agg["field_completeness"][field] / max(agg["total"], 1) * 100, 2
        )

    elapsed = time.time() - t0
    print(f"    Done: {agg['total']:,} rows in {elapsed:.0f}s")

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(agg, f, indent=2)
    print(f"    Cached to {cache_path.name}")
    return agg


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1 — ClinVar Growth & Classification Landscape (Q3, Q5)
# ═══════════════════════════════════════════════════════════════════════

def figure1(hist, out):
    """Four-panel figure: growth, proportions, absolute counts, diversity."""
    dates = sorted(hist.keys())
    x = np.array([date_to_float(d) for d in dates])
    totals = np.array([hist[d]["total"] for d in dates])

    # Build category arrays
    cat_abs = {}
    for cat in CATEGORY_ORDER:
        cat_abs[cat] = np.array([hist[d]["categories"].get(cat, 0) for d in dates])
    cat_pct = {}
    for cat in CATEGORY_ORDER:
        cat_pct[cat] = cat_abs[cat] / np.maximum(totals, 1) * 100

    diversity = np.array([hist[d]["classification_diversity"] for d in dates])

    fig, axes = plt.subplots(2, 2, figsize=(7.08, 5.0))
    fig.subplots_adjust(hspace=0.50, wspace=0.38)

    # ── (a) Database growth ──────────────────────────────────────────
    ax = axes[0, 0]
    ax.fill_between(x, 0, totals / 1e6, color=NAT["blue"], alpha=0.12)
    ax.plot(x, totals / 1e6, color=NAT["blue"], linewidth=1.3)
    # RCV→VCV transition
    ax.axvline(2020.0, color=NAT["dark_red"], linestyle="--", linewidth=0.6, alpha=0.7)
    ax.annotate("RCV\u2009\u2192\u2009VCV", xy=(2020.0, totals.max() / 1e6 * 0.93),
                fontsize=5.5, color=NAT["dark_red"], ha="left",
                xytext=(2020.3, totals.max() / 1e6 * 0.93))
    ax.annotate(f"{totals[-1] / 1e6:.1f}M", xy=(x[-1], totals[-1] / 1e6),
                fontsize=5.5, color=NAT["blue"], ha="right", va="bottom",
                xytext=(-4, 4), textcoords="offset points")
    ax.set_ylabel("Variants (millions)")
    ax.set_xlabel("Year")
    ax.set_title("Database growth")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.1f}"))
    panel_label(ax, "a")

    # ── (b) Classification proportions (stacked area) ────────────────
    ax = axes[0, 1]
    bottom = np.zeros_like(x)
    for cat in CATEGORY_ORDER:
        ax.fill_between(x, bottom, bottom + cat_pct[cat],
                        color=CATEGORY_COLORS[cat], linewidth=0.15,
                        edgecolor="white", alpha=0.88)
        bottom = bottom + cat_pct[cat]
    ax.set_ylim(0, 100)
    ax.set_ylabel("Proportion (%)")
    ax.set_xlabel("Year")
    ax.set_title("Classification distribution")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    panel_label(ax, "b")

    # ── (c) Absolute counts (stacked area) ───────────────────────────
    ax = axes[1, 0]
    bottom = np.zeros_like(x)
    for cat in CATEGORY_ORDER:
        vals = cat_abs[cat] / 1e6
        ax.fill_between(x, bottom, bottom + vals,
                        color=CATEGORY_COLORS[cat], linewidth=0.15,
                        edgecolor="white", alpha=0.88)
        bottom = bottom + vals
    ax.set_ylabel("Variants (millions)")
    ax.set_xlabel("Year")
    ax.set_title("Classification counts")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    panel_label(ax, "c")

    # ── (d) Classification diversity ─────────────────────────────────
    ax = axes[1, 1]
    ax.plot(x, diversity, color=NAT["lime"], linewidth=1.2)
    ax.fill_between(x, 0, diversity, color=NAT["lime"], alpha=0.08)
    ax.set_ylabel("Distinct classification strings")
    ax.set_xlabel("Year")
    ax.set_title("Classification vocabulary size")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    ax.annotate(f"{diversity[-1]}", xy=(x[-1], diversity[-1]),
                fontsize=5.5, color=NAT["lime"], ha="right", va="bottom",
                xytext=(-3, 4), textcoords="offset points")
    panel_label(ax, "d")

    # Shared legend below
    handles = [mpatches.Patch(color=CATEGORY_COLORS[c], label=CATEGORY_LABELS[c])
               for c in CATEGORY_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=5.5,
              bbox_to_anchor=(0.5, -0.04), frameon=False, columnspacing=1.0)

    _save(fig, out, "figure1_growth_classification_landscape")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2 — Reclassification Dynamics
# ═══════════════════════════════════════════════════════════════════════

def figure2(summary_agg, stats, out):
    """Four panels: transition matrix, top patterns, overview, P<->B reversals."""
    fig = plt.figure(figsize=(7.08, 6.0))
    gs = GridSpec(2, 2, figure=fig, hspace=0.50, wspace=0.50)

    cats = CATEGORY_ORDER  # 9 categories

    # ── (a) Transition matrix heatmap ────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    matrix = np.zeros((len(cats), len(cats)))
    for key, count in summary_agg["transition_matrix"].items():
        parts = key.split("|")
        if len(parts) == 2:
            fc, tc = parts
            if fc in cats and tc in cats:
                matrix[cats.index(fc), cats.index(tc)] = count

    # Mask zeros for clean display
    masked = np.ma.masked_where(matrix == 0, matrix)
    im = ax.imshow(masked, cmap="YlOrRd",
                   norm=LogNorm(vmin=1, vmax=max(matrix.max(), 1)),
                   aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats, fontsize=5.5)
    ax.set_xlabel("Current category", fontsize=7)
    ax.set_ylabel("First category", fontsize=7)
    ax.set_title("Classification transitions")

    # Annotate cells
    for i in range(len(cats)):
        for j in range(len(cats)):
            v = matrix[i, j]
            if v > 0:
                color = "white" if v > matrix.max() * 0.15 else "black"
                ax.text(j, i, fmt(v), ha="center", va="center",
                        fontsize=4, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.03)
    cbar.ax.tick_params(labelsize=5)
    cbar.set_label("Variants (log scale)", fontsize=6)
    panel_label(ax, "a", x=-0.22)

    # ── (b) Top reclassification patterns ────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    change_dist = stats.get("change_type_distribution", {})
    if isinstance(change_dist, dict):
        top_changes = sorted(change_dist.items(), key=lambda x: -x[1])[:18]
    else:
        top_changes = change_dist[:18]

    labels_b = [item[0] for item in reversed(top_changes)]
    values_b = [item[1] for item in reversed(top_changes)]

    # Colour-code by direction
    p_side = {"P", "LP", "P/LP"}
    b_side = {"B", "LB", "B/LB"}
    bar_colors = []
    for lbl in labels_b:
        steps = lbl.split("->")
        first_s, last_s = steps[0], steps[-1]
        if first_s in p_side and last_s in b_side:
            bar_colors.append(NAT["blue"])  # downgrade
        elif first_s in b_side and last_s in p_side:
            bar_colors.append(NAT["dark_red"])  # upgrade
        else:
            bar_colors.append(NAT["grey"])  # lateral / mixed

    y_pos = np.arange(len(labels_b))
    ax.barh(y_pos, values_b, color=bar_colors, height=0.7, edgecolor="white",
            linewidth=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels_b, fontsize=4.8)
    ax.set_xlabel("Variants")
    ax.set_title("Top reclassification patterns")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
    # Legend \u2014 only include categories actually present among these bars, so panel (b) does
    # not display the blue/red Pathogenic<->Benign colours that belong to panel (d) (R2 minor).
    present = set(bar_colors)
    cat_defs = [
        (NAT["dark_red"], "B-side \u2192 P-side"),
        (NAT["blue"], "P-side \u2192 B-side"),
        (NAT["grey"], "Lateral / mixed"),
    ]
    legend_patches = [mpatches.Patch(color=c, label=l) for c, l in cat_defs if c in present]
    if legend_patches:
        ax.legend(handles=legend_patches, fontsize=5, loc="lower right",
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "b", x=-0.25)

    # ── (c) Reclassification overview (donut chart) ──────────────────
    ax = fig.add_subplot(gs[1, 0])
    overview_vals = [
        summary_agg["stable"],
        summary_agg["changed"] - summary_agg["oscillating"],
        summary_agg["oscillating"],
    ]
    overview_labels = ["Stable", "Changed\n(non-oscillating)", "Changed\n(oscillating)"]
    overview_colors = [NAT["grey"], NAT["yellow"], NAT["dark_red"]]
    # Explode the tiny oscillating wedge so it's visible
    explode = (0, 0, 0.08)

    wedges, texts, autotexts = ax.pie(
        overview_vals, labels=None, autopct=lambda p: f"{p:.1f}%",
        colors=overview_colors, startangle=90, pctdistance=0.78,
        explode=explode,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.2),
        textprops=dict(fontsize=6),
    )
    for at in autotexts:
        at.set_fontsize(6)
        at.set_fontweight("bold")
        at.set_color("white")
    # The oscillating slice is too thin for internal label — use annotation
    # Override its autopct to empty, then add external callout
    autotexts[2].set_text("")
    osc_pct = overview_vals[2] / sum(overview_vals) * 100
    ang = (wedges[2].theta1 + wedges[2].theta2) / 2
    x_pt = np.cos(np.radians(ang))
    y_pt = np.sin(np.radians(ang))
    ax.annotate(
        f"{osc_pct:.1f}%\noscillating",
        xy=(0.75 * x_pt, 0.75 * y_pt),
        xytext=(1.35 * x_pt, 1.25 * y_pt),
        fontsize=5.5, fontweight="bold", color=NAT["dark_red"],
        ha="center", va="center",
        arrowprops=dict(arrowstyle="-|>", color=NAT["dark_red"], lw=0.8),
    )
    # External legend for categories
    ax.legend(wedges, overview_labels, loc="lower center",
              bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=5, frameon=False)
    # Centre annotation
    ax.text(0, 0, f"{summary_agg['total']:,}\ntotal", ha="center", va="center",
            fontsize=8, fontweight="bold", color="#333333")
    ax.set_title("Reclassification overview")
    panel_label(ax, "c", x=-0.12)

    # ── (d) P/LP <-> B/LB reversals ─────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    p2b = summary_agg.get("p_to_b", {})
    b2p = summary_agg.get("b_to_p", {})

    # Aggregate into broader groups
    p2b_total = sum(p2b.values())
    b2p_total = sum(b2p.values())

    # Detailed breakdown
    p2b_detail = sorted(p2b.items(), key=lambda x: -x[1])
    b2p_detail = sorted(b2p.items(), key=lambda x: -x[1])

    all_items = (
        [(lbl, val, NAT["blue"]) for lbl, val in p2b_detail[:6]] +
        [(lbl, val, NAT["dark_red"]) for lbl, val in b2p_detail[:6]]
    )
    all_items.sort(key=lambda x: -x[1])

    bar_labels = [item[0] for item in reversed(all_items)]
    bar_vals = [item[1] for item in reversed(all_items)]
    bar_cols = [item[2] for item in reversed(all_items)]

    y_pos = np.arange(len(bar_labels))
    ax.barh(y_pos, bar_vals, color=bar_cols, height=0.65, edgecolor="white",
            linewidth=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(bar_labels, fontsize=5.5)
    ax.set_xlabel("Transition events")
    ax.set_title("Pathogenic \u2194 Benign reversals")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))

    # Annotate totals — place at bottom-right to avoid overlapping bars
    ax.text(0.98, 0.30,
            f"P-side \u2192 B-side: {fmt(p2b_total)}   B-side \u2192 P-side: {fmt(b2p_total)}",
            transform=ax.transAxes, fontsize=5.5, va="center", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=NAT["lt_grey"], linewidth=0.4))
    panel_label(ax, "d", x=-0.22)

    _save(fig, out, "figure2_reclassification_dynamics")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3 — Provenance & Submission Evolution (Q4)
# ═══════════════════════════════════════════════════════════════════════

def figure3(hist, out, concordance_data=None):
    """Provenance & submission evolution.
    3-panel (1×3) without concordance data, 4-panel (2×2) with it.
    """
    dates = sorted(hist.keys())
    x = np.array([date_to_float(d) for d in dates])

    if concordance_data and concordance_data.get("snapshots"):
        fig, axes = plt.subplots(2, 2, figsize=(7.08, 5.0))
        fig.subplots_adjust(hspace=0.55, wspace=0.42)
        ax_a = axes[0, 0]
        ax_b = axes[0, 1]
        ax_c = axes[1, 0]
        ax_d = axes[1, 1]
    else:
        fig, axes = plt.subplots(1, 3, figsize=(7.08, 2.6))
        fig.subplots_adjust(wspace=0.45)
        ax_a = axes[0]
        ax_b = axes[1]
        ax_c = axes[2]
        ax_d = None

    # ── (a) Review-status distribution over time ─────────────────────
    ax = ax_a
    rs_pct = {}
    for rs_cat in REVIEW_ORDER:
        rs_pct[rs_cat] = np.array([
            hist[d]["review_statuses"].get(rs_cat, 0) / max(hist[d]["total"], 1) * 100
            for d in dates
        ])

    bottom = np.zeros_like(x)
    for rs_cat in REVIEW_ORDER:
        ax.fill_between(x, bottom, bottom + rs_pct[rs_cat],
                        color=REVIEW_COLORS[rs_cat], linewidth=0.15,
                        edgecolor="white", alpha=0.88)
        bottom = bottom + rs_pct[rs_cat]

    ax.set_ylim(0, 100)
    ax.set_ylabel("Proportion (%)")
    ax.set_xlabel("Year")
    ax.set_title("Review status evolution")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    handles = [mpatches.Patch(color=REVIEW_COLORS[c], label=c) for c in REVIEW_ORDER]
    ax.legend(handles=handles, fontsize=4.5, loc="upper left",
              bbox_to_anchor=(0.02, 0.98), ncol=1, handlelength=1.0,
              handletextpad=0.3, borderpad=0.2, frameon=True,
              facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "a")

    # ── (b) Mean submissions per variant ─────────────────────────────
    ax = ax_b
    sub_mean = np.array([hist[d]["sub_mean"] for d in dates])
    ax.plot(x, sub_mean, color=NAT["blue"], linewidth=1.2)
    ax.fill_between(x, 0, sub_mean, color=NAT["blue"], alpha=0.08)
    ax.set_ylabel("Mean submissions / variant")
    ax.set_xlabel("Year")
    ax.set_title("Submission depth")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    ax.annotate(f"{sub_mean[-1]:.2f}", xy=(x[-1], sub_mean[-1]),
                fontsize=5.5, color=NAT["blue"], ha="right", va="bottom",
                xytext=(-3, 4), textcoords="offset points")
    panel_label(ax, "b")

    # ── (c) % variants with >1 submission ────────────────────────────
    ax = ax_c
    pct_multi = np.array([hist[d]["pct_multi_sub"] for d in dates])
    ax.plot(x, pct_multi, color=NAT["lime"], linewidth=1.2)
    ax.fill_between(x, 0, pct_multi, color=NAT["lime"], alpha=0.08)
    ax.set_ylabel("Variants with >1 submission (%)")
    ax.set_xlabel("Year")
    ax.set_title("Multi-submission prevalence")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    ax.annotate(f"{pct_multi[-1]:.1f}%", xy=(x[-1], pct_multi[-1]),
                fontsize=5.5, color=NAT["lime"], ha="right", va="bottom",
                xytext=(-3, 4), textcoords="offset points")
    panel_label(ax, "c")

    # ── (d) Concordance / discordance trend (if data available) ──────
    if ax_d is not None and concordance_data:
        ax = ax_d
        conc_snaps = sorted(concordance_data["snapshots"],
                            key=lambda s: s["date"])

        # Filter to only show VCV-era data where concordance is measurable
        # (RCV-era review_status vocabulary does not distinguish concordance)
        vcv_start = "2020-01"
        vcv_snaps = [s for s in conc_snaps if s["date"] >= vcv_start]
        all_snaps = conc_snaps  # keep full for multi-sub count overlay

        # Full timeline multi-sub overlay (dotted line)
        all_x = np.array([date_to_float(s["date"]) for s in all_snaps])
        all_multi = np.array([s.get("multi_sub_variants", 0) for s in all_snaps])

        if vcv_snaps:
            v_dates = [s["date"] for s in vcv_snaps]
            v_x = np.array([date_to_float(d) for d in v_dates])
            v_conc = np.array([s.get("concordance_rate_pct", 0) for s in vcv_snaps])
            v_disc = np.array([s.get("discordance_rate_pct", 0) for s in vcv_snaps])

            # Stacked area: concordant + discordant + excluded (VCV era only)
            ax.fill_between(v_x, 0, v_conc, color=NAT["lime"], alpha=0.55,
                            label="Concordant")
            ax.fill_between(v_x, v_conc, v_conc + v_disc,
                            color=NAT["dark_red"], alpha=0.55, label="Discordant")
            ax.fill_between(v_x, v_conc + v_disc, 100,
                            color=NAT["lt_grey"], alpha=0.35, label="Excluded")

        # Grey out RCV era with annotation
        rcv_end_x = date_to_float("2019-12")
        ax.axvspan(ax.get_xlim()[0] if ax.get_xlim()[0] < rcv_end_x
                   else date_to_float("2012-11"),
                   rcv_end_x, color="#eeeeee", alpha=0.8, zorder=0)
        ax.text(date_to_float("2016-06"), 50,
                "RCV era:\nconcordance not\nmeasurable from\nReview Status",
                ha="center", va="center", fontsize=5, color=NAT["grey"],
                style="italic")
        # Boundary line
        ax.axvline(x=date_to_float("2020-01"), color="black", linestyle="--",
                   linewidth=0.6, alpha=0.4)

        ax.set_xlim(date_to_float("2012-11"), date_to_float("2025-07"))
        ax.set_ylim(0, 100)
        ax.set_ylabel("% of multi-submitter variants")
        ax.set_xlabel("Year")
        ax.set_title("Inter-submitter concordance")
        ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

        # Latest values annotation
        if vcv_snaps:
            last = vcv_snaps[-1]
            ax.text(0.98, 0.05,
                    f"Latest ({last['date']}): "
                    f"{last.get('concordance_rate_pct', 0):.1f}% concordant, "
                    f"{last.get('discordance_rate_pct', 0):.1f}% discordant",
                    transform=ax.transAxes, fontsize=4.5, ha="right", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=NAT["lt_grey"],
                              lw=0.4, alpha=0.9))

        # Secondary axis: multi-sub variants (full timeline)
        ax2 = ax.twinx()
        ax2.plot(all_x, all_multi / 1e3, "k:", linewidth=0.8, alpha=0.4,
                 label="Multi-sub variants (K)")
        ax2.set_ylabel("Multi-sub variants (\u00d71K)", fontsize=6)
        ax2.tick_params(labelsize=5.5)
        ax2.spines["right"].set_visible(True)
        ax2.spines["right"].set_linewidth(0.5)

        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, fontsize=4.5,
                  loc="upper left",
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
        panel_label(ax, "d")

    _save(fig, out, "figure3_provenance_submission_evolution")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 4 — Evidence, Linkage & Data Completeness (Q6, Q7)
# ═══════════════════════════════════════════════════════════════════════

def figure4(summary_agg, stats, out):
    """Four panels: field completeness, PMID dist, new-variants timeline, dbSNP coverage."""
    fig, axes = plt.subplots(2, 2, figsize=(7.08, 5.2))
    fig.subplots_adjust(hspace=0.55, wspace=0.40)

    # ── (a) Field completeness (horizontal bars) ─────────────────────
    ax = axes[0, 0]
    fc = summary_agg["field_completeness"]
    fields = sorted(fc.keys(), key=lambda k: fc[k], reverse=True)
    vals = [fc[f] for f in fields]
    colours = [NAT["blue"] if v >= 80 else NAT["yellow"] if v >= 50 else NAT["dark_red"]
               for v in vals]

    y_pos = np.arange(len(fields))
    ax.barh(y_pos, vals, color=colours, height=0.6, edgecolor="white", linewidth=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(fields, fontsize=6)
    ax.set_xlabel("Completeness (%)")
    ax.set_xlim(0, 105)
    ax.set_title("Field-level data completeness")
    for i, v in enumerate(vals):
        ax.text(v + 0.8, i, f"{v:.1f}%", va="center", fontsize=5, color=NAT["grey"])
    panel_label(ax, "a", x=-0.22)

    # ── (b) PubMed citation distribution ─────────────────────────────
    ax = axes[0, 1]
    pmid_order = ["0", "1", "2-5", "6-10", "11-20", ">20"]
    pmid_vals = [summary_agg["pmid_bins"].get(k, 0) for k in pmid_order]
    pmid_pct = [v / max(summary_agg["total"], 1) * 100 for v in pmid_vals]

    bar_cols_pmid = [NAT["grey"]] + [NAT["blue"]] * 5
    ax.bar(range(len(pmid_order)), pmid_pct, color=bar_cols_pmid, width=0.65,
           edgecolor="white", linewidth=0.3)
    ax.set_xticks(range(len(pmid_order)))
    ax.set_xticklabels(pmid_order, fontsize=6)
    ax.set_xlabel("PubMed citations per variant")
    ax.set_ylabel("Variants (%)")
    ax.set_title("Literature evidence distribution")
    for i, (v, pct) in enumerate(zip(pmid_vals, pmid_pct)):
        ax.text(i, pct + 0.5, fmt(v), ha="center", fontsize=4.5, color=NAT["grey"])
    panel_label(ax, "b")

    # ── (c) New variants per year ────────────────────────────────────
    ax = axes[1, 0]
    fdc = summary_agg["first_date_counts"]
    year_counts = defaultdict(int)
    for date, count in fdc.items():
        year = date.split("-")[0]
        year_counts[year] += count
    years = sorted(year_counts.keys())
    year_vals = [year_counts[y] for y in years]

    ax.bar(range(len(years)), [v / 1e6 for v in year_vals],
           color=NAT["blue"], width=0.7, edgecolor="white", linewidth=0.3)
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("New variants (millions)")
    ax.set_xlabel("Year of first ClinVar appearance")
    ax.set_title("Variant accession timeline")
    panel_label(ax, "c")

    # ── (d) dbSNP coverage by entry cohort ───────────────────────────
    ax = axes[1, 1]
    dbsnp_by_fd = summary_agg["dbsnp_by_first_date"]
    year_dbsnp = defaultdict(int)
    year_total = defaultdict(int)
    for date, count in fdc.items():
        y = date.split("-")[0]
        year_total[y] += count
        year_dbsnp[y] += dbsnp_by_fd.get(date, 0)

    dbsnp_pct = [year_dbsnp[y] / max(year_total[y], 1) * 100 for y in years]

    ax.bar(range(len(years)), dbsnp_pct, color=NAT["aqua"], width=0.7,
           edgecolor="white", linewidth=0.3)
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("Variants with dbSNP ID (%)")
    ax.set_xlabel("Year of first ClinVar appearance")
    ax.set_title("dbSNP cross-reference coverage")
    ax.set_ylim(0, 105)
    panel_label(ax, "d")

    _save(fig, out, "figure4_evidence_linkage_completeness")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 5 — XML Structural Evolution (Q1, Q2, Q6)
# ═══════════════════════════════════════════════════════════════════════

def figure5(xml_data, out, hist_agg=None):
    """Revised six-panel figure: data model evolution with format eras.
    Only generated if xml_structure_analysis.json is available.
    """
    import re as _re
    snapshots = [s for s in xml_data["snapshots"]
                 if _re.match(r"^\d{4}-\d{2}", s["snapshot_date"])]
    lifecycle = [lc for lc in xml_data["element_lifecycle"]
                 if _re.match(r"^\d{4}-\d{2}", lc["snapshot_date"])]
    fa = xml_data["field_availability"]
    fa_avail = {k: v for k, v in fa["availability"].items()
                if _re.match(r"^\d{4}-\d{2}", k)}
    fa = {"fields": fa["fields"], "availability": fa_avail}

    dates = [s["snapshot_date"] for s in snapshots]
    x = np.array([date_to_float(d) for d in dates])

    # Format-era boundary
    vcv_start = 2020.0

    def _add_era_shading(ax):
        """Add light background shading for RCV vs VCV eras."""
        ylim = ax.get_ylim()
        ax.axvspan(x.min() - 0.5, vcv_start, alpha=0.04, color=NAT["blue"], zorder=0)
        ax.axvspan(vcv_start, x.max() + 0.5, alpha=0.04, color=NAT["yellow"], zorder=0)
        ax.axvline(vcv_start, color=NAT["grey"], linestyle="--", linewidth=0.5, alpha=0.6)
        ax.set_ylim(ylim)

    fig = plt.figure(figsize=(7.08, 7.5))
    gs = GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.42)

    # ── (a) Format era timeline + variant count ───────────────────────
    ax = fig.add_subplot(gs[0, 0])
    fmt_colors_map = {"rcv": NAT["blue"], "vcv": NAT["yellow"]}
    bar_y = []
    bar_c = []
    for s in snapshots:
        d = s["snapshot_date"]
        f = s.get("format", "")
        count = hist_agg.get(d, {}).get("total", 0) if hist_agg else 0
        bar_y.append(count / 1e6)
        bar_c.append(fmt_colors_map.get(f, "gray"))

    ax.bar(x, bar_y, width=0.07, color=bar_c, edgecolor="white", linewidth=0.2)
    ax.set_ylabel("Variants (millions)")
    ax.set_xlabel("Year")
    ax.set_title("Format eras & variant growth")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    ax.axvline(vcv_start, color="black", linestyle=":", linewidth=0.5, alpha=0.5)
    ylim = ax.get_ylim()
    ax.text(2016.0, ylim[1] * 0.90, "RCV\n(submission-\ncentric)",
            fontsize=5.5, color=NAT["blue"], ha="center", fontweight="bold")
    ax.text(2022.5, ylim[1] * 0.90, "VCV\n(variant-\ncentric)",
            fontsize=5.5, color=NAT["yellow"], ha="center", fontweight="bold")

    legend_patches = [mpatches.Patch(color=NAT["blue"], label="RCV format"),
                      mpatches.Patch(color=NAT["yellow"], label="VCV format")]
    ax.legend(handles=legend_patches, fontsize=5, loc="lower right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "a")

    # ── (b) Element & attribute count (dual-axis) ─────────────────────
    ax = fig.add_subplot(gs[0, 1])
    n_elem = [s["num_unique_elements"] for s in snapshots]
    n_attr = [s["total_unique_attributes"] for s in snapshots]

    ax.plot(x, n_elem, "o-", color=NAT["blue"], markersize=2.5, linewidth=1.0,
            label="Elements")
    ax.set_ylabel("Unique elements", color=NAT["blue"])
    ax.tick_params(axis="y", labelcolor=NAT["blue"])
    ax.set_xlabel("Year")
    ax.set_title("Schema complexity")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    ax2 = ax.twinx()
    ax2.plot(x, n_attr, "s-", color=NAT["lime"], markersize=2, linewidth=0.9,
             label="Attributes")
    ax2.set_ylabel("Unique attributes", color=NAT["lime"], fontsize=6)
    ax2.tick_params(axis="y", labelcolor=NAT["lime"], labelsize=5.5)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_linewidth(0.5)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=5, loc="lower right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    _add_era_shading(ax)
    panel_label(ax, "b")

    # ── (c) Improved field availability heatmap ──────────────────────
    ax = fig.add_subplot(gs[1, 0])
    field_names = fa["fields"]
    avail = fa["availability"]
    avail_dates = sorted(avail.keys())

    # Build binary matrix
    mat = np.zeros((len(field_names), len(avail_dates)))
    for j, d in enumerate(avail_dates):
        for i, f in enumerate(field_names):
            mat[i, j] = 1 if avail[d].get(f, False) else 0

    # Only show fields present in at least one snapshot
    present_mask = mat.sum(axis=1) > 0
    mat_filtered = mat[present_mask]
    fields_filtered = [f for f, m in zip(field_names, present_mask) if m]

    # Group by semantic blocks
    FIELD_GROUPS = {
        "Record Structure": {"VariationArchive", "ClinVarSet", "ReferenceClinVarAssertion",
                             "InterpretedRecord", "ClassifiedRecord", "SimpleAllele",
                             "VariationType", "RecordStatus"},
        "Classification": {"ClinicalSignificance", "Description", "ReviewStatus",
                           "GermlineClassification", "SomaticClinicalImpact",
                           "OncogenicityClassification", "Interpretations", "Interpretation",
                           "Classifications"},
        "Evidence": {"ObservedIn", "ObservedData", "FunctionalConsequence",
                     "Citation", "CitationList", "Method", "MethodType",
                     "Sample", "Comment", "AttributeSet", "Attribute"},
        "Provenance": {"ClinicalAssertion", "ClinVarAssertion", "ClinVarSubmissionID",
                       "Submitter", "RecordHistory", "Replaced", "Accession"},
        "Nomenclature": {"HGVSlist", "HGVS", "NucleotideExpression", "ProteinExpression",
                         "Gene", "Name", "XRef", "XRefList"},
    }

    # Sort fields into groups (ungrouped fields go at end)
    ordered_fields = []
    ordered_indices = []
    group_boundaries = []
    for group_name, group_fields in FIELD_GROUPS.items():
        start = len(ordered_fields)
        for i, f in enumerate(fields_filtered):
            if f in group_fields:
                ordered_fields.append(f)
                ordered_indices.append(i)
        if len(ordered_fields) > start:
            group_boundaries.append((group_name, start, len(ordered_fields)))

    # Add ungrouped
    grouped_set = set(ordered_fields)
    start_other = len(ordered_fields)
    for i, f in enumerate(fields_filtered):
        if f not in grouped_set:
            ordered_fields.append(f)
            ordered_indices.append(i)
    if len(ordered_fields) > start_other:
        group_boundaries.append(("Other", start_other, len(ordered_fields)))

    mat_ordered = mat_filtered[ordered_indices] if ordered_indices else mat_filtered

    ax.imshow(mat_ordered, aspect="auto", cmap="Blues", interpolation="nearest",
              vmin=0, vmax=1)
    ax.set_yticks(range(len(ordered_fields)))
    ax.set_yticklabels(ordered_fields, fontsize=3.5)

    # Draw group separators
    for gname, gstart, gend in group_boundaries:
        if gstart > 0:
            ax.axhline(gstart - 0.5, color=NAT["lt_grey"], linewidth=0.5, linestyle="-")

    # VCV transition line on heatmap
    vcv_idx = None
    for j, d in enumerate(avail_dates):
        if date_to_float(d) >= vcv_start:
            vcv_idx = j
            break
    if vcv_idx is not None:
        ax.axvline(vcv_idx - 0.5, color=NAT["dark_red"], linewidth=0.6, linestyle="--", alpha=0.6)

    tick_idx = list(range(0, len(avail_dates), max(1, len(avail_dates) // 8)))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([avail_dates[i] for i in tick_idx], rotation=45,
                        ha="right", fontsize=4.5)
    ax.set_xlabel("Release", fontsize=6)
    ax.set_title("Field availability (grouped)")
    panel_label(ax, "c")

    # ── (d) Element lifecycle ─────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    lc_dates = [lc["snapshot_date"] for lc in lifecycle]
    lc_x = np.array([date_to_float(d) for d in lc_dates])
    new_e = [lc["new_elements"] for lc in lifecycle]
    dep_e = [lc["deprecated_elements"] for lc in lifecycle]
    ret_rate = [lc["retention_rate"] for lc in lifecycle]

    ax.bar(lc_x, new_e, width=0.06, color=NAT["lime"], alpha=0.8, label="New elements")
    ax.bar(lc_x, [-d for d in dep_e], width=0.06, color=NAT["dark_red"], alpha=0.8,
           label="Deprecated")
    ax.axhline(0, color="black", linewidth=0.4)
    ax.set_ylabel("Element count")
    ax.set_xlabel("Year")
    ax.set_title("Element lifecycle")
    ax.legend(fontsize=5, loc="lower left",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    ax2 = ax.twinx()
    ax2.plot(lc_x[1:], ret_rate[1:], "k--", linewidth=0.7, alpha=0.5,
             label="Retention rate")
    ax2.set_ylabel("Retention rate (%)", fontsize=6)
    ax2.tick_params(labelsize=5.5)
    ax2.set_ylim(50, 105)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_linewidth(0.5)
    _add_era_shading(ax)
    panel_label(ax, "d")

    # ── (e) Schema change taxonomy ────────────────────────────────────
    ax = fig.add_subplot(gs[2, 0])
    taxonomy = _classify_schema_changes(xml_data)
    if taxonomy:
        tax_dates = [t["date"] for t in taxonomy]
        tax_x = np.array([date_to_float(d) for d in tax_dates])
        type_colors = {"additive": NAT["lime"], "deprecation": NAT["dark_red"],
                       "mixed": NAT["yellow"], "restructuring": NAT["violet"], "none": "#e0e0e0"}

        # Stacked bars showing change type
        for ctype in ["additive", "deprecation", "mixed", "restructuring"]:
            vals = [t["new_elements"] + t["deprecated_elements"]
                    if t["change_type"] == ctype else 0 for t in taxonomy]
            if any(v > 0 for v in vals):
                ax.bar(tax_x, vals, width=0.06, color=type_colors.get(ctype, "gray"),
                       alpha=0.8, label=ctype.capitalize())

        ax.set_ylabel("Elements changed")
        ax.set_xlabel("Year")
        ax.set_title("Schema change taxonomy")
        ax.legend(fontsize=4.5, loc="lower left",
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

        # Cumulative overlay
        ax2 = ax.twinx()
        cum_new = [t["cumulative_new"] for t in taxonomy]
        cum_dep = [t["cumulative_deprecated"] for t in taxonomy]
        ax2.plot(tax_x, cum_new, "-", color=NAT["lime"], linewidth=0.8, alpha=0.6,
                 label="Cum. new")
        ax2.plot(tax_x, cum_dep, "-", color=NAT["dark_red"], linewidth=0.8, alpha=0.6,
                 label="Cum. deprecated")
        ax2.set_ylabel("Cumulative", fontsize=6)
        ax2.tick_params(labelsize=5.5)
        ax2.spines["right"].set_visible(True)
        ax2.spines["right"].set_linewidth(0.5)
    _add_era_shading(ax)
    panel_label(ax, "e")

    # ── (f) Cross-database linkage ────────────────────────────────────
    ax = fig.add_subplot(gs[2, 1])
    db_names = ["dbSNP", "MedGen", "OMIM", "Orphanet", "UniProtKB"]
    db_colors_map = {"dbSNP": NAT["blue"], "MedGen": NAT["lime"], "OMIM": NAT["dark_red"],
                     "Orphanet": NAT["violet"], "UniProtKB": NAT["yellow"]}

    for db in db_names:
        db_vals = [s["xref_databases"].get(db, 0) for s in snapshots]
        if max(db_vals) > 0:
            ax.plot(x, db_vals, "o-", color=db_colors_map.get(db, "gray"),
                    markersize=2, linewidth=0.9, label=db)

    ax.set_ylabel("XRef count (per sample)")
    ax.set_xlabel("Year")
    ax.set_title("Cross-database linkage")
    ax.legend(fontsize=4.5, ncol=1, loc="center right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    _add_era_shading(ax)
    panel_label(ax, "f")

    _save(fig, out, "figure5_xml_structural_evolution")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 6 — Supplementary: ACMG Evidence & Gene Landscape
# ═══════════════════════════════════════════════════════════════════════

def figure6(summary_agg, stats, out):
    """Supplementary: ACMG distribution + top genes."""
    fig, axes = plt.subplots(1, 2, figsize=(7.08, 3.0))
    fig.subplots_adjust(wspace=0.45)

    # ── (a) ACMG evidence code distribution ──────────────────────────
    ax = axes[0]
    acmg_order = ["0", "1-2", "3-5", "6-10", ">10"]
    acmg_vals = [summary_agg["acmg_bins"].get(k, 0) for k in acmg_order]
    acmg_pct = [v / max(summary_agg["total"], 1) * 100 for v in acmg_vals]

    ax.bar(range(len(acmg_order)), acmg_pct, color=NAT["violet"], width=0.6,
           edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_xticks(range(len(acmg_order)))
    ax.set_xticklabels(acmg_order, fontsize=6)
    ax.set_xlabel("ACMG evidence codes per variant")
    ax.set_ylabel("Variants (%)")
    ax.set_title("ACMG/AMP evidence utilisation")
    for i, (v, pct) in enumerate(zip(acmg_vals, acmg_pct)):
        ax.text(i, pct + 0.3, fmt(v), ha="center", fontsize=4.5, color=NAT["grey"])
    panel_label(ax, "a")

    # ── (b) Top 25 genes by variant count ────────────────────────────
    ax = axes[1]
    genes = stats.get("genes", {})
    if isinstance(genes, dict):
        top_genes = sorted(genes.items(), key=lambda x: -x[1])[:25]
    else:
        top_genes = genes[:25]

    g_labels = [g[0] for g in reversed(top_genes)]
    g_vals = [g[1] for g in reversed(top_genes)]

    ax.barh(range(len(g_labels)), g_vals, color=NAT["blue"], height=0.7,
            edgecolor="white", linewidth=0.2)
    ax.set_yticks(range(len(g_labels)))
    ax.set_yticklabels(g_labels, fontsize=4.5)
    ax.set_xlabel("Variants")
    ax.set_title("Top 25 genes by variant count")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
    panel_label(ax, "b", x=-0.20)

    _save(fig, out, "figure6_acmg_gene_landscape")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 7 — Supplementary: PubMed Citation Growth (Q7 temporal)
# ═══════════════════════════════════════════════════════════════════════

def figure7(summary_agg, out):
    """Temporal evolution of PubMed citation density by variant entry cohort."""
    fig, axes = plt.subplots(1, 2, figsize=(7.08, 2.8))
    fig.subplots_adjust(wspace=0.40)

    fdc = summary_agg["first_date_counts"]
    pmid_sum = summary_agg["pmid_sum_by_first_date"]

    # Aggregate by year
    year_total = defaultdict(int)
    year_pmid_sum = defaultdict(int)
    for date, count in fdc.items():
        y = date.split("-")[0]
        year_total[y] += count
        year_pmid_sum[y] += pmid_sum.get(date, 0)

    years = sorted(year_total.keys())
    mean_pmid_per_variant = [
        year_pmid_sum[y] / max(year_total[y], 1) for y in years
    ]
    total_pmids_cumul = np.cumsum([year_pmid_sum[y] for y in years])

    # ── (a) Mean PubMed citations per variant by entry year ──────────
    ax = axes[0]
    ax.bar(range(len(years)), mean_pmid_per_variant, color=NAT["dark_red"],
           width=0.7, edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("Mean PubMed citations / variant")
    ax.set_xlabel("Year of first appearance")
    ax.set_title("Literature evidence by entry cohort")
    panel_label(ax, "a")

    # ── (b) Cumulative PubMed citations ──────────────────────────────
    ax = axes[1]
    ax.fill_between(range(len(years)), 0, total_pmids_cumul / 1e6,
                    color=NAT["dark_red"], alpha=0.12)
    ax.plot(range(len(years)), total_pmids_cumul / 1e6,
            color=NAT["dark_red"], linewidth=1.2)
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("Cumulative PubMed citations (millions)")
    ax.set_xlabel("Year")
    ax.set_title("Literature evidence accumulation")
    ax.annotate(f"{total_pmids_cumul[-1] / 1e6:.2f}M",
                xy=(len(years) - 1, total_pmids_cumul[-1] / 1e6),
                fontsize=5.5, color=NAT["dark_red"], ha="right", va="bottom",
                xytext=(-3, 4), textcoords="offset points")
    panel_label(ax, "b")

    _save(fig, out, "figure7_pubmed_citation_growth")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 8 — Gap Coverage: Q2 (%), Q3 (evidence tags), Q4 (provenance),
#             Q5 (SCV/RCV/VCV data-model transition)
# ═══════════════════════════════════════════════════════════════════════

# Curated element sets for domain-specific tracking
EVIDENCE_ELEMENTS = {
    # Classification / interpretation
    "Classification", "GermlineClassification", "SomaticClinicalImpact",
    "OncogenicityClassification", "ClinicalSignificance",
    "Interpretation", "Interpretations",
    "Description", "ReviewStatus",
    # Observation / methodology
    "FunctionalConsequence", "ObservedIn", "ObservedData",
    "Sample", "Method", "MethodType",
    # Citations & comments
    "Citation", "CitationList", "Comment",
    # Trait / condition
    "TraitSet", "Trait", "TraitRelationship",
    # ACMG / assertion attributes
    "AttributeSet", "Attribute",
}

PROVENANCE_ELEMENTS = {
    # Submission identity
    "ClinVarSubmissionID", "ClinicalAssertion", "ClinVarAssertion",
    "Submitter", "SubmitterIdentifiers", "SubmissionName",
    # Versioning / history
    "RecordHistory", "Replaced", "RecordStatus",
    "ClassifiedRecord", "ClassificationComment",
    # Accession tracking
    "Accession", "ClassificationScore",
}


def figure8(xml_data, hist_agg, out):
    """Three-panel figure: additive vs transformative change,
    evidence granularity, provenance complexity.
    (RCV→VCV transition panel moved to figure5.)
    """
    import re as _re

    snapshots = [s for s in xml_data["snapshots"]
                 if _re.match(r"^\d{4}-\d{2}", s["snapshot_date"])]
    lifecycle = [lc for lc in xml_data["element_lifecycle"]
                 if _re.match(r"^\d{4}-\d{2}", lc["snapshot_date"])]

    dates = [s["snapshot_date"] for s in snapshots]
    x = np.array([date_to_float(d) for d in dates])

    fig, axes = plt.subplots(1, 3, figsize=(7.08, 2.8))
    fig.subplots_adjust(wspace=0.50)

    # ── (a) Additive vs. transformative change ────────────────────────
    ax = axes[0]
    lc_dates = [lc["snapshot_date"] for lc in lifecycle]
    lc_x = np.array([date_to_float(d) for d in lc_dates])
    pct_new = [lc["pct_new"] for lc in lifecycle]
    pct_dep = [lc["pct_deprecated"] for lc in lifecycle]
    retention = [lc["retention_rate"] for lc in lifecycle]

    ax.plot(lc_x[1:], pct_new[1:], "o-", color=NAT["lime"], markersize=2,
            linewidth=0.9, label="% new elements")
    ax.plot(lc_x[1:], pct_dep[1:], "s-", color=NAT["dark_red"], markersize=2,
            linewidth=0.9, label="% deprecated")
    ax.plot(lc_x[1:], retention[1:], "^-", color=NAT["blue"], markersize=1.8,
            linewidth=0.8, alpha=0.6, label="Retention rate")
    ax.set_ylabel("Percentage (%)")
    ax.set_xlabel("Year")
    ax.set_title("Additive vs. transformative change")
    ax.legend(fontsize=5, loc="center right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
    panel_label(ax, "a")

    # ── (b) Evidence element granularity ──────────────────────────────
    ax = axes[1]
    evidence_counts = []
    total_evidence_attrs = []
    for s in snapshots:
        elem_set = set(s["elements"])
        ev_present = elem_set & EVIDENCE_ELEMENTS
        evidence_counts.append(len(ev_present))
        attr_count = 0
        for e in ev_present:
            attr_count += len(s.get("attributes_per_element", {}).get(e, []))
        total_evidence_attrs.append(attr_count)

    ax.plot(x, evidence_counts, "o-", color=NAT["violet"], markersize=2.5,
            linewidth=1.0, label="Evidence elements")
    ax.set_ylabel("Evidence elements present", color=NAT["violet"])
    ax.tick_params(axis="y", labelcolor=NAT["violet"])
    ax.set_xlabel("Year")
    ax.set_title("Evidence representation granularity")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    ax2 = ax.twinx()
    ax2.plot(x, total_evidence_attrs, "s-", color=NAT["yellow"], markersize=2,
             linewidth=0.8, alpha=0.7, label="Evidence attributes")
    ax2.set_ylabel("Evidence attributes", color=NAT["yellow"], fontsize=6)
    ax2.tick_params(axis="y", labelcolor=NAT["yellow"], labelsize=5.5)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_linewidth(0.5)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=5, loc="lower right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "b")

    # ── (c) Provenance & versioning complexity ────────────────────────
    ax = axes[2]
    prov_elem_counts = []
    VERSION_ATTRS = {"DateCreated", "DateLastUpdated", "Version",
                     "SubmissionCount", "Accession", "DateChanged"}
    version_attr_counts = []

    for s in snapshots:
        elem_set = set(s["elements"])
        prov_present = elem_set & PROVENANCE_ELEMENTS
        prov_elem_counts.append(len(prov_present))
        va = 0
        for e, attrs in s.get("attributes_per_element", {}).items():
            attr_set = set(attrs)
            va += len(attr_set & VERSION_ATTRS)
        version_attr_counts.append(va)

    ax.plot(x, prov_elem_counts, "o-", color=NAT["aqua"], markersize=2.5,
            linewidth=1.0, label="Provenance elements")
    ax.set_ylabel("Provenance elements", color=NAT["aqua"])
    ax.tick_params(axis="y", labelcolor=NAT["aqua"])
    ax.set_xlabel("Year")
    ax.set_title("Provenance & versioning complexity")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    ax2 = ax.twinx()
    ax2.plot(x, version_attr_counts, "s-", color=NAT["pink"], markersize=2,
             linewidth=0.8, alpha=0.7, label="Versioning attributes")
    ax2.set_ylabel("Versioning attributes", color=NAT["pink"], fontsize=6)
    ax2.tick_params(axis="y", labelcolor=NAT["pink"], labelsize=5.5)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_linewidth(0.5)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=5, loc="lower right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "c")

    _save(fig, out, "figure8_data_model_evidence_provenance")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 9 — ACMG/AMP Criteria Analysis (from extractor output)
# ═══════════════════════════════════════════════════════════════════════

# Strength class colour map (warm = pathogenic, cool = benign)
ACMG_STRENGTH_COLORS = {
    "Very Strong Pathogenic": NAT["dark_red"],
    "Strong Pathogenic":      NAT["salmon"],
    "Moderate Pathogenic":    NAT["yellow"],
    "Supporting Pathogenic":  NAT["pink"],
    "Stand-alone Benign":     NAT["cyan"],
    "Strong Benign":          NAT["blue"],
    "Supporting Benign":      NAT["aqua"],
    "Other":                  NAT["grey"],
}

ACMG_STRENGTH_ORDER = [
    "Very Strong Pathogenic", "Strong Pathogenic", "Moderate Pathogenic",
    "Supporting Pathogenic", "Stand-alone Benign", "Strong Benign",
    "Supporting Benign",
]

def _acmg_strength_class(code):
    """Map an ACMG code to its strength class."""
    if code.startswith("PVS"):
        return "Very Strong Pathogenic"
    if code.startswith("PS"):
        return "Strong Pathogenic"
    if code.startswith("PM"):
        return "Moderate Pathogenic"
    if code.startswith("PP"):
        return "Supporting Pathogenic"
    if code.startswith("BA"):
        return "Stand-alone Benign"
    if code.startswith("BS"):
        return "Strong Benign"
    if code.startswith("BP"):
        return "Supporting Benign"
    return "Other"


def figure9(acmg_stats, out):
    """ACMG/AMP criteria analysis — 4-panel figure.

    (a) Top 20 individual criteria by frequency (horizontal bar, coloured by strength)
    (b) Pathogenic vs Benign criteria usage (grouped bar by strength class)
    (c) Extraction method: explicit vs inferred (stacked bar per criteria)
    (d) Criteria usage by classification category (heatmap)
    """
    fig = plt.figure(figsize=(7.08, 7.0))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.50)

    all_freq = acmg_stats.get("all_code_freq", {})
    explicit_freq = acmg_stats.get("explicit_code_freq", {})
    inferred_freq = acmg_stats.get("inferred_code_freq", {})
    code_by_cat = acmg_stats.get("code_by_category", {})

    # ── (a) Top 20 criteria by frequency ─────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    top_codes = sorted(all_freq.items(), key=lambda x: -x[1])[:20]
    if top_codes:
        labels_a = [c[0] for c in reversed(top_codes)]
        values_a = [c[1] for c in reversed(top_codes)]
        colors_a = [ACMG_STRENGTH_COLORS.get(_acmg_strength_class(c), NAT["grey"])
                     for c in labels_a]

        y_pos = np.arange(len(labels_a))
        ax.barh(y_pos, values_a, color=colors_a, height=0.7,
                edgecolor="white", linewidth=0.3)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels_a, fontsize=6)
        ax.set_xlabel("Submissions")
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
        for i, v in enumerate(values_a):
            ax.text(v + max(values_a) * 0.01, i, fmt(v),
                    va="center", fontsize=4.5, color=NAT["grey"])
    ax.set_title("Top 20 ACMG/AMP criteria")

    # Compact legend for strength classes
    seen_classes = set()
    legend_patches = []
    for code in labels_a if top_codes else []:
        sc = _acmg_strength_class(code)
        if sc not in seen_classes:
            seen_classes.add(sc)
            legend_patches.append(
                mpatches.Patch(color=ACMG_STRENGTH_COLORS.get(sc, NAT["grey"]),
                               label=sc.replace(" Pathogenic", " (P)")
                                      .replace(" Benign", " (B)")))
    if legend_patches:
        ax.legend(handles=legend_patches, fontsize=4.5, loc="lower right",
                  handlelength=1.0, handletextpad=0.3,
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "a", x=-0.22)

    # ── (b) Pathogenic vs Benign grouped bar ─────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    path_classes = ["Very Strong Pathogenic", "Strong Pathogenic",
                    "Moderate Pathogenic", "Supporting Pathogenic"]
    benign_classes = ["Stand-alone Benign", "Strong Benign", "Supporting Benign"]
    short_labels_p = ["PVS", "PS", "PM", "PP"]
    short_labels_b = ["BA", "BS", "BP"]

    # Aggregate by strength class
    path_vals = []
    for cls in path_classes:
        total = sum(v for c, v in all_freq.items() if _acmg_strength_class(c) == cls)
        path_vals.append(total)
    benign_vals = []
    for cls in benign_classes:
        total = sum(v for c, v in all_freq.items() if _acmg_strength_class(c) == cls)
        benign_vals.append(total)

    # Plot pathogenic bars
    n_p = len(path_classes)
    n_b = len(benign_classes)
    x_pos = np.arange(max(n_p, n_b))
    bar_w = 0.35

    # Pathogenic group
    p_x = np.arange(n_p)
    bars_p = ax.bar(p_x - bar_w / 2, path_vals, bar_w, color=NAT["dark_red"],
                    edgecolor="white", linewidth=0.3, label="Pathogenic criteria", alpha=0.85)
    ax.set_xticks(p_x)
    ax.set_xticklabels(short_labels_p, fontsize=5.5)

    # Add benign on right side with offset
    b_start = n_p + 0.5
    b_x = np.arange(n_b) + b_start
    bars_b = ax.bar(b_x, benign_vals, bar_w, color=NAT["blue"],
                    edgecolor="white", linewidth=0.3, label="Benign criteria", alpha=0.85)

    # Update x-ticks to include both groups
    all_ticks = list(p_x) + list(b_x)
    all_labels = short_labels_p + short_labels_b
    ax.set_xticks(all_ticks)
    ax.set_xticklabels(all_labels, fontsize=6)
    ax.set_ylabel("Submissions")
    ax.set_title("Pathogenic vs Benign criteria usage")
    # Add strength class annotations below x-axis
    strength_p = ["Very\nStrong", "Strong", "Moderate", "Supporting"]
    strength_b = ["Stand-\nalone", "Strong", "Supporting"]
    for i, (tick, s) in enumerate(zip(p_x, strength_p)):
        ax.text(tick, -0.15, s, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=3.5, color=NAT["grey"], style="italic")
    for i, (tick, s) in enumerate(zip(b_x, strength_b)):
        ax.text(tick, -0.15, s, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=3.5, color=NAT["grey"], style="italic")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
    ax.legend(fontsize=5.5, loc="upper right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)

    # Add counts on top
    for bar in list(bars_p) + list(bars_b):
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, h,
                    fmt(h), ha="center", va="bottom", fontsize=4, color=NAT["grey"])
    panel_label(ax, "b")

    # ── (c) Explicit vs Inferred stacked bar ────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    top15 = sorted(all_freq.items(), key=lambda x: -x[1])[:15]
    if top15:
        codes_c = [c[0] for c in reversed(top15)]
        expl_vals = [explicit_freq.get(c, 0) for c in codes_c]
        inf_vals = [inferred_freq.get(c, 0) for c in codes_c]

        y_pos = np.arange(len(codes_c))
        ax.barh(y_pos, expl_vals, height=0.65, color=NAT["blue"],
                edgecolor="white", linewidth=0.3, label="Explicit (regex)")
        ax.barh(y_pos, inf_vals, height=0.65, left=expl_vals,
                color=NAT["yellow"], edgecolor="white", linewidth=0.3,
                label="NLP-inferred")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(codes_c, fontsize=6)
        ax.set_xlabel("Submissions")
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
        ax.legend(fontsize=5.5, loc="lower right",
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    ax.set_title("Extraction method breakdown")
    panel_label(ax, "c", x=-0.22)

    # ── (d) Criteria x Classification heatmap ─────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    cat_order = ["P", "LP", "VUS", "LB", "B", "Conflicting", "Other"]
    # Use top 15 codes for heatmap
    top15_codes = [c[0] for c in sorted(all_freq.items(), key=lambda x: -x[1])[:15]]

    if top15_codes and code_by_cat:
        mat = np.zeros((len(top15_codes), len(cat_order)))
        for i, code in enumerate(top15_codes):
            for j, cat in enumerate(cat_order):
                mat[i, j] = code_by_cat.get(cat, {}).get(code, 0)

        # Log-scale normalisation for visibility
        mat_display = np.where(mat > 0, np.log10(mat + 1), 0)
        im = ax.imshow(mat_display, aspect="auto", cmap="YlOrRd",
                       interpolation="nearest")
        ax.set_yticks(range(len(top15_codes)))
        ax.set_yticklabels(top15_codes, fontsize=5.5)
        ax.set_xticks(range(len(cat_order)))
        ax.set_xticklabels([CATEGORY_LABELS.get(c, c) for c in cat_order],
                            rotation=45, ha="right", fontsize=5.5)

        # Annotate cells
        for i in range(len(top15_codes)):
            for j in range(len(cat_order)):
                v = mat[i, j]
                if v > 0:
                    color = "white" if mat_display[i, j] > mat_display.max() * 0.6 else "black"
                    ax.text(j, i, fmt(v), ha="center", va="center",
                            fontsize=3.8, color=color)

        cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.03)
        cbar.ax.tick_params(labelsize=5)
        cbar.set_label("log₁₀(count + 1)", fontsize=6)

    ax.set_title("Criteria by classification")
    panel_label(ax, "d")

    _save(fig, out, "figure9_acmg_criteria_analysis")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 10 — Disease Domain Analysis (NEW)
# ═══════════════════════════════════════════════════════════════════════

DOMAIN_COLORS = {
    "Hereditary Cancer":  NAT["dark_red"],
    "Cardiovascular":     NAT["blue"],
    "Neurological":       NAT["violet"],
    "Renal":              NAT["aqua"],
    "Ophthalmologic":     NAT["yellow"],
    "Other/Unclassified": NAT["grey"],
}

def figure10_disease_domain(domain_data, out):
    """Eight-panel disease domain subgroup analysis.
    Integrates domain breakdown, cancer vs non-cancer, and rare vs common
    into a single comprehensive figure.
    """
    if not domain_data:
        return

    by_domain = domain_data.get("by_domain", {})
    cancer_vs = domain_data.get("cancer_vs_noncancer", {})
    rare_vs = domain_data.get("rare_vs_common", {})
    if not by_domain:
        return

    domain_order = ["Hereditary Cancer", "Cardiovascular", "Neurological",
                    "Renal", "Ophthalmologic", "Other/Unclassified"]
    domains_present = [d for d in domain_order if d in by_domain]

    fig = plt.figure(figsize=(7.08, 9.5))
    gs = GridSpec(4, 2, figure=fig, hspace=0.65, wspace=0.48)

    # ── (a) Reclassification rate by domain + cancer vs non-cancer ────
    ax = fig.add_subplot(gs[0, 0])
    d_labels = [d for d in domains_present if d != "Other/Unclassified"]
    d_rates = [by_domain[d].get("change_rate_pct", 0) for d in d_labels]
    d_cols = [DOMAIN_COLORS.get(d, NAT["grey"]) for d in d_labels]

    y_pos = np.arange(len(d_labels))
    ax.barh(y_pos, d_rates, color=d_cols, height=0.6, edgecolor="white",
            linewidth=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(d_labels, fontsize=6)
    ax.set_xlabel("Reclassification rate (%)")
    ax.set_title("Reclassification rate by disease domain")
    for i, v in enumerate(d_rates):
        ax.text(v + 0.2, i, f"{v:.1f}%", va="center", fontsize=5, color=NAT["grey"])
    panel_label(ax, "a", x=-0.28)

    # ── (b) Cancer vs non-cancer reclassification rate comparison ─────
    ax = fig.add_subplot(gs[0, 1])
    if cancer_vs:
        groups = ["cancer", "noncancer"]
        group_labels = ["Cancer\ngenes", "Non-cancer\ngenes"]
        rates = [cancer_vs.get(g, {}).get("change_rate_pct", 0) for g in groups]
        totals = [cancer_vs.get(g, {}).get("total", 0) for g in groups]
        bar_colors = [NAT["dark_red"], NAT["blue"]]

        x_pos = np.arange(len(groups))
        bars = ax.bar(x_pos, rates, 0.5, color=bar_colors, edgecolor="white",
                      linewidth=0.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(group_labels, fontsize=6)
        ax.set_ylabel("Reclassification rate (%)")
        ax.set_title("Cancer vs non-cancer reclassification")
        for i, (v, t) in enumerate(zip(rates, totals)):
            ax.text(i, v + 0.3, f"{v:.1f}%\n(n={t:,})", ha="center",
                    va="bottom", fontsize=5.5, fontweight="bold", color=NAT["grey"])
        ax.set_ylim(0, max(rates) * 1.35)

        # Add rare vs common as inset if available
        if rare_vs:
            rare_rate = rare_vs.get("rare", {}).get("change_rate_pct", 0)
            common_rate = rare_vs.get("common", {}).get("change_rate_pct", 0)
            rare_n = rare_vs.get("rare", {}).get("total", 0)
            common_n = rare_vs.get("common", {}).get("total", 0)
            # Add text annotation
            ax.text(0.98, 0.95,
                    f"Lower-prevalence: {rare_rate:.1f}% (n={rare_n:,})\n"
                    f"Higher-prevalence: {common_rate:.1f}% (n={common_n:,})",
                    transform=ax.transAxes, fontsize=5, ha="right", va="top",
                    bbox=dict(boxstyle="round,pad=0.3", fc="#f5f5f5",
                              ec=NAT["lt_grey"], lw=0.5))
    panel_label(ax, "b")

    # ── (c) Classification distribution: Cancer vs Non-Cancer ─────────
    ax = fig.add_subplot(gs[1, 0])
    if cancer_vs:
        groups = ["cancer", "noncancer"]
        group_labels = ["Cancer genes", "Non-cancer genes"]
        x_pos = np.arange(len(groups))
        cats_to_show = ["P", "LP", "VUS", "LB", "B"]
        cat_colors = [CATEGORY_COLORS[c] for c in cats_to_show]

        bottoms = np.zeros(len(groups))
        for ci, cat in enumerate(cats_to_show):
            vals = []
            for g in groups:
                dist = cancer_vs.get(g, {}).get("classification_dist", {})
                total = cancer_vs.get(g, {}).get("total", 1)
                vals.append(dist.get(cat, 0) / max(total, 1) * 100)
            ax.bar(x_pos, vals, 0.7, bottom=bottoms,
                   color=cat_colors[ci], edgecolor="white", linewidth=0.3,
                   label=CATEGORY_LABELS.get(cat, cat))
            bottoms += np.array(vals)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(group_labels, fontsize=6)
        ax.set_ylabel("Proportion (%)")
        ax.set_ylim(0, 105)
        ax.set_title("Classification: cancer vs non-cancer", pad=14)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18),
                  ncol=5, fontsize=4.5, frameon=False)
    panel_label(ax, "c", x=-0.28)

    # ── (d) Classification distribution: Rare vs Common ───────────────
    ax = fig.add_subplot(gs[1, 1])
    if rare_vs:
        groups = ["rare", "common"]
        group_labels = ["Lower-prevalence\nMendelian", "Higher-prevalence\nMendelian"]
        x_pos = np.arange(len(groups))
        cats_to_show = ["P", "LP", "VUS", "LB", "B"]
        cat_colors = [CATEGORY_COLORS[c] for c in cats_to_show]

        bottoms = np.zeros(len(groups))
        for ci, cat in enumerate(cats_to_show):
            vals = []
            for g in groups:
                dist = rare_vs.get(g, {}).get("classification_dist", {})
                total = rare_vs.get(g, {}).get("total", 1)
                vals.append(dist.get(cat, 0) / max(total, 1) * 100)
            ax.bar(x_pos, vals, 0.7, bottom=bottoms,
                   color=cat_colors[ci], edgecolor="white", linewidth=0.3,
                   label=CATEGORY_LABELS.get(cat, cat))
            bottoms += np.array(vals)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(group_labels, fontsize=5.5)
        ax.set_ylabel("Proportion (%)")
        ax.set_ylim(0, 105)
        ax.set_title("Classification: rare vs higher-prevalence Mendelian", pad=14)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18),
                  ncol=5, fontsize=4.5, frameon=False)
    panel_label(ax, "d")

    # ── (e) P↔B reversal rates by domain ──────────────────────────────
    ax = fig.add_subplot(gs[2, 0])
    d_labels_c = [d for d in domains_present if d != "Other/Unclassified"]
    p2b_rates = []
    b2p_rates = []
    for d in d_labels_c:
        total = max(by_domain[d].get("total", 1), 1)
        p2b_rates.append(by_domain[d].get("p_to_b", 0) / total * 1000)
        b2p_rates.append(by_domain[d].get("b_to_p", 0) / total * 1000)

    y_pos = np.arange(len(d_labels_c))
    bar_h = 0.35
    ax.barh(y_pos - bar_h / 2, p2b_rates, bar_h, color=NAT["blue"],
            edgecolor="white", linewidth=0.3, label="P-side \u2192 B-side")
    ax.barh(y_pos + bar_h / 2, b2p_rates, bar_h, color=NAT["dark_red"],
            edgecolor="white", linewidth=0.3, label="B-side \u2192 P-side")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(d_labels_c, fontsize=6)
    ax.set_xlabel("Reversals per 1,000 variants")
    ax.set_title("P \u2194 B reversal rates by domain")
    ax.legend(fontsize=5, loc="lower right",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "e", x=-0.28)

    # ── (f) Top 20 genes by variant count ─────────────────────────────
    ax = fig.add_subplot(gs[2, 1])
    all_genes = []
    for d in domains_present:
        for g_entry in by_domain[d].get("top_genes", []):
            gene = g_entry.get("gene", "") if isinstance(g_entry, dict) else g_entry[0]
            count = g_entry.get("count", 0) if isinstance(g_entry, dict) else g_entry[1]
            if count >= 500:
                all_genes.append((gene, count, d))

    all_genes.sort(key=lambda x: -x[1])
    top20 = all_genes[:20]
    if top20:
        g_labels = [g[0] for g in reversed(top20)]
        g_vals = [g[1] for g in reversed(top20)]
        g_cols = [DOMAIN_COLORS.get(g[2], NAT["grey"]) for g in reversed(top20)]

        y_pos = np.arange(len(g_labels))
        ax.barh(y_pos, g_vals, color=g_cols, height=0.7, edgecolor="white",
                linewidth=0.2)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(g_labels, fontsize=5)
        ax.set_xlabel("Total variants")
        ax.set_title("Top 20 genes by variant count")
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
        domain_patches = [mpatches.Patch(color=DOMAIN_COLORS.get(d, NAT["grey"]), label=d)
                          for d in domains_present if d != "Other/Unclassified"]
        ax.legend(handles=domain_patches, fontsize=4.5, loc="lower right",
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "f", x=-0.22)

    # ── (g) Evidence depth by domain ──────────────────────────────────
    ax = fig.add_subplot(gs[3, 0])
    d_labels_e = [d for d in domains_present if d != "Other/Unclassified"]
    mean_pmids = [by_domain[d].get("mean_pmids", 0) for d in d_labels_e]
    d_cols_e = [DOMAIN_COLORS.get(d, NAT["grey"]) for d in d_labels_e]

    y_pos = np.arange(len(d_labels_e))
    ax.barh(y_pos, mean_pmids, color=d_cols_e, height=0.6, edgecolor="white",
            linewidth=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(d_labels_e, fontsize=6)
    ax.set_xlabel("Mean PubMed citations per variant")
    ax.set_title("Literature evidence depth by domain")
    for i, v in enumerate(mean_pmids):
        ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=5, color=NAT["grey"])
    panel_label(ax, "g", x=-0.28)

    # ── (h) Evidence depth: rare vs common ────────────────────────────
    ax = fig.add_subplot(gs[3, 1])
    if rare_vs and cancer_vs:
        subgroups = ["cancer", "noncancer", "rare", "common"]
        sub_labels = ["Cancer", "Non-cancer", "Lower-prev.", "Higher-prev."]
        sub_colors = [NAT["dark_red"], NAT["blue"], NAT["cyan"], NAT["salmon"]]
        subs_vals = []
        pmid_vals = []
        for sg in subgroups:
            source = cancer_vs if sg in ("cancer", "noncancer") else rare_vs
            subs_vals.append(source.get(sg, {}).get("mean_submissions", 0))
            pmid_vals.append(source.get(sg, {}).get("mean_pmids", 0))

        x_pos = np.arange(len(subgroups))
        bar_w = 0.35
        bars1 = ax.bar(x_pos - bar_w / 2, subs_vals, bar_w, color=sub_colors,
                        alpha=0.7, edgecolor="white", linewidth=0.3,
                        label="Mean submissions")
        bars2 = ax.bar(x_pos + bar_w / 2, pmid_vals, bar_w, color=sub_colors,
                        alpha=0.4, edgecolor="white", linewidth=0.3,
                        hatch="//", label="Mean PubMed IDs")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(sub_labels, fontsize=5.5)
        ax.set_ylabel("Mean per variant")
        ax.set_title("Evidence depth by subgroup", pad=18)
        for i, (s, p) in enumerate(zip(subs_vals, pmid_vals)):
            ax.text(i - bar_w / 2, s + 0.03, f"{s:.2f}", ha="center",
                    fontsize=4.5, va="bottom")
            ax.text(i + bar_w / 2, p + 0.03, f"{p:.2f}", ha="center",
                    fontsize=4.5, va="bottom")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12),
                  ncol=2, fontsize=4.5, frameon=False)
    panel_label(ax, "h")

    _save(fig, out, "figure10_disease_domain_analysis")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 11 — Reclassification Driver Analysis (schema vs knowledge)
# ═══════════════════════════════════════════════════════════════════════

DRIVER_COLORS = {
    "likely_schema_driven":    NAT["dark_red"],
    "likely_knowledge_driven": NAT["lime"],
    "transition_ambiguous":    NAT["yellow"],
    "new_submitter":           NAT["blue"],
    "reinterpretation":        NAT["aqua"],
    "conflict_emerged":        NAT["salmon"],
    "conflict_resolved":       NAT["purple"],
    "batch_event":             NAT["grey"],
}

DRIVER_LABELS = {
    "likely_schema_driven":    "Schema-driven",
    "likely_knowledge_driven": "Evidence-based",
    "transition_ambiguous":    "Ambiguous (transition)",
    "new_submitter":           "New submitter",
    "reinterpretation":        "Reinterpretation",
    "conflict_emerged":        "Conflict emerged",
    "conflict_resolved":       "Conflict resolved",
    "batch_event":             "Batch event",
}

# Ordered keys for stacking in figures
_DRIVER_STACK_ORDER = [
    "likely_schema_driven", "likely_knowledge_driven", "transition_ambiguous",
    "new_submitter", "reinterpretation", "conflict_emerged", "conflict_resolved",
    "batch_event",
]


def figure11_reclassification_drivers(driver_data, out):
    """Four-panel figure: reclassification driver analysis.

    v2: Transition window uses schema/knowledge labels; outside transition
    uses observable-mechanism labels (new_submitter, reinterpretation,
    conflict_emerged, conflict_resolved).
    """
    if not driver_data:
        return

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 8.5))
    fig.subplots_adjust(hspace=1.30, wspace=0.55, top=0.94, bottom=0.06,
                        left=0.09, right=0.96)

    summary = driver_data.get("driver_summary", {})

    # ── (a) Overall mechanism breakdown (donut) ──────────────────────
    ax = axes[0, 0]
    # Exclude transition_ambiguous (n=105, <0.02%) — visually redundant
    keys_ordered = [k for k in _DRIVER_STACK_ORDER
                    if k != "transition_ambiguous"
                    and summary.get(k, {}).get("count", 0) > 0]
    counts = [summary[k]["count"] for k in keys_ordered]
    colors = [DRIVER_COLORS[k] for k in keys_ordered]
    labels = [DRIVER_LABELS[k] for k in keys_ordered]

    if counts:
        wedges, texts, autotexts = ax.pie(
            counts, labels=None,
            autopct=lambda p: f"{p:.1f}%" if p > 1.0 else "",
            colors=colors, startangle=90, pctdistance=0.78,
            wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.0),
        )
        for at in autotexts:
            at.set_fontsize(5)
            at.set_fontweight("bold")
        ax.legend(labels, fontsize=5, loc="upper center",
                  bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False,
                  columnspacing=0.8, handletextpad=0.4)
        total_events = sum(counts)
        ax.text(0, 0, f"{total_events:,}\nevents", ha="center", va="center",
                fontsize=7, fontweight="bold", color=NAT["grey"])
    ax.set_title("Reclassification mechanism breakdown", fontsize=9, pad=14)
    panel_label(ax, "a", x=-0.05, y=1.18)

    # ── (b) Mechanism breakdown by era (stacked bar) ─────────────────
    ax = axes[0, 1]
    era_summary = driver_data.get("era_summary", {})
    era_order = ["RCV (2012-2019)", "Transition (2020-01 to 2020-06)",
                 "VCV-stable (2020-07+)"]
    eras_present = [e for e in era_order if e in era_summary]

    if eras_present:
        era_short = {"RCV (2012-2019)": "RCV\n(2012-19)",
                     "Transition (2020-01 to 2020-06)": "Transition\n(2020 H1)",
                     "VCV-stable (2020-07+)": "VCV-stable\n(2020 H2+)"}
        x_pos = np.arange(len(eras_present))
        bar_width = 0.6
        bottoms = np.zeros(len(eras_present))

        for dk in _DRIVER_STACK_ORDER:
            vals = []
            for e in eras_present:
                total_e = max(era_summary[e].get("total_events", 1), 1)
                vals.append(era_summary[e].get(dk, 0) / total_e * 100)
            if max(vals) > 0:
                ax.bar(x_pos, vals, bar_width, bottom=bottoms,
                       color=DRIVER_COLORS[dk], edgecolor="white", linewidth=0.3,
                       label=DRIVER_LABELS[dk])
                bottoms += np.array(vals)

        ax.set_xticks(x_pos)
        ax.set_xticklabels([era_short.get(e, e) for e in eras_present], fontsize=6)
        ax.set_ylabel("Proportion (%)", fontsize=7)
        ax.set_ylim(0, 105)
        ax.set_title("Mechanism breakdown by era", fontsize=9, pad=14)
        ax.legend(fontsize=4.5, loc="upper center", bbox_to_anchor=(0.5, -0.35),
                  ncol=2, frameon=False, columnspacing=0.8, handletextpad=0.4)
    panel_label(ax, "b", x=-0.14, y=1.18)

    # ── (c) Time series ──────────────────────────────────────────────
    ax = axes[1, 0]
    timeseries = driver_data.get("timeseries", [])
    if timeseries:
        ts_dates = [e["date"] for e in timeseries]
        ts_x = np.array([date_to_float(d) for d in ts_dates])

        totals = np.array([max(e.get("total_reclassifications", 1), 1)
                           for e in timeseries], dtype=float)
        bottom = np.zeros_like(ts_x)
        for dk in _DRIVER_STACK_ORDER:
            vals = np.array([e.get(dk, 0) for e in timeseries], dtype=float)
            if vals.sum() > 0:
                pcts = vals / totals * 100
                ax.fill_between(ts_x, bottom, bottom + pcts,
                                color=DRIVER_COLORS[dk], alpha=0.65,
                                linewidth=0.15, edgecolor="white",
                                label=DRIVER_LABELS[dk])
                bottom = bottom + pcts

        ax.set_ylim(0, 100)
        ax.set_ylabel("Proportion (%)", fontsize=7)
        ax.set_xlabel("Year", fontsize=7)
        ax.set_title("Reclassification mechanism over time", fontsize=9, pad=14)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(2))
        ax.axvline(x=date_to_float("2020-01"), color="black", linestyle="--",
                   linewidth=0.6, alpha=0.5)
        ax.text(date_to_float("2020-01"), 95, "RCV\u2192VCV", fontsize=5.5,
                ha="left", va="top", style="italic", alpha=0.6)
        ax.legend(fontsize=4.5, loc="upper center", bbox_to_anchor=(0.5, -0.28),
                  ncol=4, frameon=False, columnspacing=0.8, handletextpad=0.3)
    panel_label(ax, "c", x=-0.14, y=1.18)

    # ── (d) Top reclassification types by mechanism ──────────────────
    ax = axes[1, 1]
    top_changes = driver_data.get("top_changes_by_driver", {})
    # Collect top types across the main non-batch mechanisms
    mechanism_keys = ["new_submitter", "reinterpretation",
                      "conflict_emerged", "conflict_resolved"]
    all_top = []
    for mk in mechanism_keys:
        for entry in top_changes.get(mk, [])[:5]:
            all_top.append((mk, entry["change_type"], entry["count"]))
    # Sort by count, deduplicate change_types (keep highest)
    seen = set()
    unique_top = []
    for mk, ct, c in sorted(all_top, key=lambda x: -x[2]):
        if ct not in seen and len(unique_top) < 10:
            unique_top.append((mk, ct, c))
            seen.add(ct)
    unique_top.reverse()  # smallest at bottom for barh

    if unique_top:
        t_labels = [ct for _, ct, _ in unique_top]
        t_vals = [c for _, _, c in unique_top]
        t_colors = [DRIVER_COLORS[mk] for mk, _, _ in unique_top]
        y_pos = np.arange(len(t_labels))
        ax.barh(y_pos, t_vals, color=t_colors, height=0.7, alpha=0.75,
                edgecolor="white", linewidth=0.2)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(t_labels, fontsize=6)
        ax.set_xlabel("Event count", fontsize=7)
        ax.set_title("Top reclassification types by mechanism", fontsize=9, pad=14)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
        ax.tick_params(axis='x', labelsize=6)
        # Legend for mechanism colors
        handles = [mpatches.Patch(color=DRIVER_COLORS[mk], label=DRIVER_LABELS[mk])
                   for mk in mechanism_keys if any(m == mk for m, _, _ in unique_top)]
        ax.legend(handles=handles, fontsize=5, loc="lower right",
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.85,
                  handletextpad=0.4, borderpad=0.3)
    panel_label(ax, "d", x=-0.16, y=1.18)

    _save(fig, out, "figure11_reclassification_drivers")


# ═══════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY — Rare vs Common Disease Analysis
# ═══════════════════════════════════════════════════════════════════════

def figure_supp_rare_vs_common(domain_data, out):
    """Supplementary figure: rare vs common disease variant comparison."""
    if not domain_data:
        return
    rare_common = domain_data.get("rare_vs_common", {})
    if not rare_common:
        return

    fig, axes = plt.subplots(1, 3, figsize=(7.08, 2.6))
    fig.subplots_adjust(wspace=0.50)

    groups = ["rare", "common"]
    group_labels = ["Lower-prevalence\nMendelian", "Higher-prevalence\nMendelian"]
    group_colors = [NAT["blue"], NAT["yellow"]]

    # ── (a) Variant count & change rate ──────────────────────────────
    ax = axes[0]
    totals = [rare_common.get(g, {}).get("total", 0) for g in groups]
    change_rates = [rare_common.get(g, {}).get("change_rate_pct", 0) for g in groups]

    x_pos = np.arange(len(groups))
    bars = ax.bar(x_pos, totals, color=group_colors, width=0.5,
                  edgecolor="white", linewidth=0.3)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_labels, fontsize=5.5)
    ax.set_ylabel("Total variants")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: fmt(v)))
    ax.set_title("Variant count by disease prevalence")
    for i, (t, cr) in enumerate(zip(totals, change_rates)):
        ax.text(i, t + t * 0.02, f"{fmt(t)}\n({cr:.1f}% changed)",
                ha="center", va="bottom", fontsize=5, color=NAT["grey"])
    panel_label(ax, "a")

    # ── (b) Classification distribution ──────────────────────────────
    ax = axes[1]
    cats_to_show = ["P", "LP", "VUS", "LB", "B"]
    bar_width = 0.35
    x_pos = np.arange(len(groups))
    bottoms = np.zeros(len(groups))
    for cat in cats_to_show:
        vals = []
        for g in groups:
            dist = rare_common.get(g, {}).get("classification_dist", {})
            total = max(rare_common.get(g, {}).get("total", 1), 1)
            vals.append(dist.get(cat, 0) / total * 100)
        ax.bar(x_pos, vals, bar_width * 2, bottom=bottoms,
               color=CATEGORY_COLORS[cat], edgecolor="white", linewidth=0.3,
               label=CATEGORY_LABELS.get(cat, cat))
        bottoms += np.array(vals)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_labels, fontsize=6)
    ax.set_ylabel("Proportion (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Classification: by disease prevalence")
    ax.legend(fontsize=5, loc="upper right", ncol=1,
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "b")

    # ── (c) Evidence depth comparison ────────────────────────────────
    ax = axes[2]
    metrics = ["mean_submissions", "mean_pmids"]
    metric_labels = ["Mean submissions", "Mean PubMed IDs"]
    x_pos = np.arange(len(metrics))
    bar_w = 0.3
    for gi, (g, gl, gc) in enumerate(zip(groups, group_labels, group_colors)):
        vals = [rare_common.get(g, {}).get(m, 0) for m in metrics]
        ax.bar(x_pos + (gi - 0.5) * bar_w, vals, bar_w, color=gc,
               edgecolor="white", linewidth=0.3, label=gl)
        for xi, v in enumerate(vals):
            ax.text(x_pos[xi] + (gi - 0.5) * bar_w, v + 0.02,
                    f"{v:.2f}", ha="center", va="bottom", fontsize=5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(metric_labels, fontsize=6)
    ax.set_ylabel("Mean per variant")
    ax.set_title("Evidence depth: by disease prevalence")
    ax.legend(fontsize=5,
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax, "c")

    _save(fig, out, "figure_supp_rare_vs_common")


# ═══════════════════════════════════════════════════════════════════════
# METADATA STATISTICS TABLE (for Methods section)
# ═══════════════════════════════════════════════════════════════════════

def generate_metadata_table(summary_agg, hist_agg, stats, xml_data,
                            acmg_stats, out_dir):
    """Generate a metadata statistics table for the manuscript Methods section.

    Outputs plain-text, JSON, and LaTeX versions to out_dir.
    """
    from datetime import datetime

    lines = []
    data = {}

    def _s(header):
        lines.append("")
        lines.append("=" * 70)
        lines.append(header)
        lines.append("=" * 70)

    # ── 1. Dataset Overview ──────────────────────────────────────────
    _s("1. Dataset Overview")

    total = summary_agg.get("total", 0)
    changed = summary_agg.get("changed", 0)
    stable = summary_agg.get("stable", 0)
    oscillating = summary_agg.get("oscillating", 0)

    n_snapshots = len(hist_agg) if hist_agg else 0
    if hist_agg:
        dates_sorted = sorted(hist_agg.keys())
        date_range = f"{dates_sorted[0]} to {dates_sorted[-1]}"
    else:
        date_range = "N/A"

    # Format split from XML data
    n_rcv = n_vcv = 0
    if xml_data:
        for s in xml_data.get("snapshots", []):
            if s.get("format") == "rcv":
                n_rcv += 1
            elif s.get("format") == "vcv":
                n_vcv += 1

    lines.append(f"  Total unique variants:         {total:>12,}")
    lines.append(f"  Longitudinal snapshots:        {n_snapshots:>12,}")
    lines.append(f"  Date range:                    {date_range:>12}")
    lines.append(f"  RCV-format snapshots:          {n_rcv:>12,}")
    lines.append(f"  VCV-format snapshots:          {n_vcv:>12,}")

    data["dataset_overview"] = {
        "total_variants": total,
        "snapshots": n_snapshots,
        "date_range": date_range,
        "rcv_snapshots": n_rcv,
        "vcv_snapshots": n_vcv,
    }

    # ── 2. Current Classification Distribution ───────────────────────
    _s("2. Current Classification Distribution")

    cls_dist = stats.get("current_classification_distribution", {})
    if not cls_dist:
        # Fallback: reconstruct from transition matrix diagonal
        for key, count in summary_agg.get("transition_matrix", {}).items():
            parts = key.split("|")
            if len(parts) == 2:
                cur = parts[1]
                cls_dist[cur] = cls_dist.get(cur, 0) + count

    cls_data = {}
    for cat in CATEGORY_ORDER:
        count = cls_dist.get(cat, 0)
        pct = count / max(total, 1) * 100
        lines.append(f"  {cat:<20s}  {count:>10,}  ({pct:>5.1f}%)")
        cls_data[cat] = {"count": count, "pct": round(pct, 2)}

    data["classification_distribution"] = cls_data

    # ── 3. Reclassification Dynamics ─────────────────────────────────
    _s("3. Reclassification Dynamics")

    change_pct = changed / max(total, 1) * 100
    osc_pct = oscillating / max(total, 1) * 100
    p2b_total = sum(summary_agg.get("p_to_b", {}).values())
    b2p_total = sum(summary_agg.get("b_to_p", {}).values())

    lines.append(f"  Variants with reclassification: {changed:>10,}  ({change_pct:.1f}%)")
    lines.append(f"  Stable variants:                {stable:>10,}  ({100 - change_pct:.1f}%)")
    lines.append(f"  Oscillating variants:           {oscillating:>10,}  ({osc_pct:.1f}%)")
    lines.append(f"  P/LP → B/LB reversals:          {p2b_total:>10,}")
    lines.append(f"  B/LB → P/LP reversals:          {b2p_total:>10,}")

    # Top reclassification patterns
    change_dist = stats.get("change_type_distribution", {})
    if isinstance(change_dist, dict):
        top5 = sorted(change_dist.items(), key=lambda x: -x[1])[:5]
    else:
        top5 = change_dist[:5]
    lines.append("  Top reclassification patterns:")
    for pattern, count in top5:
        lines.append(f"    {pattern:<30s}  {count:>8,}")

    data["reclassification"] = {
        "changed": changed, "stable": stable, "oscillating": oscillating,
        "change_pct": round(change_pct, 2),
        "p_to_b_total": p2b_total, "b_to_p_total": b2p_total,
        "top_patterns": [{"pattern": p, "count": c} for p, c in top5],
    }

    # ── 4. Submission & Evidence Coverage ────────────────────────────
    _s("4. Submission & Evidence Coverage")

    sub_stats = stats.get("submissions_per_variant", {})
    mean_sub = sub_stats.get("mean", 0)
    median_sub = sub_stats.get("median", 0)
    max_sub = sub_stats.get("max", 0)

    # Multi-submitter % from last snapshot
    if hist_agg:
        last_date = sorted(hist_agg.keys())[-1]
        pct_multi = hist_agg[last_date].get("pct_multi_sub", 0)
    else:
        pct_multi = 0

    pmid_bins = summary_agg.get("pmid_bins", {})
    has_pmid = total - pmid_bins.get("0", 0)
    pmid_pct = has_pmid / max(total, 1) * 100

    acmg_bins = summary_agg.get("acmg_bins", {})
    has_acmg = total - acmg_bins.get("0", 0)
    acmg_pct = has_acmg / max(total, 1) * 100

    lines.append(f"  Mean submissions per variant:   {mean_sub:>12.2f}")
    lines.append(f"  Median submissions:             {median_sub:>12}")
    lines.append(f"  Max submissions:                {max_sub:>12,}")
    lines.append(f"  Variants with >1 submission:    {pct_multi:>11.1f}%")
    lines.append(f"  Variants with ≥1 PubMed ID:     {has_pmid:>10,}  ({pmid_pct:.1f}%)")
    lines.append(f"  Variants with ≥1 ACMG code:     {has_acmg:>10,}  ({acmg_pct:.1f}%)")

    data["evidence_coverage"] = {
        "mean_submissions": round(mean_sub, 2),
        "median_submissions": median_sub,
        "max_submissions": max_sub,
        "pct_multi_submission": round(pct_multi, 1),
        "has_pmid": has_pmid, "pmid_pct": round(pmid_pct, 1),
        "has_acmg": has_acmg, "acmg_pct": round(acmg_pct, 1),
    }

    # ── 5. Field Completeness ────────────────────────────────────────
    _s("5. Field Completeness")

    fc = summary_agg.get("field_completeness", {})
    fc_data = {}
    for field in sorted(fc.keys(), key=lambda k: -fc[k]):
        lines.append(f"  {field:<25s}  {fc[field]:>6.1f}%")
        fc_data[field] = fc[field]

    data["field_completeness"] = fc_data

    # ── 6. Schema Evolution ──────────────────────────────────────────
    _s("6. Schema Evolution")

    if xml_data:
        import re as _re
        valid_snaps = [s for s in xml_data.get("snapshots", [])
                       if _re.match(r"^\d{4}-\d{2}", s["snapshot_date"])]
        if valid_snaps:
            first_s = valid_snaps[0]
            last_s = valid_snaps[-1]
            lines.append(f"  Schema snapshots analyzed:      {len(valid_snaps):>12,}")
            lines.append(f"  First snapshot elements:        {first_s['num_unique_elements']:>12,}")
            lines.append(f"  Last snapshot elements:         {last_s['num_unique_elements']:>12,}")
            lines.append(f"  First snapshot attributes:      {first_s['total_unique_attributes']:>12,}")
            lines.append(f"  Last snapshot attributes:       {last_s['total_unique_attributes']:>12,}")
            lines.append(f"  Max tree depth (latest):        {last_s['max_depth']:>12}")

            # XRef databases in latest
            xrefs = last_s.get("xref_databases", {})
            if xrefs:
                lines.append("  Cross-database references (latest snapshot):")
                for db, count in sorted(xrefs.items(), key=lambda x: -x[1])[:8]:
                    lines.append(f"    {db:<20s}  {count:>8,}")

            data["schema_evolution"] = {
                "snapshots_analyzed": len(valid_snaps),
                "first_elements": first_s["num_unique_elements"],
                "last_elements": last_s["num_unique_elements"],
                "first_attributes": first_s["total_unique_attributes"],
                "last_attributes": last_s["total_unique_attributes"],
                "max_depth": last_s["max_depth"],
                "top_xrefs": dict(sorted(xrefs.items(), key=lambda x: -x[1])[:8]),
            }
    else:
        lines.append("  (XML structure analysis data not available)")

    # ── Write output files ───────────────────────────────────────────
    txt_path = out_dir / "metadata_stats.txt"
    json_path = out_dir / "metadata_stats.json"
    tex_path = out_dir / "metadata_stats.tex"

    # Plain text
    txt_content = "\n".join(lines) + "\n"
    with open(txt_path, "w") as f:
        f.write(txt_content)

    # JSON
    data["generated_at"] = datetime.now().isoformat()
    with open(json_path, "w") as f:
        json.dump(data, f, indent=2)

    # LaTeX table fragment
    tex_lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Dataset characteristics and summary statistics}",
        r"\label{tab:metadata}",
        r"\begin{tabular}{lr}",
        r"\toprule",
        r"\textbf{Metric} & \textbf{Value} \\",
        r"\midrule",
        f"Total unique variants & {total:,} \\\\",
        f"Longitudinal snapshots & {n_snapshots:,} \\\\",
        f"Date range & {date_range} \\\\",
        f"Variants reclassified & {changed:,} ({change_pct:.1f}\\%) \\\\",
        f"Oscillating variants & {oscillating:,} ({osc_pct:.1f}\\%) \\\\",
        f"P/LP $\\to$ B/LB reversals & {p2b_total:,} \\\\",
        f"B/LB $\\to$ P/LP reversals & {b2p_total:,} \\\\",
        f"Mean submissions/variant & {mean_sub:.2f} \\\\",
        f"Variants with $\\geq$1 PMID & {has_pmid:,} ({pmid_pct:.1f}\\%) \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    with open(tex_path, "w") as f:
        f.write("\n".join(tex_lines) + "\n")

    print(f"  Saved metadata_stats.txt / .json / .tex")


# ═══════════════════════════════════════════════════════════════════════
# SCHEMA CHANGE TAXONOMY HELPER (for Figure 5 revision)
# ═══════════════════════════════════════════════════════════════════════

def _classify_schema_changes(xml_data):
    """Classify each release's schema changes as additive, deprecation, or mixed.

    Returns a list of dicts with per-release taxonomy and running totals.
    """
    import re as _re
    lifecycle = [lc for lc in xml_data.get("element_lifecycle", [])
                 if _re.match(r"^\d{4}-\d{2}", lc["snapshot_date"])]

    results = []
    cum_new = 0
    cum_dep = 0

    for lc in lifecycle:
        n_new = lc.get("new_elements", 0)
        n_dep = lc.get("deprecated_elements", 0)
        depth_change = abs(lc.get("depth_change", 0))

        if n_new == 0 and n_dep == 0:
            change_type = "none"
        elif n_new > 0 and n_dep == 0:
            change_type = "additive"
        elif n_dep > 0 and n_new == 0:
            change_type = "deprecation"
        else:
            # Both new and deprecated
            if depth_change >= 2:
                change_type = "restructuring"
            else:
                change_type = "mixed"

        cum_new += n_new
        cum_dep += n_dep

        results.append({
            "date": lc["snapshot_date"],
            "new_elements": n_new,
            "deprecated_elements": n_dep,
            "change_type": change_type,
            "cumulative_new": cum_new,
            "cumulative_deprecated": cum_dep,
        })

    return results


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 12 — Time-to-Reclassification (Reviewer Comment #16)
# ═══════════════════════════════════════════════════════════════════════

def figure12_time_to_reclassification(comment_data_dir, out):
    """Two-panel figure: distribution + clinically important transition breakdown."""
    ttr_path = Path(comment_data_dir) / "time_to_reclassification.json"
    hist_path = Path(comment_data_dir) / "time_to_reclassification_histogram.json"
    if not ttr_path.exists() or not hist_path.exists():
        print("  [skip] time_to_reclassification data not found")
        return

    with open(ttr_path) as f:
        ttr = json.load(f)
    with open(hist_path) as f:
        hist_data = json.load(f)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.08, 3.0),
                                      gridspec_kw={"wspace": 0.42})

    # ── (a) Overall distribution histogram ────────────────────────
    bin_edges = np.array(hist_data["bin_edges_years"])
    counts = np.array(hist_data["counts"])
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]

    ax_a.bar(bin_centers, counts, width=bin_width * 0.92,
             color=NAT["blue"], edgecolor="white", linewidth=0.3, alpha=0.85)

    # Mark median and IQR
    overall = ttr["overall"]
    median_y = overall["median_years"]
    q25 = overall["q25_years"]
    q75 = overall["q75_years"]
    ymax = max(counts) * 1.12

    ax_a.axvline(median_y, color=NAT["dark_red"], linewidth=1.5,
                 linestyle="-", zorder=5)
    ax_a.axvspan(q25, q75, color=NAT["salmon"], alpha=0.15, zorder=0)
    ax_a.annotate(f"Median: {median_y:.2f} yr\n"
                  f"IQR: {q25:.1f}–{q75:.1f} yr",
                  xy=(median_y, ymax * 0.92),
                  xytext=(median_y + 2.5, ymax * 0.92),
                  fontsize=6, fontweight="bold", color=NAT["dark_red"],
                  ha="left", va="top",
                  arrowprops=dict(arrowstyle="-|>", color=NAT["dark_red"],
                                  lw=0.8))
    # Percentage annotations
    ax_a.text(0.97, 0.95,
              f"Within 1 yr: {overall['pct_within_1yr']:.1f}%\n"
              f"Within 2 yr: {overall['pct_within_2yr']:.1f}%\n"
              f"Within 5 yr: {overall['pct_within_5yr']:.1f}%",
              transform=ax_a.transAxes, fontsize=5.5, ha="right", va="top",
              bbox=dict(boxstyle="round,pad=0.3", fc="white",
                        ec=NAT["lt_grey"], lw=0.5))

    ax_a.set_xlabel("Years from first appearance to first reclassification")
    ax_a.set_ylabel("Number of variants")
    ax_a.set_title("Time-to-reclassification distribution")
    ax_a.set_xlim(0, 13)
    ax_a.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: fmt(v)))
    panel_label(ax_a, "a")

    # ── (b) Clinically important transitions ──────────────────────
    clin = ttr.get("clinically_important", {})
    if clin:
        labels_b = []
        medians_b = []
        iqr_lo = []
        iqr_hi = []
        n_vals = []
        colors_b = []

        transition_colors = {
            "VUS→P/LP": NAT["dark_red"],
            "VUS→B/LB": NAT["blue"],
            "VUS→Conflicting": NAT["violet"],
            "P/LP→B/LB (drastic)": NAT["salmon"],
            "B/LB→P/LP (drastic)": NAT["aqua"],
        }
        # Ordered for display
        display_order = [
            "VUS→Conflicting",
            "VUS→B/LB",
            "VUS→P/LP",
            "B/LB→P/LP (drastic)",
            "P/LP→B/LB (drastic)",
        ]

        for key in display_order:
            if key in clin:
                s = clin[key]
                labels_b.append(key)
                medians_b.append(s["median_years"])
                iqr_lo.append(s["q25_years"])
                iqr_hi.append(s["q75_years"])
                n_vals.append(s["n"])
                colors_b.append(transition_colors.get(key, NAT["grey"]))

        y_pos = np.arange(len(labels_b))
        # Plot median as point + IQR as horizontal errorbar
        for i in range(len(labels_b)):
            ax_b.barh(y_pos[i], medians_b[i], height=0.55,
                      color=colors_b[i], alpha=0.75,
                      edgecolor="white", linewidth=0.3)
            ax_b.plot([iqr_lo[i], iqr_hi[i]], [y_pos[i], y_pos[i]],
                      color=colors_b[i], linewidth=2.0, solid_capstyle="round")
            ax_b.plot(iqr_lo[i], y_pos[i], "|", color=colors_b[i],
                      markersize=6, markeredgewidth=1.5)
            ax_b.plot(iqr_hi[i], y_pos[i], "|", color=colors_b[i],
                      markersize=6, markeredgewidth=1.5)
            # n annotation
            ax_b.text(max(iqr_hi[i], medians_b[i]) + 0.15, y_pos[i],
                      f"n={n_vals[i]:,}", fontsize=5, va="center",
                      color=NAT["grey"])

        ax_b.set_yticks(y_pos)
        ax_b.set_yticklabels(labels_b, fontsize=6)
        ax_b.set_xlabel("Years (median + IQR)")
        ax_b.set_title("Time-to-reclassification by transition type")
        ax_b.set_xlim(0, max(max(iqr_hi) + 1.5, 6))

    panel_label(ax_b, "b")
    _save(fig, out, "figure12_time_to_reclassification")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 13 — Batch Event Decomposition (Reviewer Comment #17)
# ═══════════════════════════════════════════════════════════════════════

def figure13_batch_event_decomposition(comment_data_dir, out):
    """Two-panel figure: batch event characteristics and submission patterns."""
    batch_path = Path(comment_data_dir) / "batch_event_decomposition.json"
    if not batch_path.exists():
        print("  [skip] batch_event_decomposition.json not found")
        return

    with open(batch_path) as f:
        batch_data = json.load(f)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.08, 3.2),
                                      gridspec_kw={"wspace": 0.45})

    dates = sorted(batch_data.keys())
    x_pos = np.arange(len(dates))
    date_labels = [d.replace("20", "'", 1) for d in dates]

    # ── (a) Stacked bar: submission count changes ─────────────────
    subs_stable = [batch_data[d]["subs_stable_pct"] for d in dates]
    subs_inc = [batch_data[d]["subs_increased_pct"] for d in dates]
    subs_dec = [batch_data[d]["subs_decreased_pct"] for d in dates]

    bar_w = 0.55
    b1 = ax_a.bar(x_pos, subs_stable, bar_w,
                  color=NAT["blue"], edgecolor="white", linewidth=0.3,
                  label="Subs unchanged")
    b2 = ax_a.bar(x_pos, subs_inc, bar_w, bottom=subs_stable,
                  color=NAT["salmon"], edgecolor="white", linewidth=0.3,
                  label="Subs increased")
    bottoms = [s + i for s, i in zip(subs_stable, subs_inc)]
    b3 = ax_a.bar(x_pos, subs_dec, bar_w, bottom=bottoms,
                  color=NAT["yellow"], edgecolor="white", linewidth=0.3,
                  label="Subs decreased")

    # Event count annotation on top
    for i, d in enumerate(dates):
        n = batch_data[d]["total_events"]
        ax_a.text(x_pos[i], 102, f"n={n:,}", ha="center", fontsize=4.5,
                  color=NAT["grey"], fontweight="bold")
        # Source label at bottom
        source = batch_data[d]["likely_source"].split("(")[0].strip()
        short_source = {"ClinVar-wide re-processing": "ClinVar\nre-processing",
                        "Large-scale lab submission": "Lab\nsubmission",
                        "ClinVar re-aggregation": "ClinVar\nre-aggregation",
                        "Mixed": "Mixed"}.get(source, source[:12])
        ax_a.text(x_pos[i], -12, short_source, ha="center", fontsize=4,
                  color=NAT["grey"], style="italic")

    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(date_labels, fontsize=6)
    ax_a.set_ylabel("Percentage of events")
    ax_a.set_ylim(0, 110)
    ax_a.set_title("Submission count changes in batch events")
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
                ncol=3, fontsize=5, frameon=False)
    panel_label(ax_a, "a")

    # ── (b) Gene concentration + top change patterns ──────────────
    gene_counts = [batch_data[d]["unique_genes"] for d in dates]
    hhis = [batch_data[d]["gene_concentration_hhi"] for d in dates]

    color_bars = []
    for hhi in hhis:
        if hhi > 0.1:
            color_bars.append(NAT["dark_red"])   # concentrated
        elif hhi > 0.02:
            color_bars.append(NAT["yellow"])     # moderate
        else:
            color_bars.append(NAT["lime"])       # dispersed

    bars = ax_b.bar(x_pos, gene_counts, bar_w, color=color_bars,
                    edgecolor="white", linewidth=0.3)

    # HHI annotation
    for i, (gc, hhi) in enumerate(zip(gene_counts, hhis)):
        label = "concentrated" if hhi > 0.1 else \
                "moderate" if hhi > 0.02 else "dispersed"
        ax_b.text(x_pos[i], gc + 50, f"HHI={hhi:.3f}\n({label})",
                  ha="center", fontsize=4, color=NAT["grey"])

    # Top gene annotation (rotated to avoid overlap)
    for i, d in enumerate(dates):
        top_genes = batch_data[d]["top_10_genes"][:2]
        gene_str = ", ".join(g["gene"] for g in top_genes)
        ax_b.text(x_pos[i], -80, gene_str, ha="center", fontsize=4,
                  color=NAT["grey"], style="italic", rotation=0)

    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(date_labels, fontsize=6)
    ax_b.set_ylabel("Unique genes affected")
    ax_b.set_title("Gene concentration in batch events")

    # Legend for concentration categories
    from matplotlib.patches import Patch
    legend_patches = [
        Patch(color=NAT["lime"], label="Dispersed (HHI<0.02)"),
        Patch(color=NAT["yellow"], label="Moderate (0.02–0.1)"),
        Patch(color=NAT["dark_red"], label="Concentrated (HHI>0.1)"),
    ]
    ax_b.legend(handles=legend_patches, loc="upper center",
                bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=5,
                frameon=False)
    panel_label(ax_b, "b")

    plt.subplots_adjust(bottom=0.22)
    _save(fig, out, "figure13_batch_event_decomposition")


def figure14_submitter_landscape(comment_data_dir, out):
    """Four-panel figure: ClinVar submitter landscape and temporal dynamics.

    (a) Top 20 submitters by total submissions (horizontal bar)
    (b) Submission growth timeline by top submitters (stacked area)
    (c) Collection method distribution (donut)
    (d) Classification breakdown for top 10 submitters (stacked bar)
    """
    sub_overall_path = Path(comment_data_dir) / "submitter_overall.json"
    sub_temporal_path = Path(comment_data_dir) / "submitter_temporal.json"
    sub_stats_path = Path(comment_data_dir) / "submitter_stats.json"

    for p in [sub_overall_path, sub_temporal_path, sub_stats_path]:
        if not p.exists():
            print(f"  [skip] {p.name} not found — run clinvar_submitter_analysis.py first")
            return

    with open(sub_overall_path) as f:
        overall = json.load(f)
    with open(sub_temporal_path) as f:
        temporal = json.load(f)
    with open(sub_stats_path) as f:
        stats = json.load(f)

    fig = plt.figure(figsize=(7.08, 8.5))
    gs = GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.45,
                  height_ratios=[1.2, 1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, :])   # top 20 — full width
    ax_b = fig.add_subplot(gs[1, :])   # temporal — full width
    ax_c = fig.add_subplot(gs[2, 0])   # collection method donut
    ax_d = fig.add_subplot(gs[2, 1])   # classification breakdown

    # ── De-identification flag ──────────────────────────────────────
    # Set to True to replace lab names with anonymous labels (Laboratory A, B, ...)
    # Set to False to show real lab names
    DEIDENTIFY = True

    # Mapping: all known name forms -> anonymous label.
    # Includes full original names (from JSON), short display names, and truncated forms.
    # Adjust letter assignments to match your manuscript text.
    # Non-lab entities (OMIM, All of Us) are NOT in this map and keep real names.
    DEIDENTIFY_MAP = {
        # Laboratory A — Labcorp/Invitae (all name variants)
        "Labcorp/Invitae":       "Laboratory A",
        "Labcorp Genetics (formerly Invitae), Labcorp": "Laboratory A",
        "Invitae":               "Laboratory A",
        "Invitae,":              "Laboratory A",
        # Laboratory B — Ambry Genetics
        "Ambry Genetics":        "Laboratory B",
        # Laboratory C — GeneDx
        "GeneDx":                "Laboratory C",
        # Laboratory D — Illumina (all name variants)
        "Illumina Lab Services": "Laboratory D",
        "Illumina Laboratory Services, Illumina": "Laboratory D",
        "Illumina Clinical Services Laboratory,Illumina": "Laboratory D",
        # Laboratory E — Myriad Genetics
        "Myriad Genetics":       "Laboratory E",
        "Myriad Genetics, Inc.": "Laboratory E",
        # Laboratory F — PreventionGenetics (all name variants)
        "PreventionGenetics":    "Laboratory F",
        "PreventionGenetics, part of Exact Sciences": "Laboratory F",
        "PreventionGenetics,PreventionGenetics": "Laboratory F",
        # Laboratory G — CeGaT
        "CeGaT":                "Laboratory G",
        "CeGaT Center for Human Genetics Tuebingen": "Laboratory G",
        # Laboratory H — Breakthrough Genomics
        "Breakthrough Genomics": "Laboratory H",
        "Breakthrough Genomics, Breakthrough Genomics": "Laboratory H",
        # Laboratory I — Fulgent Genetics
        "Fulgent Genetics":      "Laboratory I",
        "Fulgent Genetics, Fulgent Genetics": "Laboratory I",
        # Laboratory J — LabCorp WHG
        "LabCorp WHG":           "Laboratory J",
        "Women's Health and Genetics/Laboratory Corporation of America, LabCorp": "Laboratory J",
        # Laboratory K — Color Health
        "Color Health":          "Laboratory K",
        "Color Diagnostics, LLC DBA Color Health": "Laboratory K",
        "Color Health, Inc": "Laboratory K",
        # Laboratory L — Mayo Clinic Labs
        "Mayo Clinic Labs":      "Laboratory L",
        "Mayo Clinic Laboratories, Mayo Clinic": "Laboratory L",
        # Laboratory M — Baylor Genetics
        "Baylor Genetics":       "Laboratory M",
        # Laboratory N — Natera
        "Natera":                "Laboratory N",
        "Natera, Inc.": "Laboratory N",
        # Laboratory O — Genome-Nilou Lab
        "Genome-Nilou Lab":      "Laboratory O",
        # Laboratory P — Revvity Omics
        "Revvity Omics":         "Laboratory P",
        "Revvity Omics, Revvity": "Laboratory P",
        # Laboratory Q — Mendelics
        "Mendelics":             "Laboratory Q",
        "Mendelics Analise Genomica": "Laboratory Q",
        # Laboratory R — Eurofins NTD
        "Eurofins NTD":          "Laboratory R",
        "Eurofins Ntd Llc (ga)": "Laboratory R",
        # Laboratory S — Quest Diagnostics
        "Quest Diagnostics":     "Laboratory S",
        "Quest Diagnostics Nichols Institute San Juan Capistrano": "Laboratory S",
        # Laboratory T — Counsyl (Myriad)
        "Counsyl (Myriad)":      "Laboratory T",
        "Counsyl":               "Laboratory T",
        # Laboratory U — EGL / Eurofins Clinical Diagnostics
        "EGL Genetic Diagnostics,Eurofins Clinical Diagnostics": "Laboratory U",
        "Eurofins Clinical Diagnostics": "Laboratory U",
        "Eurofins":              "Laboratory U",
        # Laboratory V — LMM (Laboratory for Molecular Medicine)
        "Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine": "Laboratory V",
        "Laboratory for Molecular Medicine, Mass General Brigham Personalized Medicine": "Laboratory V",
        # Laboratory W — Genetic Services Lab, U Chicago
        "Genetic Services Laboratory, University of Chicago": "Laboratory W",
        # Laboratory X — ARUP
        "ARUP Laboratories":     "Laboratory X",
        "ARUP Laboratories, Molecular Genetics and Genomics, ARUP Laboratories": "Laboratory X",
        # Laboratory Y — Athena
        "Athena Diagnostics":    "Laboratory Y",
        "Athena Diagnostics Inc": "Laboratory Y",
        # Laboratory Z — Blueprint
        "Blueprint Genetics":    "Laboratory Z",
        # Additional labs
        "Sinai Health System":   "Laboratory AA",
        "Department of Pathology and Laboratory Medicine, Sinai Health System": "Laboratory AA",
        "3billion":              "Laboratory AB",
        "University of Washington Department of Laboratory Medicine,University of Washington": "Laboratory AC",
        "Center for Pediatric Genomic Medicine,Children's Mercy Hospital and Clinics": "Laboratory AD",
    }

    # ── Short names for display ────────────────────────────────────
    def short_name(name, max_len=28):
        replacements = {
            "Labcorp Genetics (formerly Invitae), Labcorp": "Labcorp/Invitae",
            "PreventionGenetics, part of Exact Sciences": "PreventionGenetics",
            "CeGaT Center for Human Genetics Tuebingen": "CeGaT",
            "Illumina Laboratory Services, Illumina": "Illumina Lab Services",
            "Color Health, Inc": "Color Health",
            "Color Diagnostics, LLC DBA Color Health": "Color Health",
            "Mendelics Analise Genomica": "Mendelics",
            "ARUP Laboratories, Molecular Genetics and Genomics, ARUP Laboratories": "ARUP Laboratories",
            "Counsyl": "Counsyl (Myriad)",
            "Women's Health and Genetics/Laboratory Corporation of America, LabCorp": "LabCorp WHG",
            "Eurofins Clinical Diagnostics": "Eurofins",
            "Revvity Omics, Revvity": "Revvity Omics",
            "Centre for Mendelian Genomics, University Medical Centre Ljubljana": "UMC Ljubljana",
            "Institute of Human Genetics, University of Leipzig Medical Center": "U Leipzig",
            "Clinical Genetics Lab, RAMBAM Health Care Campus": "RAMBAM",
            "Institute of Medical Genetics and Applied Genomics, University of Tuebingen": "U Tübingen",
            "Blueprint Genetics": "Blueprint Genetics",
            "Breakthrough Genomics, Breakthrough Genomics": "Breakthrough Genomics",
            "Fulgent Genetics, Fulgent Genetics": "Fulgent Genetics",
            "Mayo Clinic Laboratories, Mayo Clinic": "Mayo Clinic Labs",
            "Natera, Inc.": "Natera",
            "Genome-Nilou Lab": "Genome-Nilou Lab",
            "Myriad Genetics, Inc.": "Myriad Genetics",
            "Eurofins Ntd Llc (ga)": "Eurofins NTD",
            "All of Us Research Program, National Institutes of Health": "All of Us (NIH)",
            "Quest Diagnostics Nichols Institute San Juan Capistrano": "Quest Diagnostics",
        }
        # Check de-identification with the ORIGINAL name first (highest priority)
        if DEIDENTIFY and name in DEIDENTIFY_MAP:
            return DEIDENTIFY_MAP[name]

        if name in replacements:
            result = replacements[name]
        else:
            result = name
            # Generic shortening
            for suffix in [", Inc", ", Inc.", " LLC", " GmbH"]:
                result = result.replace(suffix, "")

        # Check de-identification with the SHORT name
        if DEIDENTIFY and result in DEIDENTIFY_MAP:
            return DEIDENTIFY_MAP[result]

        if len(result) > max_len:
            result = result[:max_len - 1] + "…"

        return result

    # ── (a) Top 20 horizontal bar ──────────────────────────────────
    top20 = overall[:20]
    names = [short_name(e["submitter"]) for e in top20][::-1]
    counts = [e["total_submissions"] for e in top20][::-1]
    pcts = [e["pct_of_total"] for e in top20][::-1]

    # Color by institution type
    type_colors = {
        "lab": NAT["blue"],
        "consortium": NAT["aqua"],
        "resource": NAT["violet"],
        "clinic": NAT["lime"],
        "LSDB": NAT["yellow"],
        "other": NAT["grey"],
        "unknown": NAT["lt_grey"],
        "patient registry": NAT["pink"],
    }
    bar_colors = [type_colors.get(
        overall[19 - i]["institution_type"], NAT["lt_grey"]
    ) for i in range(20)]

    y_pos = np.arange(20)
    ax_a.barh(y_pos, counts, height=0.72, color=bar_colors,
              edgecolor="white", linewidth=0.3)
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(names, fontsize=5.5)
    ax_a.set_xlabel("Number of submissions (SCV records)")

    # Percentage annotations
    for i, (c, p) in enumerate(zip(counts, pcts)):
        ax_a.text(c + max(counts) * 0.01, i, f"{c:,}  ({p}%)",
                  va="center", fontsize=4.5, color=NAT["grey"])

    ax_a.set_xlim(0, max(counts) * 1.28)
    ax_a.set_title(f"Top 20 submitters ({stats['total_unique_submitters']:,} total organizations)")

    # Type legend
    seen = set()
    type_patches = []
    for e in top20:
        t = e["institution_type"]
        if t not in seen:
            seen.add(t)
            type_patches.append(mpatches.Patch(
                color=type_colors.get(t, NAT["lt_grey"]),
                label=t.capitalize() if t != "LSDB" else "LSDB"
            ))
    ax_a.legend(handles=type_patches, loc="lower right",
                fontsize=5, frameon=True, facecolor="white", edgecolor="none",
                framealpha=0.85, title="Institution type", title_fontsize=5.5)
    panel_label(ax_a, "a")

    # ── (b) Temporal stacked area ──────────────────────────────────
    years_str = temporal["years"]
    # Focus on 2013-2025 (ClinVar era)
    year_mask = [y for y in years_str if 2013 <= int(y) <= 2025]
    year_ints = [int(y) for y in year_mask]

    # Stack top 10 submitters + Other
    stack_names_raw = list(temporal["submitters"].keys())
    # Remove "_Other" and get top 10 by total in the focused period
    real_names = [n for n in stack_names_raw if n != "_Other"]

    def period_total(name):
        return sum(temporal["submitters"][name].get(y, 0) for y in year_mask)

    real_names.sort(key=period_total, reverse=True)
    stack_top10 = real_names[:10]

    # Build stacked area data
    stack_data = []
    stack_labels = []
    stack_colors_list = [NAT["dark_red"], NAT["blue"], NAT["salmon"],
                         NAT["aqua"], NAT["lime"], NAT["violet"],
                         NAT["yellow"], NAT["cyan"], NAT["pink"],
                         NAT["purple"]]
    for i, name in enumerate(stack_top10):
        vals = [temporal["submitters"][name].get(y, 0) for y in year_mask]
        stack_data.append(vals)
        stack_labels.append(short_name(name))

    # Add "Other" (all remaining submitters)
    other_vals = []
    for y in year_mask:
        total_y = stats["submissions_per_year"].get(y, 0)
        top10_y = sum(temporal["submitters"][n].get(y, 0) for n in stack_top10)
        other_y = total_y - top10_y
        other_vals.append(max(0, other_y))
    stack_data.append(other_vals)
    stack_labels.append("Other")
    stack_colors_list.append(NAT["lt_grey"])

    ax_b.stackplot(year_ints, *stack_data,
                   labels=stack_labels, colors=stack_colors_list,
                   alpha=0.85, linewidth=0.3, edgecolor="white")
    ax_b.set_xlim(2013, 2025)
    ax_b.set_xlabel("Year (DateLastEvaluated)")
    ax_b.set_ylabel("Number of submissions")
    ax_b.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{x / 1e6:.1f}M" if x >= 1e6 else (f"{x / 1e3:.0f}K" if x > 0 else "0")))
    ax_b.set_title("Submission volume by top submitters over time")

    # Legend inside the plot area (upper-left, 2 columns to avoid overlap)
    handles_b, labels_b = ax_b.get_legend_handles_labels()
    ax_b.legend(handles_b[::-1], labels_b[::-1],
                loc="upper left", bbox_to_anchor=(0.0, 1.0),
                ncol=2, fontsize=3.8, frameon=True, facecolor="white",
                edgecolor="none", framealpha=0.9, columnspacing=0.8,
                handlelength=1.0, handleheight=0.7, labelspacing=0.3)
    panel_label(ax_b, "b")

    # ── (c) Collection method donut ────────────────────────────────
    method_data = stats.get("collection_method_totals", {})
    # Group small categories
    method_groups = {
        "Clinical testing": method_data.get("clinical testing", 0),
        "Research": method_data.get("research", 0),
        "Literature only": method_data.get("literature only", 0),
        "Curation": method_data.get("curation", 0),
        "Not provided": method_data.get("not provided", 0),
    }
    method_other = sum(v for k, v in method_data.items()
                       if k not in ["clinical testing", "research",
                                    "literature only", "curation", "not provided"])
    method_groups["Other methods"] = method_other
    total_methods = sum(method_groups.values())

    m_labels = list(method_groups.keys())
    m_vals = list(method_groups.values())
    m_colors = [NAT["blue"], NAT["lime"], NAT["yellow"],
                NAT["violet"], NAT["lt_grey"], NAT["grey"]]

    wedges, _ = ax_c.pie(m_vals, colors=m_colors,
                         startangle=90, counterclock=False,
                         wedgeprops={"linewidth": 0.5, "edgecolor": "white",
                                     "width": 0.45})
    ax_c.set_title("Collection methods")
    # Center text
    ax_c.text(0, 0, f"{total_methods:,}\ntotal", ha="center", va="center",
              fontsize=6, fontweight="bold", color=NAT["grey"])
    # Legend
    m_legend = [f"{l} ({v / total_methods * 100:.1f}%)"
                for l, v in zip(m_labels, m_vals)]
    ax_c.legend(wedges, m_legend, loc="upper center",
                bbox_to_anchor=(0.5, -0.05), fontsize=4.5,
                ncol=2, frameon=False)
    panel_label(ax_c, "c", x=-0.22, y=1.12)

    # ── (d) Classification breakdown for top 10 ───────────────────
    sig_order = ["Pathogenic", "Likely pathogenic", "Uncertain significance",
                 "Likely benign", "Benign"]
    sig_colors = [NAT["dark_red"], NAT["salmon"], NAT["yellow"],
                  NAT["blue"], NAT["cyan"]]

    top10_data = overall[:10]
    d_names = [short_name(e["submitter"]) for e in top10_data]
    x_d = np.arange(len(d_names))
    bar_w = 0.65

    bottoms = np.zeros(len(d_names))
    for sig, col in zip(sig_order, sig_colors):
        vals = []
        for e in top10_data:
            bd = e.get("classification_breakdown", {})
            total = e["total_submissions"]
            cnt = bd.get(sig, 0)
            vals.append(100 * cnt / total if total > 0 else 0)
        ax_d.bar(x_d, vals, bar_w, bottom=bottoms, color=col,
                 edgecolor="white", linewidth=0.3, label=sig)
        bottoms += np.array(vals)

    ax_d.set_xticks(x_d)
    ax_d.set_xticklabels(d_names, rotation=45, ha="right", fontsize=4.5)
    ax_d.set_ylabel("Percentage")
    ax_d.set_ylim(0, 105)
    ax_d.set_title("Classification profile by submitter")
    ax_d.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35),
                ncol=3, fontsize=4, frameon=False, columnspacing=1.0)
    panel_label(ax_d, "d", x=-0.18, y=1.12)

    plt.subplots_adjust(bottom=0.14, hspace=0.45)
    _save(fig, out, "figure14_submitter_landscape")

    # ── Check for server-side enhanced data and generate figure15 ──
    # Server output may live in a submitter_analysis/ subdirectory
    server_dir = Path(comment_data_dir) / "submitter_analysis"
    if not server_dir.exists():
        server_dir = Path(comment_data_dir)
    gene_matrix_path = server_dir / "submitter_gene_matrix.json"
    reclass_link_path = server_dir / "submitter_reclassification.json"
    conflict_path = server_dir / "submitter_conflict_pairs.json"
    entry_exit_path = server_dir / "submitter_entry_exit.json"

    has_server_data = all(p.exists() for p in [
        gene_matrix_path, reclass_link_path, conflict_path, entry_exit_path
    ])

    if has_server_data:
        _figure15_submitter_deep_dive(
            gene_matrix_path, reclass_link_path, conflict_path,
            entry_exit_path, temporal, stats, short_name, out
        )
    else:
        missing = [p.name for p in [gene_matrix_path, reclass_link_path,
                                     conflict_path, entry_exit_path]
                   if not p.exists()]
        print(f"  [skip] figure15 — server-side data not found: {', '.join(missing)}")
        print(f"         Run clinvar_submitter_server.py on the server to generate these files.")


def _figure15_submitter_deep_dive(gene_matrix_path, reclass_link_path,
                                   conflict_path, entry_exit_path,
                                   temporal, stats, short_name, out):
    """Four-panel deep dive into submitter dynamics (requires server data).

    (a) Submitter × gene heatmap (top 15 submitters × top 20 genes)
    (b) Reclassification triggers: new submitter vs reinterpretation
    (c) Top conflicting submitter pairs
    (d) Submitter growth: active submitters per snapshot over time
    """
    with open(gene_matrix_path) as f:
        gene_matrix = json.load(f)
    with open(reclass_link_path) as f:
        reclass_data = json.load(f)
    with open(conflict_path) as f:
        conflict_data = json.load(f)
    with open(entry_exit_path) as f:
        entry_exit = json.load(f)

    fig = plt.figure(figsize=(7.08, 9.2))
    gs = GridSpec(2, 2, figure=fig, hspace=0.65, wspace=0.50,
                  top=0.96, bottom=0.06)
    ax_a = fig.add_subplot(gs[0, :])   # heatmap — full width
    ax_b = fig.add_subplot(gs[1, 0])   # reclassification triggers
    ax_c = fig.add_subplot(gs[1, 1])   # conflict pairs

    # ── (a) Submitter × gene heatmap ──────────────────────────────
    submitters = gene_matrix["submitters"][:15]
    genes = gene_matrix["genes"][:20]
    matrix = gene_matrix["matrix"]

    heat = np.zeros((len(submitters), len(genes)))
    for i, sub in enumerate(submitters):
        row_total = sum(matrix.get(sub, {}).values())
        for j, gene in enumerate(genes):
            val = matrix.get(sub, {}).get(gene, 0)
            heat[i, j] = 100 * val / row_total if row_total > 0 else 0

    im = ax_a.imshow(heat, aspect="auto", cmap="YlOrRd",
                     interpolation="nearest")
    ax_a.set_xticks(np.arange(len(genes)))
    ax_a.set_xticklabels(genes, rotation=60, ha="right", fontsize=4.5)
    ax_a.set_yticks(np.arange(len(submitters)))
    ax_a.set_yticklabels([short_name(s) for s in submitters], fontsize=5)
    ax_a.set_title("Submitter specialization: % of submissions per gene")

    cbar = fig.colorbar(im, ax=ax_a, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("% of submitter's total", fontsize=5)
    cbar.ax.tick_params(labelsize=4.5)
    panel_label(ax_a, "a")

    # ── (b) Reclassification triggers ─────────────────────────────
    top_new = reclass_data.get("top_30_new_submitter_triggers", [])[:15]
    top_reint = reclass_data.get("top_30_reinterpretation", [])[:15]

    # Horizontal grouped bar: for each submitter show new-sub vs reinterp
    # Combine both lists into a unified set
    all_trigger_subs = {}
    for item in top_new:
        name = item["submitter"]
        all_trigger_subs[name] = {"new": item["events_triggered"], "reinterp": 0}
    for item in top_reint:
        name = item["submitter"]
        if name in all_trigger_subs:
            all_trigger_subs[name]["reinterp"] = item["events_involved"]
        else:
            all_trigger_subs[name] = {"new": 0, "reinterp": item["events_involved"]}

    # Sort by total and take top 12
    sorted_subs = sorted(all_trigger_subs.items(),
                         key=lambda x: x[1]["new"] + x[1]["reinterp"],
                         reverse=True)[:12]

    b_names = [short_name(s[0], 22) for s in sorted_subs][::-1]
    b_new = [s[1]["new"] for s in sorted_subs][::-1]
    b_reint = [s[1]["reinterp"] for s in sorted_subs][::-1]

    y_b = np.arange(len(b_names))
    bh = 0.35
    ax_b.barh(y_b + bh / 2, b_new, bh, color=NAT["salmon"],
              edgecolor="white", linewidth=0.3, label="New submitter")
    ax_b.barh(y_b - bh / 2, b_reint, bh, color=NAT["blue"],
              edgecolor="white", linewidth=0.3, label="Reinterpretation")
    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(b_names, fontsize=4.5)
    ax_b.set_xlabel("Reclassification events", fontsize=6)
    ax_b.set_title("Reclassification triggers by submitter", fontsize=8)
    ax_b.legend(loc="lower right", fontsize=5,
                frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
    panel_label(ax_b, "b")

    # ── (c) Conflict pairs ────────────────────────────────────────
    pairs = conflict_data.get("top_50_conflict_pairs", [])[:12]
    c_labels = [f"{short_name(p['submitter_a'], 16)} vs\n{short_name(p['submitter_b'], 16)}"
                for p in pairs][::-1]
    c_vals = [p["variant_count"] for p in pairs][::-1]

    y_c = np.arange(len(c_labels))
    ax_c.barh(y_c, c_vals, height=0.65, color=NAT["violet"],
              edgecolor="white", linewidth=0.3)
    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels(c_labels, fontsize=4.2, linespacing=0.85)
    ax_c.set_xlabel("Variants with conflicting classifications", fontsize=6)
    ax_c.set_title("Most frequent inter-submitter conflicts", fontsize=8)
    ax_c.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v/1000:.0f}K" if v >= 1000 else f"{v:.0f}"))
    ax_c.tick_params(axis='x', labelsize=5)

    total_conflicts = conflict_data.get("total_conflict_variant_instances", 0)
    n_pairs = conflict_data.get("unique_conflict_pairs", 0)
    ax_c.text(0.95, 0.05,
              f"{n_pairs:,} unique pairs\n{total_conflicts:,} total instances",
              transform=ax_c.transAxes, fontsize=4.5, ha="right",
              va="bottom", color=NAT["grey"])
    panel_label(ax_c, "c")

    _save(fig, out, "figure15_submitter_deep_dive")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 16 — Case Study: Landmark Variant Reclassification Trajectories
# ═══════════════════════════════════════════════════════════════════════

def figure16_case_study_trajectories(comment_data_dir, out):
    """Classification trajectories of three landmark reclassified variants.

    Panels:
      (a) TNNI3 p.Pro82Ser   — Manrai et al. NEJM 2016 (HCM/diversity)
      (b) MYBPC3 p.Gly278Glu — Manrai et al. NEJM 2016 (HCM/diversity)
      (c) UTP4/CIRH1A p.Arg565Trp — Lek et al. Nature 2016 (ExAC)
    """
    import re as _re

    print("\n── Figure 16: Case Study — Variant Reclassification Trajectories ──")

    data_dir = Path(comment_data_dir)
    reversal_path = data_dir / "plp_blb_reversal_variants.csv"
    events_path  = data_dir / "reclassification_events_detail.csv"

    for p, label in [(reversal_path, "plp_blb_reversal_variants.csv"),
                     (events_path,  "reclassification_events_detail.csv")]:
        if not p.exists():
            print(f"  [skip] {label} not found — skipping figure 16")
            return

    # ── Helpers ────────────────────────────────────────────────────────
    CAT_Y = OrderedDict([
        ("P", 8), ("LP", 7), ("P/LP", 6), ("VUS", 5),
        ("Conflicting", 4), ("Other", 3),
        ("LB", 2), ("B/LB", 1), ("B", 0),
    ])
    CAT_DISPLAY = {
        "P": "Pathogenic", "LP": "Likely path.", "P/LP": "Path./L.path.",
        "VUS": "VUS", "Conflicting": "Conflicting", "Other": "Other",
        "LB": "Likely benign", "B/LB": "Ben./L.ben.", "B": "Benign",
    }

    def _df(d):
        """'YYYY-MM' → float year."""
        y, m = d.split("-")
        return int(y) + (int(m) - 1) / 12

    _EV_RE = _re.compile(
        r'(\d{4}-\d{2}):([^→]+)→(\d{4}-\d{2}):([^\[]+)\[([^\]]+)\]'
    )

    def _parse_eh(s):
        """Parse Event_History field into list of transition dicts."""
        return [{"fd": m.group(1), "fc": m.group(2).strip(),
                 "td": m.group(3), "tc": m.group(4).strip(),
                 "drv": m.group(5)}
                for m in _EV_RE.finditer(s)]

    # ── Read reversal CSV for TNNI3 & MYBPC3 ─────────────────────────
    var_data = {}
    want = {"VCV000012421", "VCV000042797"}
    with open(reversal_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            vid = row["VCV_ID"]
            if vid not in want:
                continue
            evs = _parse_eh(row["Event_History"])
            init_cat = row["Change_Type"].split("->")[0]
            traj = [(_df(row["First_Date"]), init_cat)]
            for e in evs:
                traj.append((_df(e["td"]), e["tc"]))
            traj.append((2026.25, traj[-1][1]))
            var_data[vid] = {
                "traj": traj,
                "n_events": int(row["Num_Reclass_Events"]),
                "n_subs": int(row["Num_Submissions"]),
                "current_cat": row["Current_Category"],
            }
            want.discard(vid)
            if not want:
                break

    # ── Read events CSV for UTP4/CIRH1A variant 3194 ──────────────────
    utp4_rows = []
    with open(events_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["Variation_ID"] == "3194" and row["Gene"] == "UTP4":
                utp4_rows.append(row)

    if utp4_rows:
        traj = [(2012 + 10 / 12, utp4_rows[0]["Prev_Category"])]
        for ev in utp4_rows:
            traj.append((_df(ev["Curr_Date"]), ev["Curr_Category"]))
        traj.append((2026.25, traj[-1][1]))
        var_data["UTP4_3194"] = {
            "traj": traj,
            "n_events": len(utp4_rows),
            "n_subs": int(utp4_rows[-1]["Curr_Num_Subs"]),
            "current_cat": traj[-2][1],
        }

    if len(var_data) < 3:
        print(f"  [skip] Found only {len(var_data)}/3 target variants")
        return

    # ── Panel specifications ──────────────────────────────────────────
    panel_specs = [
        {"key": "VCV000012421", "gene": "TNNI3", "var": "p.Pro82Ser",
         "paper_label": "Manrai et al. NEJM 2016",
         "pub": 2016 + 7 / 12},
        {"key": "VCV000042797", "gene": "MYBPC3", "var": "p.Gly278Glu",
         "paper_label": "Manrai et al. NEJM 2016",
         "pub": 2016 + 7 / 12},
        {"key": "UTP4_3194", "gene": "UTP4", "var": "p.Arg565Trp",
         "alias": "CIRH1A",
         "paper_label": "Lek et al. Nature 2016 (ExAC)",
         "pub": 2016 + 7 / 12},
    ]

    # ── Key annotations per panel (date_str, text, y_offset) ─────────
    KEY_ANNOTS = {
        "VCV000012421": [
            ("2014-09", "First\ndowngrade", -1.6),
            ("2019-04", "Settled\nB/LB", -1.4),
        ],
        "VCV000042797": [
            ("2013-12", "Batch\ndowngrade", -1.6),
            ("2017-11", "Settled\nBenign", -1.4),
        ],
        "UTP4_3194": [
            ("2016-10", "Population\nfrequency", -1.6),
            ("2020-11", "Second\nsubmitter", -1.4),
            ("2022-04", "New\ndispute", 1.4),
        ],
    }

    # ── Create figure ─────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(7.08, 7.5), sharex=True)
    plt.subplots_adjust(hspace=0.45, top=0.94, bottom=0.07,
                        left=0.14, right=0.90)

    for idx, (ax, spec) in enumerate(zip(axes, panel_specs)):
        d = var_data[spec["key"]]
        traj = d["traj"]

        # ---- background zones ----
        ax.axhspan(5.5, 8.5, color="#FFEBEE", alpha=0.45, zorder=0, lw=0)
        ax.axhspan(2.5, 5.5, color="#FFF8E1", alpha=0.35, zorder=0, lw=0)
        ax.axhspan(-0.5, 2.5, color="#E3F2FD", alpha=0.45, zorder=0, lw=0)

        # Zone labels on far right
        for ztxt, zy, zcol in [("Pathogenic", 7.0, NAT["dark_red"]),
                                ("Uncertain",  4.5, NAT["violet"]),
                                ("Benign",     1.0, NAT["blue"])]:
            ax.text(2026.0, zy, ztxt, fontsize=3.5, color=zcol,
                    alpha=0.5, ha="left", va="center", style="italic")

        # ---- trajectory: coloured horizontal segments ----
        for j in range(len(traj) - 1):
            t0, c0 = traj[j]
            t1, _  = traj[j + 1]
            y0 = CAT_Y.get(c0, 3)
            ax.plot([t0, t1], [y0, y0],
                    color=CATEGORY_COLORS.get(c0, NAT["grey"]),
                    linewidth=2.5, solid_capstyle="butt", zorder=3)
            # vertical connector
            if j < len(traj) - 2:
                _, c1 = traj[j + 1]
                y1 = CAT_Y.get(c1, 3)
                if y0 != y1:
                    ax.plot([t1, t1], [y0, y1], color="#BBBBBB",
                            linewidth=0.4, zorder=2)

        # ---- transition markers ----
        for j in range(1, len(traj) - 1):
            t, c = traj[j]
            ax.plot(t, CAT_Y.get(c, 3), "o",
                    color=CATEGORY_COLORS.get(c, NAT["grey"]),
                    markersize=3, zorder=4,
                    markeredgecolor="white", markeredgewidth=0.3)

        # ---- paper publication line ----
        ax.axvline(spec["pub"], color="#333333", lw=0.8, ls="--",
                   zorder=2, alpha=0.55)
        ax.annotate(spec["paper_label"],
                    xy=(spec["pub"], 8.35),
                    xytext=(spec["pub"] + 0.25, 8.35),
                    fontsize=4.5, va="top", ha="left", style="italic",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              ec="#AAAAAA", lw=0.3, alpha=0.92),
                    arrowprops=dict(arrowstyle="-", color="#AAAAAA",
                                   lw=0.3))

        # ---- key transition annotations ----
        for date_str, text, y_off in KEY_ANNOTS.get(spec["key"], []):
            tx = _df(date_str)
            # find the category at this date in the trajectory
            cat_at = "VUS"
            for jj in range(len(traj) - 1):
                if traj[jj][0] <= tx <= traj[jj + 1][0]:
                    cat_at = traj[jj + 1][1] if tx >= traj[jj + 1][0] - 0.01 else traj[jj][1]
                    break
            # use the destination category (transition dot is at the new category)
            for jj in range(len(traj)):
                if abs(traj[jj][0] - tx) < 0.02:
                    cat_at = traj[jj][1]
                    break
            base_y = CAT_Y.get(cat_at, 3)
            ax.annotate(text, xy=(tx, base_y),
                        xytext=(tx, base_y + y_off),
                        fontsize=4, ha="center", color="#444444",
                        arrowprops=dict(arrowstyle="->",
                                        color="#999999", lw=0.4),
                        bbox=dict(boxstyle="round,pad=0.1", fc="white",
                                  ec="none", alpha=0.8))

        # ---- current status ----
        cur = d["current_cat"]
        cur_y = CAT_Y.get(cur, 3)
        ax.text(2026.0, cur_y, f" {CAT_DISPLAY.get(cur, cur)}",
                fontsize=4.5, va="center", ha="left",
                color=CATEGORY_COLORS.get(cur, NAT["grey"]),
                fontweight="bold", alpha=0.8,
                bbox=dict(boxstyle="round,pad=0.08", fc="white",
                          ec="none", alpha=0.7))

        # ---- title ----
        alias = f" ({spec['alias']})" if "alias" in spec else ""
        ax.set_title(
            f"{spec['gene']}{alias} {spec['var']}  \u2014  "
            f"{d['n_events']} reclassification events, "
            f"{d['n_subs']} submitters",
            fontsize=7, fontweight="bold", loc="left", pad=4)

        # ---- y-axis ----
        ax.set_yticks(list(CAT_Y.values()))
        ax.set_yticklabels([CAT_DISPLAY[c] for c in CAT_Y], fontsize=5.5)
        ax.set_ylim(-0.8, 8.8)
        ax.yaxis.grid(True, lw=0.15, alpha=0.3, color="#CCCCCC")
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        panel_label(ax, chr(ord("a") + idx))

    # ---- shared x-axis ----
    axes[-1].set_xlim(2012.5, 2026.5)
    axes[-1].set_xlabel("Year", fontsize=7)
    xt = list(range(2013, 2027))
    axes[-1].set_xticks(xt)
    axes[-1].set_xticklabels([str(y) for y in xt], fontsize=5.5,
                              rotation=45, ha="right")

    _save(fig, Path(out), "figure16_case_study_trajectories")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 17 — Reclassification Driver Classification Framework
# ═══════════════════════════════════════════════════════════════════════

def figure17_conceptual_framework(comment_data_dir, out):
    """Decision-tree flowchart showing how reclassification events are
    classified into driver categories, with downstream outcome summary."""
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
    from matplotlib.collections import PatchCollection

    print("Figure 17: Reclassification driver framework ...")

    # ── load data ───────────────────────────────────────────────────
    driver_path = Path(comment_data_dir) / "reclassification_drivers.json"
    ttr_path = Path(comment_data_dir) / "time_to_reclassification.json"
    rev_path = Path(comment_data_dir) / "plp_blb_reversal_summary.json"

    drivers, ttr, rev = {}, {}, {}
    for p, d in [(driver_path, "drivers"), (ttr_path, "ttr"),
                 (rev_path, "rev")]:
        if p.exists():
            with open(p) as f:
                locals()[d] = json.load(f)

    if driver_path.exists():
        with open(driver_path) as f:
            drivers = json.load(f)
    if ttr_path.exists():
        with open(ttr_path) as f:
            ttr = json.load(f)
    if rev_path.exists():
        with open(rev_path) as f:
            rev = json.load(f)

    ds = drivers.get("driver_summary", {})
    total = drivers.get("total_reclassification_events", 697698)

    # ── figure setup ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9.2, 14))
    ax.set_xlim(0, 9.2)
    ax.set_ylim(0, 14)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.01)

    # ── colours ─────────────────────────────────────────────────────
    C_INPUT  = "#2c3e50"    # dark header
    C_DIAM   = "#f0f4f8"    # diamond fill
    C_DIAM_E = "#5b7a9d"    # diamond edge
    C_BATCH  = "#7f8c8d"    # grey
    C_SCHEMA = "#c0764a"    # brown/orange
    C_AMBIG  = "#95a5a6"    # light grey
    C_KNOWL  = "#5b9ea6"    # teal
    C_CONFL  = "#c0392b"    # red
    C_RESOL  = "#2980b9"    # blue
    C_NEWSUB = "#27ae60"    # green
    C_REINT  = "#8e44ad"    # purple
    C_Q_BG   = "#f0f4f8"    # question box bg
    C_SCORE  = "#f8f9fa"    # scoring box bg
    C_YES    = "#2d8a4e"
    C_GREY   = "#666666"
    C_OUT    = NAT["aqua"]

    def _rbox(x, y, w, h, fc, ec=None, lw=1.2, zorder=3):
        """Rounded rectangle patch."""
        ec = ec or fc
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                           fc=fc, ec=ec, lw=lw, zorder=zorder)
        ax.add_patch(b)
        return b

    def _txt(x, y, text, **kw):
        kw.setdefault("ha", "center")
        kw.setdefault("va", "center")
        kw.setdefault("fontsize", 8)
        kw.setdefault("zorder", 5)
        return ax.text(x, y, text, **kw)

    def _diamond(cx, cy, hw, hh):
        """Draw a diamond centred at (cx,cy) with half-width hw, half-height hh."""
        verts = [(cx, cy + hh), (cx + hw, cy), (cx, cy - hh), (cx - hw, cy)]
        p = Polygon(verts, closed=True, fc=C_DIAM, ec=C_DIAM_E, lw=1.5, zorder=3)
        ax.add_patch(p)

    def _varrow(x, y1, y2, color=C_GREY, lw=1.3):
        """Vertical arrow from (x,y1) to (x,y2)."""
        ax.annotate("", xy=(x, y2), xytext=(x, y1),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=lw, mutation_scale=10),
                    zorder=2)

    def _harrow(x1, x2, y, color=C_GREY, lw=1.3):
        """Horizontal arrow from (x1,y) to (x2,y)."""
        ax.annotate("", xy=(x2, y), xytext=(x1, y),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=lw, mutation_scale=10),
                    zorder=2)

    CX = 4.6   # centre x for the main flow

    # ================================================================
    # INPUT BOX
    # ================================================================
    _rbox(2.0, 13.2, 5.2, 0.55, C_INPUT)
    _txt(CX, 13.58, "Reclassification Events Detected",
         fontsize=10, fontweight="bold", color="white")
    _txt(CX, 13.36, f"N = {total:,} classification changes across 148 monthly snapshots",
         fontsize=7.5, color="white")

    # ================================================================
    # STEP 1 — BATCH CHECK
    # ================================================================
    _varrow(CX, 13.2, 12.85)
    _txt(CX, 12.75, "STEP 1", fontsize=7, fontweight="bold", color="#888",
         fontstyle="italic")
    _diamond(CX, 12.25, 1.3, 0.4)
    _txt(CX, 12.3, "Did >5% of all tracked variants", fontsize=7)
    _txt(CX, 12.15, "change in this snapshot?", fontsize=7)

    # YES → left: Batch Event
    _harrow(CX - 1.3, 1.7, 12.25, color=C_YES)
    _txt(2.8, 12.37, "YES", fontsize=7, fontweight="bold", color=C_YES)

    _rbox(0.15, 11.7, 1.55, 0.85, C_BATCH)
    _txt(0.92, 12.35, "BATCH EVENT", fontsize=8, fontweight="bold", color="white")
    _txt(0.92, 12.15, f"n = {ds.get('batch_event',{}).get('count',89219):,} ({ds.get('batch_event',{}).get('pct',12.8)}%)",
         fontsize=7, color="white")
    _txt(0.92, 11.92, "5 episodes: Dec 2013,", fontsize=5.5,
         color="white", fontstyle="italic")
    _txt(0.92, 11.8, "Jan\u201314, Feb\u201314, Jun\u201318, Jan\u201320",
         fontsize=5.5, color="white", fontstyle="italic")

    # NO → down
    _varrow(CX, 11.85, 11.45)
    _txt(CX + 0.15, 11.6, "NO", fontsize=7, fontweight="bold", color=C_GREY,
         ha="left")

    # ================================================================
    # STEP 2 — TRANSITION WINDOW
    # ================================================================
    _txt(CX, 11.35, "STEP 2", fontsize=7, fontweight="bold", color="#888",
         fontstyle="italic")
    _diamond(CX, 10.85, 1.15, 0.35)
    _txt(CX, 10.88, "Event within", fontsize=7)
    _txt(CX, 10.75, "Jan\u2013Jun 2020?", fontsize=7)

    # YES → left: Weighted Scoring
    _harrow(CX - 1.15, 2.5, 10.85, color=C_YES)
    _txt(3.3, 10.97, "YES", fontsize=7, fontweight="bold", color=C_YES)

    # Scoring box
    _rbox(0.15, 9.9, 2.35, 0.95, C_SCORE, ec=C_DIAM_E, lw=1.0)
    _txt(1.32, 10.73, "Weighted Signal Scoring", fontsize=7.5,
         fontweight="bold", color="#333")
    # Schema signals
    _txt(1.32, 10.52, "Schema signals (+1 each):", fontsize=6,
         fontweight="bold", color="#333")
    _txt(1.32, 10.38, "S1: Event at primary boundary (Dec\u2192Jan)",
         fontsize=5.5, color="#444")
    _txt(1.32, 10.26, "S2: Submission count changed",
         fontsize=5.5, color="#444")
    _txt(1.32, 10.14, "S3: Review status vocabulary changed",
         fontsize=5.5, color="#444")
    # Evidence signals
    _txt(1.32, 9.99, "Evidence signals: E1\u2013E3 (negation of S1\u2013S3)",
         fontsize=5.5, fontweight="bold", color="#333")

    # Arrow down from scoring box
    _varrow(1.32, 9.9, 9.5)

    # Three transition outcomes — counts INSIDE boxes
    bw = 0.73   # box width for transition outcomes
    bh = 0.7    # taller to fit count line
    gap = 0.06
    x0 = 0.15
    _rbox(x0, 8.7, bw, bh, C_SCHEMA)
    _txt(x0 + bw/2, 9.27, "LIKELY", fontsize=6, fontweight="bold", color="white")
    _txt(x0 + bw/2, 9.14, "SCHEMA-", fontsize=6, fontweight="bold", color="white")
    _txt(x0 + bw/2, 9.01, "DRIVEN", fontsize=6, fontweight="bold", color="white")
    _txt(x0 + bw/2, 8.84, f"n={ds.get('likely_schema_driven',{}).get('count',17907):,}",
         fontsize=5.5, color="white")
    _txt(x0 + bw/2, 8.74, "(2.6%)", fontsize=5.5, color="white")

    x1 = x0 + bw + gap
    _rbox(x1, 8.7, bw, bh, C_AMBIG)
    _txt(x1 + bw/2, 9.27, "TRANSITION", fontsize=5.5, fontweight="bold", color="white")
    _txt(x1 + bw/2, 9.14, "AMBIGUOUS", fontsize=5.5, fontweight="bold", color="white")
    _txt(x1 + bw/2, 8.84, "n=105", fontsize=5.5, color="white")
    _txt(x1 + bw/2, 8.74, "(<0.1%)", fontsize=5.5, color="white")

    x2 = x1 + bw + gap
    _rbox(x2, 8.7, bw, bh, C_KNOWL)
    _txt(x2 + bw/2, 9.27, "LIKELY", fontsize=6, fontweight="bold", color="white")
    _txt(x2 + bw/2, 9.14, "EVIDENCE-", fontsize=6, fontweight="bold", color="white")
    _txt(x2 + bw/2, 9.01, "BASED", fontsize=6, fontweight="bold", color="white")
    _txt(x2 + bw/2, 8.84, f"n={ds.get('likely_knowledge_driven',{}).get('count',21525):,}",
         fontsize=5.5, color="white")
    _txt(x2 + bw/2, 8.74, "(3.1%)", fontsize=5.5, color="white")

    # NO → down from Step 2
    _varrow(CX, 10.50, 10.05)
    _txt(CX + 0.15, 10.25, "NO", fontsize=7, fontweight="bold", color=C_GREY,
         ha="left")

    # ================================================================
    # STEP 3 — OBSERVABLE MECHANISM CLASSIFICATION
    # ================================================================
    _txt(CX, 9.95, "STEP 3: OBSERVABLE MECHANISM CLASSIFICATION",
         fontsize=7, fontweight="bold", color="#888", fontstyle="italic")

    # ── Q1: Conflict Emerged ────────────────────────────────────────
    qw, qh = 2.0, 0.38
    qx = CX - qw / 2
    q1y = 9.3
    _rbox(qx, q1y, qw, qh, C_Q_BG, ec=C_DIAM_E, lw=1.0)
    _txt(CX, q1y + 0.25, 'Aggregate classification changed', fontsize=7, color="#333")
    _txt(CX, q1y + 0.1, 'TO "conflicting interpretations"?', fontsize=7, color="#333")

    _varrow(CX, 9.75, q1y + qh)

    # YES → right
    _harrow(CX + qw/2, 7.1, q1y + qh/2, color=C_YES)
    _txt(CX + qw/2 + 0.35, q1y + qh/2 + 0.1, "YES", fontsize=7,
         fontweight="bold", color=C_YES)

    _rbox(7.1, q1y - 0.02, 1.95, 0.42, C_CONFL)
    _txt(8.08, q1y + 0.27, "CONFLICT EMERGED", fontsize=8,
         fontweight="bold", color="white")
    _txt(8.08, q1y + 0.1, f"n = {ds.get('conflict_emerged',{}).get('count',161722):,} (23.2%)",
         fontsize=7, color="white")

    # NO → down
    _varrow(CX, q1y, 8.52)
    _txt(CX + 0.15, q1y - 0.2, "NO", fontsize=7, fontweight="bold",
         color=C_GREY, ha="left")

    # ── Q2: Conflict Resolved ──────────────────────────────────────
    q2y = 8.05
    _rbox(qx, q2y, qw, qh, C_Q_BG, ec=C_DIAM_E, lw=1.0)
    _txt(CX, q2y + 0.25, 'Aggregate classification changed', fontsize=7, color="#333")
    _txt(CX, q2y + 0.1, 'FROM "conflicting interpretations"?', fontsize=7, color="#333")

    _harrow(CX + qw/2, 7.1, q2y + qh/2, color=C_YES)
    _txt(CX + qw/2 + 0.35, q2y + qh/2 + 0.1, "YES", fontsize=7,
         fontweight="bold", color=C_YES)

    _rbox(7.1, q2y - 0.02, 1.95, 0.42, C_RESOL)
    _txt(8.08, q2y + 0.27, "CONFLICT RESOLVED", fontsize=8,
         fontweight="bold", color="white")
    _txt(8.08, q2y + 0.1, f"n = {ds.get('conflict_resolved',{}).get('count',47909):,} (6.9%)",
         fontsize=7, color="white")

    # NO → down
    _varrow(CX, q2y, 7.27)
    _txt(CX + 0.15, q2y - 0.2, "NO", fontsize=7, fontweight="bold",
         color=C_GREY, ha="left")

    # ── Q3: New Submitter ──────────────────────────────────────────
    q3y = 6.8
    _rbox(qx, q3y, qw, 0.32, C_Q_BG, ec=C_DIAM_E, lw=1.0)
    _txt(CX, q3y + 0.16, 'Submission count increased?', fontsize=7, color="#333")

    _harrow(CX + qw/2, 7.1, q3y + 0.16, color=C_YES)
    _txt(CX + qw/2 + 0.35, q3y + 0.26, "YES", fontsize=7,
         fontweight="bold", color=C_YES)

    _rbox(7.1, q3y - 0.05, 1.95, 0.42, C_NEWSUB)
    _txt(8.08, q3y + 0.24, "NEW SUBMITTER", fontsize=8,
         fontweight="bold", color="white")
    _txt(8.08, q3y + 0.07, f"n = {ds.get('new_submitter',{}).get('count',127652):,} (18.3%)",
         fontsize=7, color="white")

    # NO → down
    _varrow(CX, q3y, 6.3)
    _txt(CX + 0.15, q3y - 0.2, "NO", fontsize=7, fontweight="bold",
         color=C_GREY, ha="left")

    # ── Default: Reinterpretation ──────────────────────────────────
    re_y = 5.8
    _rbox(qx, re_y, qw, 0.42, C_REINT)
    _txt(CX, re_y + 0.28, "REINTERPRETATION", fontsize=9,
         fontweight="bold", color="white")
    _txt(CX, re_y + 0.1, f"n = {ds.get('reinterpretation',{}).get('count',231659):,} (33.2%)",
         fontsize=7.5, color="white")

    # ================================================================
    # SUMMARY STACKED BAR
    # ================================================================
    _txt(CX, 5.45, "Driver Distribution Summary", fontsize=9,
         fontweight="bold", color="#333")

    bar_y = 5.0
    bar_h = 0.25
    bar_left = 0.6
    bar_total_w = 8.0

    # Proportions → widths
    cats = [
        (C_BATCH,  12.8, "Batch\n12.8%"),
        (C_SCHEMA,  2.6, ""),
        (C_AMBIG,   0.02, ""),
        (C_KNOWL,   3.1, ""),
        (C_CONFL,  23.2, "Conflict\nEmerged\n23.2%"),
        (C_RESOL,   6.9, "Resolved\n6.9%"),
        (C_NEWSUB, 18.3, "New Sub.\n18.3%"),
        (C_REINT,  33.2, "Reinterpretation\n33.2%"),
    ]
    cx = bar_left
    for color, pct, label in cats:
        w = bar_total_w * pct / 100.0
        ax.add_patch(FancyBboxPatch((cx, bar_y), w, bar_h,
                     boxstyle="square,pad=0", fc=color, ec="none", zorder=3))
        if label:
            _txt(cx + w/2, bar_y - 0.22, label, fontsize=5.5,
                 color=color, fontweight="bold")
        cx += w

    # Transition group label
    trans_x = bar_left + bar_total_w * 12.8 / 100
    trans_w = bar_total_w * 5.72 / 100
    _txt(trans_x + trans_w/2, bar_y - 0.2, "Transition\n5.7%",
         fontsize=5.5, color="#888", fontweight="bold")

    # Bar outline
    ax.add_patch(FancyBboxPatch((bar_left, bar_y), bar_total_w, bar_h,
                 boxstyle="round,pad=0.01", fc="none", ec="#ccc",
                 lw=0.5, zorder=4))

    # ================================================================
    # DOWNSTREAM OUTCOMES
    # ================================================================
    _txt(CX, 4.25, "Downstream Outcomes", fontsize=9,
         fontweight="bold", color=C_OUT)

    n_reclass = ttr.get("overall", {}).get("n", 342057)
    med_yr = ttr.get("overall", {}).get("median_years", 2.58)
    n_rev = rev.get("total_reversal_variants", 6340)

    ow, oh = 2.65, 0.7
    oy = 3.35

    _rbox(0.3, oy, ow, oh, "#E8F8E8", ec=C_OUT, lw=1.2)
    _txt(0.3 + ow/2, oy + oh - 0.17,
         "Variant-Level Impact", fontsize=7.5, fontweight="bold", color=C_OUT)
    _txt(0.3 + ow/2, oy + 0.3,
         f"{n_reclass:,} reclassified (9.3%)", fontsize=6.5, color="#333")
    _txt(0.3 + ow/2, oy + 0.13,
         f"Median {med_yr:.1f} yr  |  63,495 oscillating", fontsize=6, color="#555")

    _rbox(3.3, oy, ow, oh, "#E8F8E8", ec=C_OUT, lw=1.2)
    _txt(3.3 + ow/2, oy + oh - 0.17,
         "System-Level Impact", fontsize=7.5, fontweight="bold", color=C_OUT)
    _txt(3.3 + ow/2, oy + 0.3,
         "Concordance: 72.1% (Jul 2025)", fontsize=6.5, color="#333")
    _txt(3.3 + ow/2, oy + 0.13,
         "Cancer 16.6% vs non-cancer 8.2%", fontsize=6, color="#555")

    _rbox(6.3, oy, ow, oh, "#E8F8E8", ec=C_OUT, lw=1.2)
    _txt(6.3 + ow/2, oy + oh - 0.17,
         "Evidence & Clinical Impact", fontsize=7.5, fontweight="bold", color=C_OUT)
    _txt(6.3 + ow/2, oy + 0.3,
         "Only 4.5% have ACMG codes", fontsize=6.5, color="#333")
    _txt(6.3 + ow/2, oy + 0.13,
         f"{n_rev:,} P\u2194B reversals", fontsize=6, color="#555")

    # Arrows from bar to outcomes
    _varrow(CX, bar_y, 4.35, color="#bbbbbb", lw=1.0)
    _varrow(CX, 4.15, oy + oh, color="#bbbbbb", lw=1.0)

    # ── TITLE ───────────────────────────────────────────────────────
    ax.set_title(
        "Reclassification Driver Classification Framework",
        fontsize=12, fontweight="bold", pad=10)

    _save(fig, Path(out), "figure17_conceptual_framework")


# ═══════════════════════════════════════════════════════════════════════
# SAVE HELPER
# ═══════════════════════════════════════════════════════════════════════

def _save(fig, out_dir, name):
    """Save figure as PDF (vector), high-res PNG, and TIFF for journal submission."""
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  Saved {name}.pdf / .png (600dpi) / .tiff (600dpi)")


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="ClinVar Longitudinal Evolution — Publication-Quality Visualizations"
    )
    parser.add_argument("--data-dir", default="~/Downloads",
                        help="Directory containing ClinVar output files")
    parser.add_argument("--output-dir", "-o", default="~/Downloads/clinvar_figures",
                        help="Directory for output figures")
    parser.add_argument("--skip-history", action="store_true",
                        help="Skip history processing (produces only figs 2, 4, 6)")
    parser.add_argument("--concordance-file", type=str, default=None,
                        help="Path to concordance_timeseries.json")
    parser.add_argument("--domain-file", type=str, default=None,
                        help="Path to domain_summary.json")
    parser.add_argument("--driver-file", type=str, default=None,
                        help="Path to reclassification_drivers.json")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).expanduser()
    out_dir = Path(args.output_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("ClinVar Longitudinal Evolution — Visualisation Pipeline")
    print("=" * 70)
    print(f"Data directory:   {data_dir}")
    print(f"Output directory: {out_dir}")
    print()

    # ── Locate input files ───────────────────────────────────────────
    summary_path = find_file(data_dir, "summary_*.csv")
    history_path = find_file(data_dir, "classification_history_*.csv")
    stats_path = find_file(data_dir, "statistics_*.json")
    xml_path = find_file(data_dir, "xml_structure_analysis.json")
    acmg_stats_path = find_file(data_dir, "acmg_criteria_detailed.json")

    if not summary_path:
        print("ERROR: summary_*.csv not found in data directory")
        sys.exit(1)
    if not stats_path:
        print("ERROR: statistics_*.json not found in data directory")
        sys.exit(1)

    print("Input files:")
    print(f"  Summary:  {summary_path.name}")
    print(f"  History:  {history_path.name if history_path else 'NOT FOUND (skipping figs 1, 3)'}")
    print(f"  Stats:    {stats_path.name}")
    print(f"  XML:      {xml_path.name if xml_path else 'not available (skip fig 5)'}")
    print(f"  ACMG:     {acmg_stats_path.name if acmg_stats_path else 'not available (skip fig 9)'}")
    print()

    # ── Load data ────────────────────────────────────────────────────
    set_nature_style()

    print("Loading statistics ...")
    with open(stats_path) as f:
        stats = json.load(f)

    print("Aggregating summary CSV ...")
    summary_agg = aggregate_summary(summary_path, cache_dir / "summary_agg.json")

    hist_agg = None
    if history_path and not args.skip_history:
        print("Aggregating classification history CSV ...")
        hist_agg = aggregate_history(history_path, cache_dir / "history_agg.json")

    xml_data = None
    if xml_path:
        print("Loading XML structure analysis ...")
        with open(xml_path) as f:
            xml_data = json.load(f)

    acmg_stats = None
    if acmg_stats_path:
        print("Loading ACMG criteria statistics ...")
        with open(acmg_stats_path) as f:
            acmg_stats = json.load(f)

    # ── Load optional new data files ─────────────────────────────────
    # Also search in ClinVar_Research_Project/data/ sub-directory
    project_data_dir = data_dir / "ClinVar_Research_Project" / "data"

    def _load_optional_json(cli_arg, *search_patterns):
        """Try CLI arg first, then search for each pattern across ALL dirs
        before falling back to the next pattern. This ensures v2 files are
        preferred over v1 even if v1 exists in an earlier-searched directory."""
        if cli_arg:
            p = Path(cli_arg).expanduser()
            if p.exists():
                with open(p) as f:
                    return json.load(f)
        dirs = [d for d in [data_dir, project_data_dir] if d.exists()]
        for pat in search_patterns:
            for parent in dirs:
                found = find_file(parent, pat)
                if found:
                    with open(found) as f:
                        return json.load(f)
        return None

    concordance_data = _load_optional_json(
        args.concordance_file,
        "concordance_timeseries_v2.json",
        "concordance_timeseries.json",
    )
    if concordance_data:
        print("Loading concordance data ... OK")

    domain_data = _load_optional_json(
        args.domain_file,
        "domain_summary.json",
    )
    if domain_data:
        print("Loading domain analysis data ... OK")

    driver_data = _load_optional_json(
        args.driver_file,
        "reclassification_drivers.json",
    )
    if driver_data:
        print("Loading reclassification driver data ... OK")

    # ── Generate metadata stats ──────────────────────────────────────
    print()
    print("Generating metadata statistics ...")
    generate_metadata_table(summary_agg, hist_agg, stats, xml_data,
                            acmg_stats, out_dir)

    # ── Generate figures ─────────────────────────────────────────────
    print()
    print("Generating figures ...")

    if hist_agg:
        figure1(hist_agg, out_dir)                           # Q3, Q5
        figure3(hist_agg, out_dir, concordance_data)         # Q4 + concordance

    figure2(summary_agg, stats, out_dir)    # Reclassification dynamics
    figure4(summary_agg, stats, out_dir)    # Q6, Q7
    figure6(summary_agg, stats, out_dir)    # Supplementary
    figure7(summary_agg, out_dir)           # Q7 temporal

    if xml_data:
        figure5(xml_data, out_dir, hist_agg)                 # Q1, Q2, Q6 (revised)
        figure8(xml_data, hist_agg, out_dir)                 # Q2, Q3, Q4 (3-panel)

    if acmg_stats:
        figure9(acmg_stats, out_dir)        # ACMG/AMP criteria analysis

    if domain_data:
        figure10_disease_domain(domain_data, out_dir)        # Disease domain analysis
        figure_supp_rare_vs_common(domain_data, out_dir)     # Rare vs common supplement

    if driver_data:
        figure11_reclassification_drivers(driver_data, out_dir)  # Schema vs knowledge

    # Comment-driven supplementary figures (from clinvar_comment_analyses.py)
    comment_data_dir = data_dir / "ClinVar_Research_Project" / "data"
    if not comment_data_dir.exists():
        comment_data_dir = data_dir  # fallback
    figure12_time_to_reclassification(comment_data_dir, out_dir)
    figure13_batch_event_decomposition(comment_data_dir, out_dir)
    figure14_submitter_landscape(comment_data_dir, out_dir)
    figure16_case_study_trajectories(comment_data_dir, out_dir)
    figure17_conceptual_framework(comment_data_dir, out_dir)

    print()
    print("=" * 70)
    print("All figures saved to:", out_dir)
    print("=" * 70)

    if not hist_agg:
        print("\nNote: Figures 1 and 3 were skipped because classification_history")
        print("was not found or --skip-history was used.")
    if not xml_data:
        print("\nNote: Figure 5 (XML structural evolution) was skipped because")
        print("xml_structure_analysis.json was not found.")
        print("To generate it, run clinvar_xml_structure_analyzer.py on your server,")
        print("then copy the output JSON here and re-run this script.")
    if not acmg_stats:
        print("\nNote: Figure 9 (ACMG criteria analysis) was skipped because")
        print("acmg_criteria_detailed.json was not found.")
        print("To generate it, run clinvar_acmg_extractor.py on your server,")
        print("then copy the output JSON here and re-run this script.")
    if concordance_data:
        print("\nNote: Concordance data loaded — Figure 3 includes panel (d).")
    else:
        print("\nNote: concordance_timeseries.json not found — Figure 3 uses 3-panel layout.")
        print("To generate it, run clinvar_concordance_analyzer.py on your server.")
    if domain_data:
        print("\nNote: Figure 10 (disease domain) + rare-vs-common supplement generated.")
    else:
        print("\nNote: domain_summary.json not found — Figure 10 skipped.")
        print("To generate it, run clinvar_domain_analyzer.py.")
    if driver_data:
        print("\nNote: Figure 11 (reclassification drivers) generated.")
    else:
        print("\nNote: reclassification_drivers.json not found — Figure 11 skipped.")
        print("To generate it, run clinvar_reclassification_driver.py.")


if __name__ == "__main__":
    main()
