#!/usr/bin/env python3
"""
Generate a supplementary table and figure showing conflict severity stratification.

Categorizes reclassification events into:
  - Medically significant: P/LP <-> VUS/B/LB (clinical actionability changes)
  - Non-medically significant: VUS <-> B/LB
  - Adjacent (within-tier): P <-> LP, B <-> LB
  - Other: transitions involving Conflicting, P/LP, Other

Uses reclassification_events_detail.csv from the data directory.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pathlib import Path
import json

# ── Paths ──────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
CSV_PATH = DATA_DIR / "reclassification_events_detail.csv"

# ── Classification tiers (Larry Babb's 3-range framework) ─────────────
SUPPORTING = {"P", "LP"}         # Pathogenic / Likely Pathogenic
UNCERTAIN  = {"VUS"}             # Variant of Uncertain Significance
NOT_SUPPORTING = {"B", "LB"}    # Benign / Likely Benign
COMPOSITE  = {"Conflicting", "P/LP", "B/LB", "Other"}

def classify_transition(prev, curr):
    """Classify a single reclassification event by clinical severity."""
    # Skip if same category
    if prev == curr:
        return None

    # Adjacent within-tier transitions
    if {prev, curr} == {"P", "LP"}:
        return "Adjacent (P↔LP)"
    if {prev, curr} == {"B", "LB"}:
        return "Adjacent (B↔LB)"

    # Medically significant: crosses actionability boundary
    if (prev in SUPPORTING and curr in (UNCERTAIN | NOT_SUPPORTING)) or \
       (curr in SUPPORTING and prev in (UNCERTAIN | NOT_SUPPORTING)):
        return "Medically significant"

    # Non-medically significant: VUS <-> B/LB
    if (prev in UNCERTAIN and curr in NOT_SUPPORTING) or \
       (curr in UNCERTAIN and prev in NOT_SUPPORTING):
        return "Non-medically significant"

    # Involves composite categories (Conflicting, P/LP, B/LB, Other)
    if prev in COMPOSITE or curr in COMPOSITE:
        return "Involves composite category"

    return "Other"


def main():
    print("Loading reclassification events...")
    df = pd.read_csv(CSV_PATH)
    print(f"  Total events: {len(df):,}")

    # Classify each transition
    df["severity"] = df.apply(
        lambda r: classify_transition(r["Prev_Category"], r["Curr_Category"]),
        axis=1
    )

    # Remove null (same-category, shouldn't exist but just in case)
    df = df.dropna(subset=["severity"])

    # ── Summary Table ──────────────────────────────────────────────────
    summary = df.groupby("severity").agg(
        count=("Variation_ID", "size"),
        unique_variants=("Variation_ID", "nunique"),
    ).reset_index()
    summary["pct_of_events"] = (summary["count"] / summary["count"].sum() * 100).round(1)
    summary = summary.sort_values("count", ascending=False)

    print("\n=== CONFLICT SEVERITY STRATIFICATION ===")
    print(summary.to_string(index=False))

    # Save table
    summary.to_csv(DATA_DIR / "conflict_severity_stratification.csv", index=False)
    print(f"\nSaved: {DATA_DIR / 'conflict_severity_stratification.csv'}")

    # ── Detailed transition pairs within each severity ─────────────────
    detail = df.groupby(["severity", "Prev_Category", "Curr_Category"]).size() \
               .reset_index(name="count") \
               .sort_values(["severity", "count"], ascending=[True, False])
    detail.to_csv(DATA_DIR / "conflict_severity_detail.csv", index=False)
    print(f"Saved: {DATA_DIR / 'conflict_severity_detail.csv'}")

    # ── By mechanism ───────────────────────────────────────────────────
    mech = df.groupby(["severity", "Driver_Label"]).size().reset_index(name="count")
    mech_pivot = mech.pivot_table(
        index="severity", columns="Driver_Label", values="count", fill_value=0
    )
    mech_pivot.to_csv(DATA_DIR / "conflict_severity_by_mechanism.csv")
    print(f"Saved: {DATA_DIR / 'conflict_severity_by_mechanism.csv'}")

    # ── Figure ─────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(7.08, 3.0))
    fig.subplots_adjust(wspace=0.45)

    # Nature-style colors
    colors = {
        "Medically significant": "#E64B35",      # Red
        "Non-medically significant": "#F39B7F",   # Salmon
        "Adjacent (P↔LP)": "#4DBBD5",            # Cyan
        "Adjacent (B↔LB)": "#91D1C2",            # Teal
        "Involves composite category": "#B09C85",  # Tan
        "Other": "#CCCCCC",
    }

    # Panel (a): Bar chart of severity categories
    ax = axes[0]
    severity_order = [
        "Medically significant",
        "Non-medically significant",
        "Adjacent (P↔LP)",
        "Adjacent (B↔LB)",
        "Involves composite category",
        "Other",
    ]
    plot_data = summary.set_index("severity").reindex(severity_order).dropna()
    bar_colors = [colors.get(s, "#CCCCCC") for s in plot_data.index]

    y_pos = np.arange(len(plot_data))
    bars = ax.barh(y_pos, plot_data["count"], color=bar_colors,
                   edgecolor="white", linewidth=0.3, height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_data.index, fontsize=6)
    ax.set_xlabel("Number of reclassification events", fontsize=7)
    ax.set_title("(a) Reclassification events\nby clinical severity", fontsize=8)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v/1000:.0f}K" if v >= 1000 else f"{v:.0f}"
    ))

    for i, (c, p) in enumerate(zip(plot_data["count"], plot_data["pct_of_events"])):
        ax.text(c + max(plot_data["count"]) * 0.02, i,
                f"{c:,.0f} ({p}%)", va="center", fontsize=5, color="#666666")
    ax.set_xlim(0, max(plot_data["count"]) * 1.35)

    # Panel (b): Top transition pairs for medically significant only
    ax = axes[1]
    med_sig = detail[detail["severity"] == "Medically significant"].head(10)
    if len(med_sig) > 0:
        labels = [f"{r['Prev_Category']}→{r['Curr_Category']}" for _, r in med_sig.iterrows()]
        y_pos = np.arange(len(med_sig))
        ax.barh(y_pos, med_sig["count"].values, color="#E64B35",
                edgecolor="white", linewidth=0.3, height=0.6)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_xlabel("Number of events", fontsize=7)
        ax.set_title("(b) Top medically significant\ntransition pairs", fontsize=8)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(
            lambda v, _: f"{v/1000:.0f}K" if v >= 1000 else f"{v:.0f}"
        ))
        for i, c in enumerate(med_sig["count"].values):
            ax.text(c + max(med_sig["count"].values) * 0.02, i,
                    f"{c:,}", va="center", fontsize=5, color="#666666")
        ax.set_xlim(0, max(med_sig["count"].values) * 1.25)

    plt.tight_layout()

    # Save
    for ext in ["png", "pdf"]:
        fig.savefig(FIG_DIR / f"figure_supp_conflict_severity.{ext}",
                    dpi=600, bbox_inches="tight")
    print(f"\nSaved figure: {FIG_DIR / 'figure_supp_conflict_severity.png'}")

    plt.close()

    # ── Print key stats for manuscript ─────────────────────────────────
    print("\n=== KEY STATS FOR MANUSCRIPT ===")
    total = summary["count"].sum()
    for _, row in summary.iterrows():
        print(f"  {row['severity']}: {row['count']:,.0f} events ({row['pct_of_events']}%), "
              f"{row['unique_variants']:,.0f} unique variants")
    print(f"  Total: {total:,} events")


if __name__ == "__main__":
    main()
