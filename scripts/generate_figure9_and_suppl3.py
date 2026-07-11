#!/usr/bin/env python3
"""
Generate:
  - Main manuscript Figure 9 = panel (d) of figure3 (Inter-submitter concordance)
  - Supplementary Figure S3 = panel (a) of figure3 (Review status evolution)

Both as standalone single-panel figures, 600 DPI TIFF + PDF + PNG.
Output directly into the Dropbox submission package folders.
"""

import sys
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches

# Import shared visualization constants from the main script
sys.path.insert(0, str(Path(__file__).parent))
from clinvar_visualizations import (
    NAT, REVIEW_ORDER, REVIEW_COLORS, panel_label, date_to_float,
)

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT = Path(__file__).resolve().parent.parent
CACHE_DIR = PROJECT / "figures" / ".cache"
DATA_DIR = PROJECT / "data"

SUBMISSION = Path.home() / "Library/CloudStorage/Dropbox/ClinVar_Genome_Medicine_Submission"
FIG_OUT = SUBMISSION / "figures"
SUPP_OUT = SUBMISSION / "additional_files"

# ── Load cached history and concordance data ───────────────────────────
with open(CACHE_DIR / "history_agg.json") as f:
    hist = json.load(f)
# history_agg is already a dict keyed by snapshot date (e.g. "2012-11")

# Concordance data (VCV-era concordance/discordance rates)
concordance_path = None
for candidate in [
    DATA_DIR / "concordance_timeseries_v2.json",
    DATA_DIR / "concordance_timeseries_v1.json",
]:
    if candidate.exists():
        concordance_path = candidate
        break

concordance_data = None
if concordance_path:
    with open(concordance_path) as f:
        concordance_data = json.load(f)


def _save_all(fig, out_dir: Path, name: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    # TIFF with LZW compression (BMC requirement) + PDF + PNG
    fig.savefig(out_dir / f"{name}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.png", dpi=600, bbox_inches="tight")
    print(f"Saved {name}.tiff / .pdf / .png to {out_dir}")


# ══════════════════════════════════════════════════════════════════════
# MAIN FIGURE 9 — Inter-submitter concordance (formerly figure3 panel d)
# ══════════════════════════════════════════════════════════════════════

def make_figure_9():
    fig, ax = plt.subplots(figsize=(3.5, 2.8))

    if not concordance_data or not concordance_data.get("snapshots"):
        print("ERROR: concordance data not available; cannot generate Figure 9")
        return

    conc_snaps = sorted(concordance_data["snapshots"], key=lambda s: s["date"])
    vcv_start = "2020-01"
    vcv_snaps = [s for s in conc_snaps if s["date"] >= vcv_start]
    all_snaps = conc_snaps

    all_x = np.array([date_to_float(s["date"]) for s in all_snaps])
    all_multi = np.array([s.get("multi_sub_variants", 0) for s in all_snaps])

    if vcv_snaps:
        v_dates = [s["date"] for s in vcv_snaps]
        v_x = np.array([date_to_float(d) for d in v_dates])
        v_conc = np.array([s.get("concordance_rate_pct", 0) for s in vcv_snaps])
        v_disc = np.array([s.get("discordance_rate_pct", 0) for s in vcv_snaps])

        ax.fill_between(v_x, 0, v_conc, color=NAT["lime"], alpha=0.55,
                        label="Concordant")
        ax.fill_between(v_x, v_conc, v_conc + v_disc,
                        color=NAT["dark_red"], alpha=0.55, label="Discordant")
        ax.fill_between(v_x, v_conc + v_disc, 100,
                        color=NAT["lt_grey"], alpha=0.35, label="Excluded")

    # Grey RCV era
    rcv_end_x = date_to_float("2019-12")
    ax.set_xlim(date_to_float("2012-11"), date_to_float("2025-07"))
    ax.axvspan(date_to_float("2012-11"), rcv_end_x,
               color="#eeeeee", alpha=0.8, zorder=0)
    ax.text(date_to_float("2016-06"), 50,
            "RCV era:\nconcordance not\nmeasurable from\nReview Status",
            ha="center", va="center", fontsize=6, color=NAT["grey"],
            style="italic")
    ax.axvline(x=date_to_float("2020-01"), color="black", linestyle="--",
               linewidth=0.6, alpha=0.4)

    ax.set_ylim(0, 100)
    ax.set_ylabel("% of multi-submitter variants")
    ax.set_xlabel("Year")
    ax.set_title("Inter-submitter concordance")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    if vcv_snaps:
        last = vcv_snaps[-1]
        ax.text(0.98, 0.05,
                f"Latest ({last['date']}): "
                f"{last.get('concordance_rate_pct', 0):.1f}% concordant, "
                f"{last.get('discordance_rate_pct', 0):.1f}% discordant",
                transform=ax.transAxes, fontsize=5, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec=NAT["lt_grey"], lw=0.4, alpha=0.9))

    # Secondary axis: multi-sub variants (×1K) across full timeline
    ax2 = ax.twinx()
    ax2.plot(all_x, all_multi / 1e3, "k:", linewidth=0.8, alpha=0.4,
             label="Multi-sub variants (K)")
    ax2.set_ylabel("Multi-sub variants (\u00d71K)", fontsize=7)
    ax2.tick_params(labelsize=6)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_linewidth(0.5)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=5,
              loc="upper left",
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)

    plt.tight_layout()
    _save_all(fig, FIG_OUT, "Figure_9")
    plt.close()


# ══════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY FIGURE S3 — Review status evolution (formerly figure3 panel a)
# ══════════════════════════════════════════════════════════════════════

def make_supplementary_s3():
    fig, ax = plt.subplots(figsize=(3.5, 2.8))

    dates = sorted(hist.keys())
    x = np.array([date_to_float(d) for d in dates])

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
    ax.legend(handles=handles, fontsize=5, loc="upper left",
              bbox_to_anchor=(0.02, 0.98), ncol=1, handlelength=1.0,
              handletextpad=0.3, borderpad=0.2, frameon=True,
              facecolor="white", edgecolor="none", framealpha=0.85)

    plt.tight_layout()
    _save_all(fig, SUPP_OUT, "Additional_file_3_Supplementary_Figure_S3")
    plt.close()


if __name__ == "__main__":
    make_figure_9()
    make_supplementary_s3()
    print("\nDone.")
