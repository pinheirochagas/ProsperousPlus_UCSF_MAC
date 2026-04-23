#!/usr/bin/env python3
"""
Genome-browser style visualization of a completed pipeline run.

All 4 tracks share the same tau-protein position x-axis so any position
can be traced from raw cleavage scores down to biomarker candidates.

Usage
-----
    conda activate prosperousplus
    python plot_results.py --run-dir results/tau_only_2026-04-23_2038
    # saves: results/tau_only_2026-04-23_2038/pipeline_summary.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


# ── Colours ───────────────────────────────────────────────────────────────────
C_CATHD    = "#d62728"    # red
C_LEGUMAIN = "#1f77b4"    # blue
C_CATHL    = "#2ca02c"    # green
C_CATHB    = "#ff7f0e"    # orange
C_REGION   = "steelblue"  # region bands
C_ABOVE    = "#e6550d"    # mutation stems above threshold
C_BELOW    = "#cccccc"    # mutation stems below threshold
C_THRESH   = "#888888"    # threshold dashed lines


def _shade_regions(ax, regions, ymin=0, ymax=1):
    """Draw low-cleavage region bands on an axis."""
    for _, r in regions.iterrows():
        ax.axvspan(r["start_position"], r["end_position"],
                   alpha=0.13, color=C_REGION, zorder=0, linewidth=0)


def _add_region_labels(ax, regions, y, fontsize=7):
    """Print region index above each band on one axis only."""
    for i, (_, r) in enumerate(regions.iterrows(), 1):
        mid = (r["start_position"] + r["end_position"]) / 2
        ax.text(mid, y, str(i), ha="center", va="bottom",
                fontsize=fontsize, color=C_REGION, fontweight="bold")


def plot_pipeline(run_dir: Path, out_path: Path):
    # ── Load data ──────────────────────────────────────────────────────────
    def load(rel):
        p = run_dir / rel
        if not p.exists():
            sys.exit(f"ERROR: missing file {p}\n(Is the run complete?)")
        return pd.read_csv(p)

    cathd    = load("step1_predictions/cathepsinD_cleavage.csv")
    legumain = load("step1_predictions/legumain_cleavage.csv")
    regions  = load("step2_regions/low_cleavage_regions.csv")
    scores   = load("step5_scores/mutation_scores.csv")
    cands    = load("step5_scores/biomarker_candidates.csv")

    prot_len = int(cathd["position"].max())
    n_cands  = len(cands)

    # ── Figure layout ──────────────────────────────────────────────────────
    fig = plt.figure(figsize=(17, 12))
    fig.patch.set_facecolor("white")

    gs = gridspec.GridSpec(
        4, 1,
        height_ratios=[2.5, 2.5, 1.5, 0.75],
        hspace=0.06,
        left=0.07, right=0.97, top=0.93, bottom=0.07,
    )
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)

    xmin, xmax = 1, prot_len

    # ── TRACK 1 — Cleavage scores ──────────────────────────────────────────
    _shade_regions(ax1, regions)
    ax1.plot(cathd["position"],    cathd["pro"],    color=C_CATHD,
             lw=0.9, alpha=0.85, label="Cathepsin D")
    ax1.plot(legumain["position"], legumain["pro"], color=C_LEGUMAIN,
             lw=0.9, alpha=0.85, label="Legumain")
    ax1.axhline(0.3, color=C_THRESH, lw=0.8, ls="--", label="Threshold 0.3")
    _add_region_labels(ax1, regions, y=1.02)
    ax1.set_ylim(0, 1.12)
    ax1.set_ylabel("Cleavage score", fontsize=9)
    ax1.set_title(
        f"Tau biomarker pipeline — {run_dir.name}   "
        f"({n_cands} candidates, {len(regions)} protected regions)",
        fontsize=11, fontweight="bold", pad=6,
    )
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.7)
    ax1.tick_params(labelbottom=False)
    ax1.set_xlim(xmin, xmax)

    # ── TRACK 2 — Mutation sums ────────────────────────────────────────────
    _shade_regions(ax2, regions)
    above = scores[scores["mutation_sum"] >= 0.475]
    below = scores[scores["mutation_sum"] <  0.475]

    # stems for below-threshold
    ax2.vlines(below["position"], 0, below["mutation_sum"],
               color=C_BELOW, lw=0.7, zorder=1)
    ax2.scatter(below["position"], below["mutation_sum"],
                color=C_BELOW, s=8, zorder=2)

    # stems for above-threshold (candidates)
    ax2.vlines(above["position"], 0, above["mutation_sum"],
               color=C_ABOVE, lw=1.2, zorder=3)
    ax2.scatter(above["position"], above["mutation_sum"],
                color=C_ABOVE, s=22, zorder=4, label="Candidates")

    ax2.axhline(0.475, color=C_THRESH, lw=0.8, ls="--", label="Threshold 0.475")
    ax2.set_ylabel("Mutation sum", fontsize=9)
    ax2.set_ylim(bottom=scores["mutation_sum"].min() - 0.05)
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.7)
    ax2.tick_params(labelbottom=False)

    # label top 6 candidates
    top = cands.nlargest(6, "mutation_sum")
    for _, c in top.iterrows():
        ax2.annotate(
            f"{c['peptide']}\n(pos {int(c['position'])})",
            xy=(c["position"], c["mutation_sum"]),
            xytext=(0, 8), textcoords="offset points",
            ha="center", fontsize=6.5, color=C_ABOVE,
            arrowprops=dict(arrowstyle="-", color=C_ABOVE, lw=0.5),
        )

    # ── TRACK 3 — CathL / CathB contributions ─────────────────────────────
    _shade_regions(ax3, regions)
    ax3.bar(scores["position"], scores["cathL_diff"].clip(lower=0),
            width=1.2, color=C_CATHL, alpha=0.8, label="CathL diff")
    ax3.bar(scores["position"], scores["cathB_diff"].clip(lower=0),
            width=1.2, bottom=scores["cathL_diff"].clip(lower=0),
            color=C_CATHB, alpha=0.8, label="CathB diff")
    ax3.axhline(0, color="#aaaaaa", lw=0.5)
    ax3.set_ylabel("Contribution", fontsize=9)
    ax3.legend(loc="upper right", fontsize=8, framealpha=0.7)
    ax3.tick_params(labelbottom=False)

    # ── TRACK 4 — Candidate markers ───────────────────────────────────────
    _shade_regions(ax4, regions)
    norm = mcolors.Normalize(
        vmin=cands["mutation_sum"].min(),
        vmax=cands["mutation_sum"].max(),
    )
    cmap = cm.get_cmap("YlOrRd")
    for _, c in cands.iterrows():
        color = cmap(norm(c["mutation_sum"]))
        ax4.scatter(c["position"], 0.5, marker="D", s=60,
                    color=color, zorder=5, edgecolors="white", linewidths=0.4)

    # label top 5 only to avoid clutter
    for _, c in cands.nlargest(5, "mutation_sum").iterrows():
        ax4.text(c["position"], 0.82, c["peptide"],
                 ha="center", va="bottom", fontsize=6,
                 color="black", rotation=45)

    ax4.set_ylim(0, 1.5)
    ax4.set_yticks([])
    ax4.set_ylabel("Candidates", fontsize=9)
    ax4.set_xlabel("Tau protein position (aa)", fontsize=10)

    # colorbar for candidate diamonds
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax4, orientation="vertical",
                        fraction=0.015, pad=0.01, aspect=8)
    cbar.set_label("Mut sum", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # ── Shared formatting ──────────────────────────────────────────────────
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim(xmin, xmax)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)

    # region legend on ax1
    region_patch = plt.matplotlib.patches.Patch(
        color=C_REGION, alpha=0.3, label="Low-cleavage region"
    )
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles + [region_patch], labels + ["Low-cleavage region"],
               loc="upper right", fontsize=8, framealpha=0.7)

    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Genome-browser plot of a ProsperousPlus pipeline run"
    )
    parser.add_argument(
        "--run-dir", required=True,
        help="Path to a completed run directory (e.g. results/tau_only_2026-04-23_2038)",
    )
    parser.add_argument(
        "--out", default=None,
        help="Output PNG path (default: <run-dir>/pipeline_summary.png)",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        sys.exit(f"ERROR: run directory not found: {run_dir}")

    out_path = Path(args.out) if args.out else run_dir / "pipeline_summary.png"
    plot_pipeline(run_dir, out_path)


if __name__ == "__main__":
    main()
