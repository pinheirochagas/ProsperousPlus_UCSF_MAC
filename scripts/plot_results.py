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


def plot_residue_detail(run_dir: Path, out_path: Path):
    """Per-residue contribution bar chart for every unique candidate peptide.

    Each panel shows one candidate 8-mer: green bars = CathL contribution,
    orange bars = CathB contribution, stacked per S/T position.
    Non-S/T positions are shown as grey zero-height markers so the sequence
    context is visible.  Negative contributions are drawn below zero.
    """
    def load(rel):
        p = run_dir / rel
        if not p.exists():
            sys.exit(f"ERROR: missing file {p}")
        return pd.read_csv(p)

    cands   = load("step5_scores/biomarker_candidates.csv")
    residue = load("step5_scores/per_residue_contributions.csv")

    # ── Deduplicate: keep one representative window per unique mutation_sum
    # group within each sequence_id (overlapping windows share the same sum)
    cands_sorted = cands.sort_values("position")
    seen = set()
    unique_cands = []
    for _, row in cands_sorted.iterrows():
        key = (row["sequence_id"], round(row["mutation_sum"], 4))
        if key not in seen:
            seen.add(key)
            unique_cands.append(row)
    unique_cands = sorted(unique_cands,
                          key=lambda r: r["mutation_sum"], reverse=True)

    n = len(unique_cands)
    if n == 0:
        print("No candidates to plot.")
        return

    # Build a lookup: abs_position → (contrib_L, contrib_B)
    res_lookup = {
        row["abs_position"]: (row["contrib_L"], row["contrib_B"])
        for _, row in residue.iterrows()
    }

    # Shared y-scale across all panels
    ymax = max(0.1, max(a + max(0, b) for a, b in res_lookup.values()) * 1.15)
    ymin = min(0.0, min(min(a, 0) + min(b, 0) for a, b in res_lookup.values()) * 1.15)

    # Cap at top-10 candidates and arrange in at most 2 columns of 5
    import math
    MAX_CANDS = 10
    MAX_ROWS  = 5
    unique_cands = unique_cands[:MAX_CANDS]
    n        = len(unique_cands)
    n_cols   = math.ceil(n / MAX_ROWS)
    n_rows   = min(n, MAX_ROWS)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(10 * n_cols, 3.2 * n_rows),
        squeeze=False,
    )
    fig.patch.set_facecolor("white")
    fig.suptitle(
        "Per-residue S/T contribution to cleavage resistance",
        fontsize=13, fontweight="bold", y=1.01,
    )

    # Hide any unused axes (row-major order)
    for idx in range(n, n_rows * n_cols):
        r, c = idx // n_cols, idx % n_cols
        axes[r][c].set_visible(False)

    for rank, cand in enumerate(unique_cands, 1):
        row = (rank - 1) // n_cols   # left→right, top→bottom
        col = (rank - 1) % n_cols
        ax  = axes[row][col]

        peptide  = cand["peptide"]
        center   = int(cand["position"])   # cleavage-site position (1-indexed)
        seq_id   = cand["sequence_id"]
        mut_sum  = cand["mutation_sum"]
        # region label: last two underscore-delimited tokens, e.g. TAU_HUMAN_402-461
        region_label = "_".join(seq_id.split("_")[-3:]) if "_" in seq_id else seq_id

        # 8-mer spans center-3 … center+4  (window center = pos 4 of 8-mer)
        pep_start = center - 3
        xs = list(range(8))           # bar x-positions 0–7

        contrib_L = []
        contrib_B = []
        colors_x  = []
        for i, aa in enumerate(peptide):
            abs_pos = pep_start + i
            if abs_pos in res_lookup:
                cL, cB = res_lookup[abs_pos]
            else:
                cL, cB = 0.0, 0.0
            contrib_L.append(cL)
            contrib_B.append(cB)
            colors_x.append(C_CATHL if aa in ("S", "T") else "#cccccc")

        contrib_L = np.array(contrib_L)
        contrib_B = np.array(contrib_B)

        # Split positive and negative portions for correct stacking
        pos_L = np.clip(contrib_L, 0, None)
        neg_L = np.clip(contrib_L, None, 0)
        pos_B = np.clip(contrib_B, 0, None)
        neg_B = np.clip(contrib_B, None, 0)

        bar_w = 0.55

        # Positive stack: CathL (green) then CathB (orange) on top
        ax.bar(xs, pos_L, width=bar_w, color=C_CATHL, alpha=0.85, label="CathL (S/T→P)")
        ax.bar(xs, pos_B, width=bar_w, bottom=pos_L, color=C_CATHB, alpha=0.85, label="CathB (S/T→E)")

        # Negative stack: below zero
        ax.bar(xs, neg_L, width=bar_w, color=C_CATHL, alpha=0.55)
        ax.bar(xs, neg_B, width=bar_w, bottom=neg_L, color=C_CATHB, alpha=0.55)

        ax.axhline(0, color="#888888", lw=0.7)

        # X-axis: residue letters, S/T in bold colour
        ax.set_xticks(xs)
        labels = []
        for i, aa in enumerate(peptide):
            abs_pos = pep_start + i
            labels.append(f"{aa}\n{abs_pos}")
        ax.set_xticklabels(labels, fontsize=14)
        for tick, aa in zip(ax.get_xticklabels(), peptide):
            if aa in ("S", "T"):
                tick.set_color(C_ABOVE)
                tick.set_fontweight("bold")
            else:
                tick.set_color("#666666")

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(-0.6, 7.6)
        ax.set_ylabel("Contribution", fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=10)

        ax.set_title(
            f"#{rank}  {peptide}  |  center pos {center}  |  {region_label}  |  "
            f"mut sum = {mut_sum:.3f}",
            fontsize=10, fontweight="bold", loc="left", pad=4,
        )

        # Legend only on first panel
        if rank == 1:  # noqa: SIM102
            ax.legend(fontsize=10, loc="upper right", framealpha=0.7)

    fig.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {out_path}")


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

    # label the single top candidate only
    top1 = cands.nlargest(1, "mutation_sum").iloc[0]
    ax4.text(top1["position"], 0.82,
             f"{top1['peptide']} (pos {int(top1['position'])})",
             ha="center", va="bottom", fontsize=7,
             color=C_ABOVE, fontweight="bold")

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

    out_path    = Path(args.out) if args.out else run_dir / "pipeline_summary.png"
    detail_path = out_path.with_name("residue_detail.png")

    plot_pipeline(run_dir, out_path)
    plot_residue_detail(run_dir, detail_path)


if __name__ == "__main__":
    main()
