#!/usr/bin/env python3
"""
Fast cohort-level postprocessing and visualization.

Reads only the merged cohort files:
  - funnel.csv
  - summary.csv
  - master_biomarker_candidates.csv
  - step6_ptm/annotated_summary.csv  (auto-detected; override with --ptm)

Optionally reads the input FASTA (--fasta) for exact sequence lengths.
If not provided, sequence length is estimated from the max region-end
coordinate in master_biomarker_candidates.csv.

Usage
-----
    python plot_cohort.py \\
        --cohort-dir results/cohort_tissue_category_rna_brain_Tissue_Tissue_enriched_2026-04-24_0525 \\
        [--fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta] \\
        [--ptm results/.../step6_ptm/annotated_summary.csv]

Outputs
-------
    <cohort-dir>/cohort_publication_summary.png
    <cohort-dir>/cohort_exploration_detail.png
    <cohort-dir>/cohort_summary.txt
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


THRESHOLD = 0.475

# Restrained publication palette
BLACK = "#111111"
GREY_900 = "#262626"
GREY_700 = "#525252"
GREY_500 = "#737373"
GREY_300 = "#D4D4D4"
GREY_200 = "#E5E5E5"
GREY_100 = "#F5F5F5"
BLUE = "#2B6CB0"
RED = "#B23A48"
GREEN = "#2F6F4E"
ORANGE = "#B56A2A"

# Evidence strength palette (0 = no phospho hit, 1–3 = database count)
PHOSPHO_COLORS = {0: "#D4D4D4", 1: "#F59E0B", 2: "#EF4444", 3: "#7C3AED"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: missing required file: {path}")
    return pd.read_csv(path)


def _as_optional_int(value) -> Optional[int]:
    if pd.isna(value):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _read_funnel(cohort_dir: Path) -> Tuple[Dict[str, Optional[int]], List[str]]:
    funnel_path = cohort_dir / "funnel.csv"
    if not funnel_path.exists():
        return {
            "initial": None,
            "with_regions": None,
            "with_candidates": None,
            "total_regions": None,
            "total_residues": None,
        }, [
            "WARNING: funnel.csv missing — run run_cohort.py --merge-only to generate it."
        ]
    df = pd.read_csv(funnel_path)
    counts = dict(zip(df["stage"], df["count"]))
    return {
        "initial": _as_optional_int(counts.get("Initial proteins")),
        "with_regions": _as_optional_int(counts.get("Proteins with regions (step 2)")),
        "with_candidates": _as_optional_int(counts.get("Proteins with candidates")),
        "total_regions": _as_optional_int(counts.get("Total candidate regions")),
        "total_residues": _as_optional_int(counts.get("Total candidate residues")),
    }, []


def _load_seq_lengths(fasta_path: Optional[Path], master: pd.DataFrame) -> Dict[str, int]:
    """Return {uniprot_id: seq_len}. Uses FASTA parse if available, else region-end proxy."""
    if fasta_path is not None and fasta_path.exists():
        lengths: Dict[str, int] = {}
        current_id: Optional[str] = None
        current_len = 0
        with open(str(fasta_path)) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if current_id is not None:
                        lengths[current_id] = current_len
                    parts = line[1:].split("|")
                    current_id = parts[1] if len(parts) >= 2 else parts[0].split()[0]
                    current_len = 0
                else:
                    current_len += len(line)
        if current_id is not None:
            lengths[current_id] = current_len
        return lengths

    # Fallback: infer from max region-end coordinate per protein
    def _region_end(region_str: str) -> int:
        try:
            return int(str(region_str).split("-")[-1])
        except (ValueError, AttributeError):
            return 0

    proxy = master.copy()
    proxy["_rend"] = proxy["region"].apply(_region_end)
    return proxy.groupby("uniprot_id")["_rend"].max().to_dict()


def _load_ptm(ptm_path: Optional[Path], warn_missing: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load annotated_summary.csv. Returns (per_residue_df, per_protein_df).

    Empty DataFrames are returned silently when ptm_path is None, or with a
    warning when the file was explicitly requested but not found.
    """
    _cols_res = ["uniprot_id", "position", "aa", "contrib_total",
                 "n_sources", "sources", "n_pubmed", "pubmed_ids"]
    _cols_prot = ["uniprot_id", "n_phospho_residues", "max_n_sources", "max_n_pubmed"]
    _empty_res = pd.DataFrame(columns=_cols_res)
    _empty_prot = pd.DataFrame(columns=_cols_prot)

    if ptm_path is None:
        return _empty_res, _empty_prot
    if not ptm_path.exists():
        if warn_missing:
            print(f"WARNING: PTM file not found: {ptm_path}")
        return _empty_res, _empty_prot

    ptm = pd.read_csv(ptm_path)
    per_protein = ptm.groupby("uniprot_id").agg(
        n_phospho_residues=("position", "count"),
        max_n_sources=("n_sources", "max"),
        max_n_pubmed=("n_pubmed", "max"),
    ).reset_index()
    return ptm, per_protein


def load_cohort_data(
    cohort_dir: Path,
    fasta_path: Optional[Path] = None,
    ptm_path: Optional[Path] = None,
    warn_ptm_missing: bool = False,
):
    """Load all cohort data; enrich summary and master with seq_len and PTM columns."""
    master = _load_csv(cohort_dir / "master_biomarker_candidates.csv")
    summary = _load_csv(cohort_dir / "summary.csv")
    funnel, warnings = _read_funnel(cohort_dir)

    if funnel["with_candidates"] is None:
        funnel["with_candidates"] = int(summary["uniprot_id"].nunique())
    if funnel["total_regions"] is None:
        funnel["total_regions"] = int(summary["n_regions"].sum())
    if funnel["total_residues"] is None:
        funnel["total_residues"] = int(len(master))

    # Sequence lengths (FASTA or region-end proxy)
    seq_len_dict = _load_seq_lengths(fasta_path, master)
    summary["seq_len"] = summary["uniprot_id"].map(seq_len_dict).fillna(0).astype(int)

    # PTM annotations
    ptm_res, ptm_prot = _load_ptm(ptm_path, warn_missing=warn_ptm_missing)

    # Join PTM per-protein stats onto summary
    if not ptm_prot.empty:
        summary = summary.merge(
            ptm_prot[["uniprot_id", "n_phospho_residues", "max_n_sources", "max_n_pubmed"]],
            on="uniprot_id",
            how="left",
        )
        summary["n_phospho_residues"] = summary["n_phospho_residues"].fillna(0).astype(int)
        summary["max_n_sources"] = summary["max_n_sources"].fillna(0).astype(int)
        summary["max_n_pubmed"] = summary["max_n_pubmed"].fillna(0).astype(int)
    else:
        summary["n_phospho_residues"] = 0
        summary["max_n_sources"] = 0
        summary["max_n_pubmed"] = 0

    # Mark per-residue phospho status in master
    if not ptm_res.empty:
        phospho_keys = set(
            zip(ptm_res["uniprot_id"].astype(str), ptm_res["position"].astype(str))
        )
        master["is_phospho"] = [
            (str(u), str(p)) in phospho_keys
            for u, p in zip(master["uniprot_id"], master["position"])
        ]
    else:
        master["is_phospho"] = False

    # Add protein_name to ptm_res for Panel H labels (column may already exist in the CSV)
    if not ptm_res.empty:
        if "protein_name" not in ptm_res.columns:
            name_map = summary[["uniprot_id", "protein_name"]].drop_duplicates()
            ptm_res = ptm_res.merge(name_map, on="uniprot_id", how="left")
        ptm_res["protein_name"] = ptm_res["protein_name"].fillna(ptm_res["uniprot_id"])

    return master, summary, funnel, warnings, ptm_res, ptm_prot


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _count_text(value: Optional[int]) -> str:
    return "n/a" if value is None else f"{value:,}"


def _pct(value, total: Optional[int]) -> str:
    if value is None or total in (None, 0):
        return "n/a"
    return f"{int(value) / int(total) * 100:.1f}%"


def _format_name(name: str) -> str:
    return str(name).replace("_HUMAN", "")


def _candidate_cmap():
    try:
        return cm.get_cmap("mako")
    except ValueError:
        return cm.get_cmap("viridis")


def _style_axis(ax, grid: Optional[str] = None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GREY_300)
    ax.spines["bottom"].set_color(GREY_300)
    ax.tick_params(colors=GREY_700, length=2.5, width=0.6)
    ax.xaxis.label.set_color(GREY_700)
    ax.yaxis.label.set_color(GREY_700)
    if grid == "x":
        ax.xaxis.grid(True, color=GREY_200, lw=0.5, zorder=0)
    elif grid == "y":
        ax.yaxis.grid(True, color=GREY_200, lw=0.5, zorder=0)
    elif grid == "both":
        ax.grid(True, color=GREY_200, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def _panel_letter(ax, letter: str, label: str):
    ax.text(
        -0.12, 1.06, letter, transform=ax.transAxes,
        ha="left", va="bottom", fontsize=11, fontweight="bold", color=BLACK,
    )
    ax.set_title(label, loc="left", fontsize=9, fontweight="bold", color=BLACK, pad=6)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def build_report_lines(
    cohort_dir: Path,
    master: pd.DataFrame,
    summary: pd.DataFrame,
    funnel: Dict[str, Optional[int]],
    warnings: List[str],
    ptm_res: Optional[pd.DataFrame] = None,
) -> List[str]:
    initial = funnel["initial"]
    lines = [
        "=" * 72,
        f"Cohort summary: {cohort_dir.name}",
        "=" * 72,
    ]
    if warnings:
        lines.extend(warnings)
        lines.append("")

    lines.extend([
        "Selection funnel",
        f"  Initial proteins             : {_count_text(initial)}",
        "  Proteins with regions        : "
        f"{_count_text(funnel['with_regions'])}  ({_pct(funnel['with_regions'], initial)})",
        "  Proteins with candidates     : "
        f"{_count_text(funnel['with_candidates'])}  ({_pct(funnel['with_candidates'], initial)})",
        f"  Total candidate regions      : {_count_text(funnel['total_regions'])}",
        f"  Total candidate residues     : {_count_text(funnel['total_residues'])}",
    ])

    if ptm_res is not None and not ptm_res.empty:
        n_prots = int(ptm_res["uniprot_id"].nunique())
        n_multi = int((ptm_res.groupby("uniprot_id")["n_sources"].max() >= 2).sum())
        with_cands = funnel["with_candidates"]
        lines.extend([
            f"  Proteins with phospho hit    : {n_prots}  ({_pct(n_prots, with_cands)} of candidates)",
            f"  Proteins with ≥2 DB support  : {n_multi}",
        ])

    lines.extend([
        "",
        "Top 10 proteins by strongest candidate residue",
        f"  {'#':>3}  {'Protein':<18}  {'UniProt':>10}  "
        f"{'top_total':>9}  {'n_cands':>7}  {'n_regions':>9}",
        "  " + "-" * 66,
    ])

    top = summary.sort_values("top_contrib_total", ascending=False).head(10)
    for rank, (_, row) in enumerate(top.iterrows(), 1):
        lines.append(
            f"  {rank:>3}  {row['protein_name']:<18}  {row['uniprot_id']:>10}  "
            f"{row['top_contrib_total']:>9.4f}  {int(row['n_candidates']):>7}  "
            f"{int(row['n_regions']):>9}"
        )

    lines.extend([
        "",
        "Residue-level candidate composition",
        f"  Serine candidates            : {(master['aa'] == 'S').sum():,}",
        f"  Threonine candidates         : {(master['aa'] == 'T').sum():,}",
        f"  Median contrib_total         : {master['contrib_total'].median():.4f}",
        f"  Max contrib_total            : {master['contrib_total'].max():.4f}",
    ])

    if ptm_res is not None and not ptm_res.empty:
        lines.extend([
            "",
            "Phospho-annotation detail (step 6)",
            f"  Total annotated residues     : {len(ptm_res):,}",
            f"  Supported by 3 databases     : {(ptm_res['n_sources'] == 3).sum():,}",
            f"  Supported by 2 databases     : {(ptm_res['n_sources'] == 2).sum():,}",
            f"  Supported by 1 database      : {(ptm_res['n_sources'] == 1).sum():,}",
            "",
            "Top 10 by phospho evidence (n_sources → contrib_total)",
            f"  {'#':>3}  {'Residue':<24}  {'contrib':>7}  {'n_src':>6}  {'n_pubs':>7}  sources",
            "  " + "-" * 72,
        ])
        top_ptm = ptm_res.sort_values(
            ["n_sources", "contrib_total"], ascending=[False, False]
        ).head(10)
        for rank, (_, row) in enumerate(top_ptm.iterrows(), 1):
            pname = _format_name(str(row.get("protein_name", row["uniprot_id"])))
            residue_lbl = f"{pname} {row['aa']}{int(row['position'])}"
            lines.append(
                f"  {rank:>3}  {residue_lbl:<24}  {row['contrib_total']:>7.4f}  "
                f"{int(row['n_sources']):>6}  {int(row['n_pubmed']):>7}  {row['sources']}"
            )

    lines.append("=" * 72)
    return lines


def save_and_print_report(lines: List[str], report_path: Path):
    text = "\n".join(lines) + "\n"
    print(text)
    report_path.write_text(text)
    print(f"Saved: {report_path}")


# ---------------------------------------------------------------------------
# Publication figure panels
# ---------------------------------------------------------------------------

def _pub_funnel(ax, funnel: Dict[str, Optional[int]], ptm_res: pd.DataFrame):
    """Panel A — selection funnel, extended to include step-6 phospho stages."""
    n_initial = funnel["initial"] or 0
    counts = [n_initial, funnel["with_regions"] or 0, funnel["with_candidates"] or 0]
    labels = ["Input\nFASTA", "With\nregions", "Candidates\n(step 5)"]

    if not ptm_res.empty:
        n_phospho_prots = int(ptm_res["uniprot_id"].nunique())
        n_multi = int((ptm_res.groupby("uniprot_id")["n_sources"].max() >= 2).sum())
        counts += [n_phospho_prots, n_multi]
        labels += ["Phospho-\nannotated", "Multi-DB\nevidence"]

    x = np.arange(len(labels))
    stage_colors = [GREY_500, BLUE, ORANGE, GREEN, RED][: len(labels)]

    ax.plot(x, counts, color=GREY_300, lw=1.2, zorder=1)
    ax.scatter(x, counts, s=160, color=stage_colors, edgecolor="white", linewidth=1.2, zorder=3)

    ymax = max(counts) if counts else 1
    for xi, count, _color in zip(x, counts, stage_colors):
        pct = _pct(count, n_initial)
        ax.text(xi, count + ymax * 0.09, f"{count:,}",
                ha="center", va="bottom", fontsize=9, fontweight="bold", color=BLACK)
        if pct != "n/a" and xi > 0:
            ax.text(xi, count + ymax * 0.025, f"({pct})",
                    ha="center", va="bottom", fontsize=6.5, color=GREY_700)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7, color=GREY_700)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(0, ymax * 1.35)
    ax.set_ylabel("# proteins")
    _panel_letter(ax, "A", "Selection funnel")
    _style_axis(ax, grid="y")


def _pub_score_distribution(ax, master: pd.DataFrame, ptm_res: pd.DataFrame):
    """Panel B — contrib_total histogram split by S/T with phospho overlay."""
    all_scores = master["contrib_total"]
    bins = np.linspace(THRESHOLD, all_scores.max() * 1.02, 28)

    s_all = master.loc[master["aa"] == "S", "contrib_total"]
    t_all = master.loc[master["aa"] == "T", "contrib_total"]

    ax.hist(s_all, bins=bins, color=GREEN, alpha=0.20, zorder=1)
    ax.hist(t_all, bins=bins, color=ORANGE, alpha=0.20, zorder=1)

    ps_count = int(((master["aa"] == "S") & master["is_phospho"]).sum())
    pt_count = int(((master["aa"] == "T") & master["is_phospho"]).sum())

    if ps_count > 0:
        ps_scores = master.loc[(master["aa"] == "S") & master["is_phospho"], "contrib_total"]
        ax.hist(ps_scores, bins=bins, color=GREEN, alpha=0.82, zorder=2)
    if pt_count > 0:
        pt_scores = master.loc[(master["aa"] == "T") & master["is_phospho"], "contrib_total"]
        ax.hist(pt_scores, bins=bins, color=ORANGE, alpha=0.82, zorder=2)

    legend_handles = [
        Patch(facecolor=GREEN, alpha=0.25, label=f"Ser all (n={len(s_all)})"),
        Patch(facecolor=ORANGE, alpha=0.25, label=f"Thr all (n={len(t_all)})"),
    ]
    if ps_count or pt_count:
        legend_handles += [
            Patch(facecolor=GREEN, alpha=0.85, label=f"Ser + phospho (n={ps_count})"),
            Patch(facecolor=ORANGE, alpha=0.85, label=f"Thr + phospho (n={pt_count})"),
        ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=6, loc="upper right")

    ax.axvline(THRESHOLD, color=GREY_500, lw=0.9, ls=(0, (3, 3)))
    ax.text(0.03, 0.97,
            f"median {all_scores.median():.3f}  ·  max {all_scores.max():.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.5, color=GREY_700)
    ax.set_xlabel("contrib_total")
    ax.set_ylabel("# residues")
    _panel_letter(ax, "B", "Score distribution")
    _style_axis(ax, grid="y")


def _pub_ranked_strength(ax, summary: pd.DataFrame):
    """Panel C — all proteins ranked by top_contrib_total, dots colored by phospho evidence."""
    df = summary.sort_values("top_contrib_total", ascending=False).reset_index(drop=True)
    x = np.arange(1, len(df) + 1)
    vals = df["top_contrib_total"].values
    src = df["max_n_sources"].values

    ax.fill_between(x, THRESHOLD, vals, color=BLUE, alpha=0.07, zorder=0)
    ax.plot(x, vals, color=GREY_300, lw=0.9, zorder=2)

    for n_src in [0, 1, 2, 3]:
        mask = src == n_src
        if not mask.any():
            continue
        color = PHOSPHO_COLORS[n_src]
        size = 10 if n_src > 0 else 6
        alpha = 0.90 if n_src > 0 else 0.50
        lbl = (f"{n_src} DB{'s' if n_src > 1 else ''}" if n_src > 0 else "no phospho")
        ax.scatter(x[mask], vals[mask], s=size, color=color, alpha=alpha,
                   edgecolors="none", zorder=3 + n_src, label=lbl)

    ax.axhline(THRESHOLD, color=GREY_500, lw=0.8, ls=(0, (3, 3)))
    ax.text(len(df) * 0.97, THRESHOLD + 0.012, f"threshold {THRESHOLD}",
            ha="right", va="bottom", fontsize=6.5, color=GREY_700)

    top1 = df.iloc[0]
    ax.annotate(
        _format_name(top1["protein_name"]),
        xy=(1, top1["top_contrib_total"]),
        xytext=(18, -4), textcoords="offset points",
        fontsize=7, color=RED, fontweight="bold",
        arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY_500),
    )

    top5_lines = []
    for i, row in df.head(5).iterrows():
        mark = "●" if row["max_n_sources"] > 0 else " "
        top5_lines.append(
            f"  {i + 1}. {mark} {_format_name(row['protein_name'])}  {row['top_contrib_total']:.3f}"
        )
    ax.text(0.97, 0.96, "Top 5\n" + "\n".join(top5_lines),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=6, color=GREY_700, family="monospace",
            bbox=dict(facecolor="white", edgecolor=GREY_200, linewidth=0.6, pad=3))

    ax.legend(frameon=False, fontsize=6.5, loc="upper right", markerscale=1.5,
              bbox_to_anchor=(0.97, 0.68))
    ax.set_xlim(0, len(df) + 3)
    ax.set_ylim(0.44, max(vals.max() * 1.08, 1.06))
    ax.set_xlabel("Protein rank")
    ax.set_ylabel("Top contrib_total")
    _panel_letter(ax, "C", "Candidate strength — all proteins")
    _style_axis(ax, grid="y")


def _pub_contribution_balance(ax, master: pd.DataFrame):
    """Panel D — CathL vs CathB scatter, axes clipped to 99th percentile."""
    xlo = -0.05
    xhi = master["contrib_L"].quantile(0.99) * 1.08
    ylo = -0.05
    yhi = master["contrib_B"].quantile(0.99) * 1.08

    for aa, color, label in [("S", GREEN, "Serine"), ("T", ORANGE, "Threonine")]:
        sub = master[master["aa"] == aa]
        ax.scatter(sub["contrib_L"], sub["contrib_B"], s=9, color=color, alpha=0.35,
                   edgecolors="none", label=f"{label} (n={len(sub):,})", rasterized=True)

    xs = np.linspace(xlo, xhi, 200)
    ax.plot(xs, THRESHOLD - xs, color=GREY_500, lw=0.85, ls=(0, (3, 3)),
            label=f"L+B = {THRESHOLD}")
    ax.axhline(0, color=GREY_200, lw=0.6)
    ax.axvline(0, color=GREY_200, lw=0.6)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel("CathL contribution")
    ax.set_ylabel("CathB contribution")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    _panel_letter(ax, "D", "Protease contribution balance")
    _style_axis(ax, grid="both")


def _pub_burden_vs_strength(ax, summary: pd.DataFrame):
    """Panel E — n_candidates vs top_contrib_total, colored by max phospho evidence."""
    cmap = mcolors.ListedColormap([PHOSPHO_COLORS[i] for i in [0, 1, 2, 3]])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    sc = ax.scatter(
        summary["top_contrib_total"], summary["n_candidates"],
        c=summary["max_n_sources"], cmap=cmap, norm=norm,
        s=24, alpha=0.82, edgecolors="none",
    )
    ax.axvline(THRESHOLD, color=GREY_500, lw=0.8, ls=(0, (3, 3)))

    # Label top 5 by evidence then candidate count
    top5 = summary.sort_values(
        ["max_n_sources", "n_candidates"], ascending=[False, False]
    ).head(5)
    for _, row in top5.iterrows():
        ax.annotate(
            _format_name(row["protein_name"]),
            xy=(row["top_contrib_total"], row["n_candidates"]),
            xytext=(5, 3), textcoords="offset points",
            fontsize=6, color=GREY_900,
        )

    cbar = plt.colorbar(sc, ax=ax, fraction=0.055, pad=0.02, aspect=22)
    cbar.set_label("# phospho\ndatabases", fontsize=6.5, labelpad=2)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.ax.tick_params(labelsize=6)

    ax.set_xlabel("Top contrib_total")
    ax.set_ylabel("# candidate residues")
    _panel_letter(ax, "E", "Candidate burden vs strength")
    _style_axis(ax, grid="both")


def _pub_length_vs_density(ax, summary: pd.DataFrame):
    """Panel F — sequence length vs candidate count; color = density (candidates/100 aa)."""
    df = summary[summary["seq_len"] > 0].copy()
    if df.empty:
        ax.text(
            0.5, 0.5,
            "Sequence length data unavailable.\n"
            "Pass --fasta for exact lengths;\n"
            "region-end proxy is used otherwise.",
            ha="center", va="center", transform=ax.transAxes,
            color=GREY_700, fontsize=8,
        )
        _panel_letter(ax, "F", "Candidate density vs protein length")
        _style_axis(ax)
        return

    df["density"] = df["n_candidates"] / (df["seq_len"] / 100.0)
    norm = mcolors.Normalize(vmin=df["density"].min(), vmax=df["density"].max())

    score_range = df["top_contrib_total"].max() - THRESHOLD
    sizes = 18 + 55 * (df["top_contrib_total"] - THRESHOLD) / (score_range if score_range > 0 else 1.0)

    sc = ax.scatter(
        df["seq_len"], df["n_candidates"],
        c=df["density"], cmap="YlOrRd", norm=norm,
        s=sizes, alpha=0.78, edgecolors="none",
    )

    # Median density reference line
    median_dens = df["density"].median()
    x_ref = np.linspace(df["seq_len"].min(), df["seq_len"].max(), 200)
    ax.plot(x_ref, median_dens * x_ref / 100.0, color=GREY_300, lw=0.9, ls="--",
            label=f"median density {median_dens:.2f}/100 aa")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")

    # Label top-5 highest-density proteins
    top5_dens = df.nlargest(5, "density")
    for _, row in top5_dens.iterrows():
        ax.annotate(
            f"{_format_name(row['protein_name'])} ({row['density']:.1f}/100aa)",
            xy=(row["seq_len"], row["n_candidates"]),
            xytext=(6, 3), textcoords="offset points",
            fontsize=6, color=GREY_900,
        )

    cbar = plt.colorbar(sc, ax=ax, fraction=0.014, pad=0.01, aspect=32)
    cbar.set_label("candidates / 100 aa", fontsize=6.5, labelpad=2)
    cbar.ax.tick_params(labelsize=6)

    ax.set_xlabel("Protein sequence length (aa)")
    ax.set_ylabel("# candidate residues")
    _panel_letter(ax, "F", "Candidate density vs protein length")
    _style_axis(ax, grid="both")


def _pub_phospho_evidence(ax, ptm_res: pd.DataFrame):
    """Panel G — horizontal bar chart of annotated residues by source combination."""
    if ptm_res.empty:
        ax.text(0.5, 0.5, "No PTM data\n(run annotate_ptm.py first)",
                ha="center", va="center", transform=ax.transAxes, color=GREY_700, fontsize=8)
        _panel_letter(ax, "G", "Phospho evidence by source combination")
        _style_axis(ax)
        return

    source_counts = ptm_res["sources"].value_counts().sort_values(ascending=True)

    def _n_src(s: str) -> int:
        return len(str(s).split(";"))

    bar_colors = [PHOSPHO_COLORS.get(_n_src(s), GREY_500) for s in source_counts.index]
    y = np.arange(len(source_counts))
    ax.barh(y, source_counts.values, color=bar_colors, alpha=0.85, height=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([str(s) for s in source_counts.index], fontsize=7)
    ax.set_xlabel("# candidate residues")

    for yi, count in zip(y, source_counts.values):
        ax.text(count + 0.3, yi, str(int(count)), va="center", fontsize=7, color=GREY_700)

    legend_handles = [
        Patch(facecolor=PHOSPHO_COLORS[3], alpha=0.85, label="3 databases"),
        Patch(facecolor=PHOSPHO_COLORS[2], alpha=0.85, label="2 databases"),
        Patch(facecolor=PHOSPHO_COLORS[1], alpha=0.85, label="1 database"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=6.5, loc="lower right")
    _panel_letter(ax, "G", "Phospho evidence by source combination")
    _style_axis(ax, grid="x")


def _pub_top_annotated(ax, ptm_res: pd.DataFrame, n: int = 20):
    """Panel H — top N residues from annotated_summary sorted by n_sources then contrib_total."""
    if ptm_res.empty:
        ax.text(0.5, 0.5, "No PTM data\n(run annotate_ptm.py first)",
                ha="center", va="center", transform=ax.transAxes, color=GREY_700, fontsize=8)
        _panel_letter(ax, "H", f"Top {n} phospho-prioritized candidates")
        _style_axis(ax)
        return

    df = (
        ptm_res
        .sort_values(["n_sources", "contrib_total"], ascending=[False, False])
        .head(n)
        .sort_values("contrib_total", ascending=True)
        .reset_index(drop=True)
    )

    y = np.arange(len(df))

    if "protein_name" in df.columns:
        y_labels = [
            f"{_format_name(str(r.protein_name))} {r.aa}{int(r.position)}"
            for r in df.itertuples()
        ]
    else:
        y_labels = [
            f"{r.uniprot_id} {r.aa}{int(r.position)}"
            for r in df.itertuples()
        ]

    dot_colors = [PHOSPHO_COLORS.get(int(ns), GREY_500) for ns in df["n_sources"]]
    n_pubmed_vals = pd.to_numeric(df["n_pubmed"], errors="coerce").fillna(1)
    dot_sizes = (22 + 8 * n_pubmed_vals.clip(upper=10)).values

    ax.hlines(y, THRESHOLD, df["contrib_total"], color=GREY_300, lw=1.0, zorder=1)
    ax.scatter(df["contrib_total"], y, s=dot_sizes, color=dot_colors,
               edgecolor="white", linewidth=0.5, zorder=3)
    ax.axvline(THRESHOLD, color=GREY_500, lw=0.8, ls=(0, (3, 3)))

    ax.set_yticks(y)
    ax.set_yticklabels(y_labels, fontsize=6.5)
    ax.set_xlabel("contrib_total")

    legend_handles = [
        Patch(facecolor=PHOSPHO_COLORS[3], label="3 DBs"),
        Patch(facecolor=PHOSPHO_COLORS[2], label="2 DBs"),
        Patch(facecolor=PHOSPHO_COLORS[1], label="1 DB"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=6.5, loc="lower right")
    _panel_letter(ax, "H", f"Top {len(df)} phospho-prioritized candidates")
    _style_axis(ax, grid="x")


# ---------------------------------------------------------------------------
# Publication figure — 8-panel layout
# ---------------------------------------------------------------------------

def plot_publication_summary(
    cohort_dir: Path,
    master: pd.DataFrame,
    summary: pd.DataFrame,
    funnel: Dict[str, Optional[int]],
    ptm_res: pd.DataFrame,
    ptm_prot: pd.DataFrame,
    out_path: Path,
):
    has_ptm = not ptm_res.empty
    n_cands_str = _count_text(funnel["with_candidates"])
    n_res_str = _count_text(funnel["total_residues"])
    subtitle = f"{n_cands_str} candidate proteins  ·  {n_res_str} candidate residues"
    if has_ptm:
        n_ph = int(ptm_res["uniprot_id"].nunique())
        subtitle += f"  ·  {n_ph} with phospho evidence (step 6)"

    fig = plt.figure(figsize=(10.0, 15.5))
    fig.patch.set_facecolor("white")
    fig.text(0.055, 0.984, "Brain-enriched protein biomarker screen",
             ha="left", va="top", fontsize=13, fontweight="bold", color=BLACK)
    fig.text(0.055, 0.966, subtitle,
             ha="left", va="top", fontsize=8, color=GREY_700)

    # 5-row × 2-col grid
    # Rows 2 and 4 (C and F) span both columns
    gs = gridspec.GridSpec(
        5, 2, figure=fig,
        left=0.10, right=0.97, top=0.940, bottom=0.032,
        wspace=0.48, hspace=0.62,
        height_ratios=[1.10, 0.82, 1.18, 0.82, 2.25],
    )

    ax_a = fig.add_subplot(gs[0, 0])   # A — funnel
    ax_b = fig.add_subplot(gs[0, 1])   # B — score distribution
    ax_c = fig.add_subplot(gs[1, :])   # C — ranked strength (full width)
    ax_d = fig.add_subplot(gs[2, 0])   # D — CathL vs CathB
    ax_e = fig.add_subplot(gs[2, 1])   # E — burden vs strength
    ax_f = fig.add_subplot(gs[3, :])   # F — length vs density (full width)
    ax_g = fig.add_subplot(gs[4, 0])   # G — phospho evidence breakdown
    ax_h = fig.add_subplot(gs[4, 1])   # H — top annotated candidates

    _pub_funnel(ax_a, funnel, ptm_res)
    _pub_score_distribution(ax_b, master, ptm_res)
    _pub_ranked_strength(ax_c, summary)
    _pub_contribution_balance(ax_d, master)
    _pub_burden_vs_strength(ax_e, summary)
    _pub_length_vs_density(ax_f, summary)
    _pub_phospho_evidence(ax_g, ptm_res)
    _pub_top_annotated(ax_h, ptm_res)

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Exploration detail figure (unchanged structure)
# ---------------------------------------------------------------------------

def _explore_ranked_all(ax, summary: pd.DataFrame):
    df = summary.sort_values("top_contrib_total", ascending=False).reset_index(drop=True)
    x = np.arange(1, len(df) + 1)
    norm = mcolors.Normalize(vmin=df["n_candidates"].min(), vmax=df["n_candidates"].max())
    cmap = _candidate_cmap()
    ax.scatter(x, df["top_contrib_total"], c=df["n_candidates"], cmap=cmap, norm=norm,
               s=20, alpha=0.85, edgecolors="none")
    ax.plot(x, df["top_contrib_total"], color=GREY_300, lw=0.8, zorder=0)
    ax.axhline(THRESHOLD, color=GREY_500, lw=0.8, ls=(0, (3, 3)))
    ax.set_xlabel("Protein rank")
    ax.set_ylabel("Top contrib_total")
    ax.set_title("All candidate proteins by strongest residue",
                 loc="left", fontsize=10, fontweight="bold")
    _style_axis(ax, grid="y")
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.012)
    cbar.set_label("# candidate residues", fontsize=7)
    cbar.ax.tick_params(labelsize=6)


def _explore_burden(ax, summary: pd.DataFrame, top_n: int):
    df = summary.nlargest(top_n, "n_candidates").sort_values("n_candidates", ascending=True)
    y = np.arange(len(df))
    ax.barh(y, df["n_candidates"], color=BLUE, alpha=0.78)
    ax.set_yticks(y)
    ax.set_yticklabels([_format_name(v) for v in df["protein_name"]], fontsize=6.5)
    ax.set_xlabel("# candidate residues")
    ax.set_title(f"Top {top_n} proteins by candidate burden",
                 loc="left", fontsize=10, fontweight="bold")
    _style_axis(ax, grid="x")
    for yi, row in zip(y, df.itertuples()):
        ax.text(row.n_candidates + 0.35, yi, f"{int(row.n_candidates)}",
                va="center", fontsize=6.5, color=GREY_700)


def _explore_regions_vs_candidates(ax, summary: pd.DataFrame):
    ax.scatter(summary["n_regions"], summary["n_candidates"], s=22,
               c=summary["top_contrib_total"], cmap="YlOrRd", alpha=0.75, edgecolors="none")
    ax.set_xlabel("# candidate regions")
    ax.set_ylabel("# candidate residues")
    ax.set_title("Regional spread vs candidate burden", loc="left", fontsize=10, fontweight="bold")
    _style_axis(ax, grid="both")
    sm = cm.ScalarMappable(
        cmap=cm.get_cmap("YlOrRd"),
        norm=mcolors.Normalize(
            vmin=summary["top_contrib_total"].min(),
            vmax=summary["top_contrib_total"].max(),
        ),
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.035, pad=0.015)
    cbar.set_label("top contrib_total", fontsize=7)
    cbar.ax.tick_params(labelsize=6)


def _explore_all_labelled(ax, summary: pd.DataFrame):
    df = summary.sort_values("top_contrib_total", ascending=True).reset_index(drop=True)
    y = np.arange(len(df))
    norm = mcolors.Normalize(vmin=df["n_candidates"].min(), vmax=df["n_candidates"].max())
    cmap = _candidate_cmap()
    colors = cmap(norm(df["n_candidates"]))
    ax.barh(y, df["top_contrib_total"], color=colors, height=0.72)
    ax.axvline(THRESHOLD, color=GREY_500, lw=0.8, ls=(0, (3, 3)))
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{_format_name(r.protein_name)} ({r.uniprot_id})" for r in df.itertuples()],
        fontsize=4.9,
    )
    ax.set_xlabel("Top contrib_total")
    ax.set_title("Every candidate protein", loc="left", fontsize=10, fontweight="bold")
    _style_axis(ax, grid="x")


def plot_exploration_detail(
    master: pd.DataFrame,
    summary: pd.DataFrame,
    top_n: int,
    out_path: Path,
):
    height = max(14.0, 0.13 * len(summary))
    fig = plt.figure(figsize=(11.5, height))
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(
        4, 2, figure=fig,
        left=0.08, right=0.98, top=0.96, bottom=0.035,
        hspace=0.55, wspace=0.35,
        height_ratios=[1.0, 1.0, 1.0, max(2.6, 0.017 * len(summary))],
    )
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2, :])
    ax5 = fig.add_subplot(gs[3, :])

    _explore_ranked_all(ax1, summary)
    _explore_burden(ax2, summary, top_n)
    _explore_regions_vs_candidates(ax3, summary)
    _pub_contribution_balance(ax4, master)
    _explore_all_labelled(ax5, summary)

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fast cohort summary and visualization from merged CSV files only."
    )
    parser.add_argument(
        "--cohort-dir", required=True,
        help="Path to cohort run directory containing funnel.csv, summary.csv, master CSV.",
    )
    parser.add_argument(
        "--fasta", default=None,
        help="Input FASTA for exact protein sequence lengths (optional; region-end proxy used otherwise).",
    )
    parser.add_argument(
        "--ptm", default=None,
        help=(
            "Path to step6_ptm/annotated_summary.csv "
            "(default: auto-detected from <cohort-dir>/step6_ptm/)."
        ),
    )
    parser.add_argument(
        "--top-n", type=int, default=30,
        help="Proteins shown in exploration burden panels (default: 30).",
    )
    parser.add_argument(
        "--out", default=None,
        help="Publication PNG path (default: <cohort-dir>/cohort_publication_summary.png).",
    )
    parser.add_argument(
        "--explore-out", "--all-out", dest="explore_out", default=None,
        help="Exploration PNG path (default: <cohort-dir>/cohort_exploration_detail.png).",
    )
    parser.add_argument(
        "--report-out", default=None,
        help="Text report path (default: <cohort-dir>/cohort_summary.txt).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    cohort_dir = Path(args.cohort_dir)
    if not cohort_dir.is_dir():
        sys.exit(f"ERROR: cohort directory not found: {cohort_dir}")

    fasta_path = Path(args.fasta) if args.fasta else None

    # Auto-detect PTM file; warn only if user explicitly passed --ptm
    if args.ptm:
        ptm_path: Optional[Path] = Path(args.ptm)
        warn_ptm = True
    else:
        auto_ptm = cohort_dir / "step6_ptm" / "annotated_summary.csv"
        ptm_path = auto_ptm if auto_ptm.exists() else None
        warn_ptm = False

    pub_path = Path(args.out) if args.out else cohort_dir / "cohort_publication_summary.png"
    explore_path = (
        Path(args.explore_out) if args.explore_out else cohort_dir / "cohort_exploration_detail.png"
    )
    report_path = (
        Path(args.report_out) if args.report_out else cohort_dir / "cohort_summary.txt"
    )

    master, summary, funnel, warnings, ptm_res, ptm_prot = load_cohort_data(
        cohort_dir, fasta_path=fasta_path, ptm_path=ptm_path, warn_ptm_missing=warn_ptm
    )

    report_lines = build_report_lines(cohort_dir, master, summary, funnel, warnings, ptm_res)
    save_and_print_report(report_lines, report_path)
    plot_publication_summary(cohort_dir, master, summary, funnel, ptm_res, ptm_prot, pub_path)
    plot_exploration_detail(master, summary, args.top_n, explore_path)


if __name__ == "__main__":
    main()
