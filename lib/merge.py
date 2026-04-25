"""Merge per-chunk biomarker candidate results from a cohort run."""

import re
from pathlib import Path

import pandas as pd


# sequence_id format: sp|P10636|TAU_HUMAN_402-461
_SEQ_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|(\S+?)_(\d+-\d+)$")


def _parse_sequence_id(seq_id: str) -> tuple:
    """Return (uniprot_id, protein_name, region) from a sequence_id string."""
    m = _SEQ_RE.match(seq_id)
    if m:
        return m.group(1), m.group(2), m.group(3)
    # Fallback: return the whole thing
    return seq_id, seq_id, ""


def merge_cohort_results(cohort_dir: Path) -> tuple:
    """Collect all per-chunk biomarker_candidates CSVs and produce master outputs.

    Globs ``chunk_*/step5_scores/biomarker_candidates.csv`` under *cohort_dir*,
    enriches each row with parsed protein metadata, concatenates, deduplicates,
    and writes:

    - ``<cohort_dir>/master_biomarker_candidates.csv``
    - ``<cohort_dir>/summary.csv``

    Parameters
    ----------
    cohort_dir : Path
        Root output directory for the cohort run.

    Returns
    -------
    master : pd.DataFrame
        All candidate rows, sorted by ``mutation_sum`` descending.
    summary : pd.DataFrame
        One row per protein with candidate counts and top hit.
    """
    cohort_dir = Path(cohort_dir)
    candidate_files = sorted(cohort_dir.glob("chunk_*/step5_scores/biomarker_candidates.csv"))

    if not candidate_files:
        print("  [merge] No biomarker_candidates.csv files found — nothing to merge.")
        return pd.DataFrame(), pd.DataFrame()

    frames = []
    for csv_path in candidate_files:
        chunk_name = csv_path.parts[-3]  # e.g. chunk_000
        df = pd.read_csv(csv_path)
        if df.empty:
            continue
        df["chunk"] = chunk_name
        frames.append(df)

    if not frames:
        print("  [merge] All chunk candidate files are empty.")
        return pd.DataFrame(), pd.DataFrame()

    master = pd.concat(frames, ignore_index=True)

    # Parse sequence_id into protein metadata columns
    parsed = master["sequence_id"].apply(lambda s: pd.Series(
        _parse_sequence_id(s), index=["uniprot_id", "protein_name", "region"]
    ))
    master = pd.concat([master, parsed], axis=1)

    # Re-order columns for readability
    col_order = [
        "uniprot_id", "protein_name", "region",
        "sequence_id", "position", "peptide",
        "cathL_original", "cathB_original",
        "cathL_diff", "cathB_diff", "mutation_sum",
        "chunk",
    ]
    master = master[[c for c in col_order if c in master.columns]]
    master = master.sort_values("mutation_sum", ascending=False).reset_index(drop=True)

    # ── Summary: one row per protein ─────────────────────────────────────────
    summary_rows = []
    for (uid, pname), grp in master.groupby(["uniprot_id", "protein_name"], sort=False):
        top = grp.iloc[0]
        summary_rows.append({
            "uniprot_id":       uid,
            "protein_name":     pname,
            "n_candidates":     len(grp),
            "n_regions":        grp["region"].nunique(),
            "top_peptide":      top["peptide"],
            "top_position":     top["position"],
            "top_region":       top["region"],
            "top_mutation_sum": round(top["mutation_sum"], 4),
        })
    summary = (
        pd.DataFrame(summary_rows)
        .sort_values("top_mutation_sum", ascending=False)
        .reset_index(drop=True)
    )

    # ── Save ─────────────────────────────────────────────────────────────────
    master_path  = cohort_dir / "master_biomarker_candidates.csv"
    summary_path = cohort_dir / "summary.csv"
    master.to_csv(master_path,  index=False)
    summary.to_csv(summary_path, index=False)

    # ── Print summary to stdout ───────────────────────────────────────────────
    n_proteins   = len(summary)
    n_candidates = len(master)
    n_with_hits  = (summary["n_candidates"] > 0).sum()

    print(f"\n  Merged {len(candidate_files)} chunk(s):")
    print(f"    Proteins with candidates : {n_with_hits} / {n_proteins}")
    print(f"    Total candidate rows     : {n_candidates:,}")
    print(f"\n  Top 10 candidates across all proteins:")
    print(f"  {'#':>3}  {'UniProt':>10}  {'Protein':<16}  {'Peptide':<10}  "
          f"{'Pos':>6}  {'Mut Sum':>8}")
    print(f"  {'-'*65}")
    for i, row in master.head(10).iterrows():
        print(f"  {i+1:>3}  {row['uniprot_id']:>10}  {row['protein_name']:<16}  "
              f"{row['peptide']:<10}  {row['position']:>6}  {row['mutation_sum']:>8.4f}")

    print(f"\n  Saved: {master_path}")
    print(f"  Saved: {summary_path}")

    return master, summary
