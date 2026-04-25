#!/usr/bin/env python3
"""
PTM / Phosphorylation Site Annotation  (Step 6)
=================================================
Annotates protease-resistant biomarker candidates from a completed
ProsperousPlus pipeline run with known human phosphorylation evidence from
three open databases:

  - EPSD 2.0  (Eukaryotic Phosphorylation Site Database)
  - dbPTM     (experimental phosphorylation sites, all species → filtered to human)
  - UniProt   (Swiss-Prot curated Modified residue features)

Starting point
--------------
The primary input is the cohort-level ``master_biomarker_candidates.csv``,
which aggregates every individual S/T residue (one row per residue) that
crossed the protease-resistance threshold (contrib_total ≥ 0.475) across
all proteins in the cohort run.

For each of those 916 (protein, position) pairs the script asks:
"Is this specific amino acid position known to be phosphorylated in humans?"

Results are ranked by n_sources (number of databases confirming the site)
then by contrib_total (strength of protease resistance), so the top rows are
residues that are BOTH strongly protease-resistant AND well-supported as
human phosphorylation sites.

Usage
-----
    conda activate prosperousplus

    # Recommended: point to summary.csv — script auto-resolves to
    # master_biomarker_candidates.csv in the same cohort folder
    python annotate_ptm.py \\
        --candidates results/cohort_<name>/summary.csv

    # Or point directly to master_biomarker_candidates.csv
    python annotate_ptm.py \\
        --candidates results/cohort_<name>/master_biomarker_candidates.csv

    # Single-protein run (per-protein biomarker_candidates.csv also works)
    python annotate_ptm.py \\
        --candidates results/<run>/step5_scores/biomarker_candidates.csv

Default source file paths (relative to the working directory) match the layout
documented in CLAUDE.md and are used when --epsd / --dbptm / --uniprot are
omitted.  Any source whose file is missing is skipped with a warning.
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from lib.ptm_annotate import (
    load_epsd,
    load_dbptm,
    load_uniprot_modres,
    annotate_candidates,
    summarize_annotations,
)

# ── Default data paths ────────────────────────────────────────────────────────
_DEFAULT_EPSD   = Path("data/reference/ptm_sources/epsd/Homo sapiens.txt")
_DEFAULT_DBPTM  = Path("data/reference/ptm_sources/dbptm/Phosphorylation")
_DEFAULT_UNIPROT = Path("data/reference/ptm_sources/uniprot/human_reviewed_modres.tsv")


def _elapsed(t0: float) -> str:
    s = time.time() - t0
    return f"{s:.1f}s" if s < 60 else f"{s / 60:.1f}min"


def main():
    parser = argparse.ArgumentParser(
        description="Annotate ProsperousPlus biomarker candidates with phospho-site databases"
    )
    parser.add_argument(
        "--candidates", required=True,
        help=(
            "Path to biomarker_candidates.csv or per_residue_contributions.csv "
            "from a completed pipeline run (step5_scores/)"
        ),
    )
    parser.add_argument(
        "--epsd", default=None,
        help=(
            f"Path to EPSD 'Homo sapiens.txt' "
            f"(default: {_DEFAULT_EPSD})"
        ),
    )
    parser.add_argument(
        "--dbptm", default=None,
        help=(
            f"Path to dbPTM 'Phosphorylation' file "
            f"(default: {_DEFAULT_DBPTM})"
        ),
    )
    parser.add_argument(
        "--uniprot", default=None,
        help=(
            f"Path to UniProt human_reviewed_modres.tsv "
            f"(default: {_DEFAULT_UNIPROT})"
        ),
    )
    parser.add_argument(
        "--out", default=None,
        help=(
            "Output directory for step6_ptm/ results "
            "(default: <run_dir>/step6_ptm/ inferred from --candidates path)"
        ),
    )
    args = parser.parse_args()

    candidates_path = Path(args.candidates)
    if not candidates_path.exists():
        sys.exit(f"ERROR: candidates file not found: {candidates_path}")

    epsd_path    = Path(args.epsd)    if args.epsd    else _DEFAULT_EPSD
    dbptm_path   = Path(args.dbptm)   if args.dbptm   else _DEFAULT_DBPTM
    uniprot_path = Path(args.uniprot) if args.uniprot else _DEFAULT_UNIPROT

    # Resolve the actual candidates file and output directory.
    #
    # Three accepted inputs:
    #   (a) cohort summary.csv          → use master_biomarker_candidates.csv
    #                                      from the same directory
    #   (b) master_biomarker_candidates.csv (cohort root)
    #                                    → use directly; output → sibling step6_ptm/
    #   (c) step5_scores/biomarker_candidates.csv (single-protein run)
    #                                    → use directly; output → <run>/step6_ptm/

    if candidates_path.name == "summary.csv":
        master_path = candidates_path.parent / "master_biomarker_candidates.csv"
        if not master_path.exists():
            sys.exit(
                f"ERROR: Expected master_biomarker_candidates.csv next to summary.csv "
                f"but it was not found:\n  {master_path}"
            )
        print(f"  [info] summary.csv detected — using {master_path.name} for annotation")
        candidates_path = master_path

    # Output dir: sibling of wherever candidates_path lives
    if args.out:
        out_dir = Path(args.out)
    elif candidates_path.name == "master_biomarker_candidates.csv":
        # Cohort root — put step6_ptm/ next to master file
        out_dir = candidates_path.parent / "step6_ptm"
    else:
        # Per-protein run: candidates_path is <run>/step5_scores/<file>.csv
        out_dir = candidates_path.parent.parent / "step6_ptm"

    out_dir.mkdir(parents=True, exist_ok=True)

    wall_start = time.time()

    print(f"\n{'='*60}")
    print(f" Candidates: {candidates_path}")
    print(f" Output:     {out_dir}")
    print(f"{'='*60}")

    # ── Load candidates ───────────────────────────────────────────────────────
    candidates = pd.read_csv(candidates_path)
    print(f"\n  Loaded {len(candidates):,} candidate residue(s) from {candidates_path.name}")
    if candidates.empty:
        print("  No candidates to annotate — done.")
        return

    # ── Load UniProt (always first — provides human_accessions for dbPTM) ─────
    uniprot_sites = None
    human_accessions = set()

    if uniprot_path.exists():
        t0 = time.time()
        print(f"\n[Source 1] UniProt  ({uniprot_path})")
        uniprot_sites, human_accessions = load_uniprot_modres(uniprot_path)
        print(f"  {len(human_accessions):,} human accessions indexed")
        print(f"  {len(uniprot_sites):,} phospho-site rows parsed")
        print(f"  done in {_elapsed(t0)}")
    else:
        print(f"\n[Source 1] UniProt  — file not found, skipping ({uniprot_path})")

    # ── Load EPSD ─────────────────────────────────────────────────────────────
    epsd_sites = None
    if epsd_path.exists():
        t0 = time.time()
        print(f"\n[Source 2] EPSD     ({epsd_path})")
        epsd_sites = load_epsd(epsd_path)
        print(f"  {len(epsd_sites):,} phosphorylation site rows")
        print(f"  done in {_elapsed(t0)}")
    else:
        print(f"\n[Source 2] EPSD     — file not found, skipping ({epsd_path})")

    # ── Load dbPTM (filtered to human using UniProt accessions) ──────────────
    dbptm_sites = None
    if dbptm_path.exists():
        if not human_accessions:
            print(
                f"\n[Source 3] dbPTM    — WARNING: UniProt accessions not available; "
                f"dbPTM will not be filtered to human proteins. Load UniProt first."
            )
        t0 = time.time()
        print(f"\n[Source 3] dbPTM    ({dbptm_path})")
        dbptm_sites = load_dbptm(dbptm_path, human_accessions)
        print(f"  {len(dbptm_sites):,} human phosphorylation site rows (after species filter)")
        print(f"  done in {_elapsed(t0)}")
    else:
        print(f"\n[Source 3] dbPTM    — file not found, skipping ({dbptm_path})")

    # ── Annotate ──────────────────────────────────────────────────────────────
    t0 = time.time()
    print(f"\n[Annotate] Joining candidates with PTM databases...")
    hits = annotate_candidates(
        candidates,
        epsd=epsd_sites,
        dbptm=dbptm_sites,
        uniprot=uniprot_sites,
    )
    print(f"  {len(hits):,} candidate × source hit(s)")
    print(f"  done in {_elapsed(t0)}")

    # ── Summarize ─────────────────────────────────────────────────────────────
    summary = summarize_annotations(hits)
    n_annotated = len(summary)
    n_multi     = (summary["n_sources"] > 1).sum() if not summary.empty else 0

    # ── Save ──────────────────────────────────────────────────────────────────
    hits_path    = out_dir / "annotated_candidates.csv"
    summary_path = out_dir / "annotated_summary.csv"

    if not hits.empty:
        hits.to_csv(hits_path, index=False)
    if not summary.empty:
        summary.to_csv(summary_path, index=False)

    # ── Report ────────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  Candidates annotated : {n_annotated} / {len(candidates)}")
    print(f"  Supported by ≥2 sources : {n_multi}")

    if not summary.empty:
        print(f"\n  Top annotated residues:")
        hdr = (f"  {'#':>3}  {'UniProt':>10}  {'Protein':<14}  "
               f"{'Pos':>5}  {'AA':>2}  {'Contrib':>8}  {'Src':>3}  {'Sources':<20}  PMIDs")
        print(hdr)
        print(f"  {'-'*80}")
        for i, row in summary.head(15).iterrows():
            print(
                f"  {i+1:>3}  {row.get('uniprot_id',''):>10}  "
                f"{row.get('protein_name',''):<14}  "
                f"{row.get('position', ''):>5}  {row.get('aa', ''):>2}  "
                f"{row.get('contrib_total', 0):>8.4f}  "
                f"{row.get('n_sources', 0):>3}  "
                f"{row.get('sources',''):<20}  {row.get('n_pubmed', 0)}"
            )

    print(f"\n  Saved: {hits_path}  ({len(hits):,} rows)")
    print(f"  Saved: {summary_path}  ({n_annotated} rows)")
    print(f"\n  Total time: {_elapsed(wall_start)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
