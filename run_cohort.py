#!/usr/bin/env python3
"""
Cohort Biomarker Pipeline
=========================
Two-phase architecture that minimises ML model loads:

  Phase 1 — Cohort Step 1
      Run CathD + Legumain on the full cohort FASTA in batches of
      STEP1_BATCH_SIZE proteins, with both proteases running in parallel
      (each on half the available cores).  Results are split per-protein
      and written to each protein's step1_predictions/ directory.

  Phase 2 — Per-protein Steps 2-5
      Each protein's Steps 2-5 (region detection → mutant generation →
      mutant predictions → scoring) run independently in parallel
      (--parallel proteins simultaneously).  Step 1 is skipped because
      the CSVs already exist from Phase 1.

No intermediate per-protein FASTA files are saved to disk — each worker
writes a small temp file and deletes it once the subprocess exits.

Usage
-----
    conda activate prosperousplus
    cd /shared/macdata/groups/ppc/projects/ProsperousPlus

    # Test: first 5 proteins, 5 parallel workers for steps 2-5
    python run_cohort.py \\
        --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \\
        --parallel 5 --cores 20 --limit 5

    # Full run
    python run_cohort.py \\
        --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \\
        --parallel 5 --cores 20

    # Resume (both phases detect already-complete work and skip it)
    python run_cohort.py --fasta ... --parallel 5 --cores 20

    # Merge only (regenerate master files from completed proteins)
    python run_cohort.py --fasta ... --parallel 5 --cores 20 --merge-only

    # Tail a protein's log for detail
    tail -f results/cohort_.../sp|P10636|TAU_HUMAN/pipeline.log

Output
------
    results/cohort_<fasta_stem>_YYYY-MM-DD_HHMM/
        _step1_batches/          ← shared FASTA chunks (50 proteins each)
        _step1_cathD/            ← raw CathD batch results (cached for resume)
        _step1_legumain/         ← raw Legumain batch results
        master_biomarker_candidates.csv
        summary.csv
        sp|P19622|HME2_HUMAN/
            pipeline.log
            step1_predictions/
            step2_regions/
            step3_fastas/
            step4_mutant_predictions/
            step5_scores/
                biomarker_candidates.csv
        sp|O15079|SNPH_HUMAN/
            ...
"""

import argparse
import datetime
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))

from lib import run_prediction_batched

# Proteins per ProsperousPlus call in Phase 1 (Step 1 cohort run)
STEP1_BATCH_SIZE = 50

# Protease model IDs for Step 1
_CATHEPSIN_D = "A01.009"
_LEGUMAIN    = "C13.004"


# ── Helpers ──────────────────────────────────────────────────────────────────

def _elapsed(t0: float) -> str:
    s = time.time() - t0
    if s < 60:
        return f"{s:.0f}s"
    elif s < 3600:
        return f"{s/60:.1f}min"
    else:
        return f"{s/3600:.1f}hr"


def _eta(done: int, total: int, t0: float) -> str:
    if done == 0:
        return "?"
    elapsed = time.time() - t0
    remaining = (total - done) / (done / elapsed)
    if remaining < 60:
        return f"{remaining:.0f}s"
    elif remaining < 3600:
        return f"{remaining/60:.0f}min"
    else:
        return f"{remaining/3600:.1f}hr"


# Sentinel files that mark each pipeline step complete
_STEP_SENTINELS = [
    ("step1", "step1_predictions/cathepsinD_cleavage.csv"),
    ("step2", "step2_regions/low_cleavage_regions.csv"),
    ("step3", "step3_fastas/original.fasta"),
    ("step4", "step4_mutant_predictions/cathepsinB_single_E.csv"),
    ("step5", "step5_scores/biomarker_candidates.csv"),
]


def _chunk_status(chunk_dir: Path) -> str:
    """Return a short status string for one protein's output directory."""
    if not chunk_dir.exists():
        return "queued"
    # pipeline.log is written by the worker the moment it starts — without it
    # the protein has step1 CSVs from Phase 1 but hasn't been picked up yet.
    if not (chunk_dir / "pipeline.log").exists():
        return "queued"
    last_done: Optional[str] = None
    for step_name, sentinel in _STEP_SENTINELS:
        if (chunk_dir / sentinel).exists():
            last_done = step_name
        else:
            if last_done is None:
                return "step1 running..."
            return f"{last_done} done  |  {step_name} running..."
    return "COMPLETE"


def _progress_monitor(
    protein_dirs: List[tuple],    # [(short_name, out_path), ...]
    t0: float,
    stop_event: threading.Event,
    total_proteins: int,
    completed_counter: list,      # mutable [int] bumped by main thread
    max_display: int = 5,         # cap lines shown to the parallel worker count
):
    """Background thread: print per-protein status every 30 seconds."""
    while not stop_event.wait(30):
        done    = completed_counter[0]
        elapsed = _elapsed(t0)
        eta     = _eta(done, total_proteins, t0)
        print(f"\n[Progress  {elapsed}  |  {done}/{total_proteins} done  |  ETA {eta}]")
        running = [
            (name, _chunk_status(out_path))
            for name, out_path in protein_dirs
            if _chunk_status(out_path) not in ("COMPLETE", "queued")
        ]
        for name, status in running[:max_display]:
            print(f"  {name:<35}  {status}")
        if len(running) > max_display:
            print(f"  (+ {len(running) - max_display} more in queue)")
        print(flush=True)


# ── Phase 1: cohort-level Step 1 ─────────────────────────────────────────────

def _run_step1_cohort(
    fasta_file: Path,
    cohort_dir: Path,
    all_proteins: List[tuple],   # [(protein_id, sequence), ...]
    batch_size: int = STEP1_BATCH_SIZE,
    total_cores: int = 20,
) -> None:
    """Run CathD + Legumain on the full cohort FASTA (batched), then split per protein.

    Writes per-protein CSVs:
        cohort_dir/<safe_id>/step1_predictions/cathepsinD_cleavage.csv
        cohort_dir/<safe_id>/step1_predictions/legumain_cleavage.csv

    Proteins that already have both CSVs are skipped at the split step.
    Batch-level resume is handled by run_prediction_batched internally.
    """
    # Determine which proteins still need their step1 CSVs written
    need_split = []
    for protein_id, _ in all_proteins:
        safe_id  = protein_id.replace("/", "_").replace("\\", "_")
        pred_dir = cohort_dir / safe_id / "step1_predictions"
        if not (
            (pred_dir / "cathepsinD_cleavage.csv").exists() and
            (pred_dir / "cathepsinD_cleavage.csv").stat().st_size > 0 and
            (pred_dir / "legumain_cleavage.csv").exists() and
            (pred_dir / "legumain_cleavage.csv").stat().st_size > 0
        ):
            need_split.append(protein_id)

    if not need_split:
        print(f"  [step1] All {len(all_proteins)} proteins already have step1 CSVs — skipping")
        return

    n_batches = -(-len(all_proteins) // batch_size)   # ceiling division
    pn = max(1, total_cores // 2)                      # cores per protease
    print(f"  {len(need_split)}/{len(all_proteins)} protein(s) need step1 CSVs")
    print(f"  {n_batches} batch(es) of ≤{batch_size} proteins  |  "
          f"CathD + Legumain in parallel  |  {pn} cores each")

    # Intermediate dirs (kept for resume — batch results are cached here)
    batch_dir    = cohort_dir / "_step1_batches"
    cathD_dir    = cohort_dir / "_step1_cathD"
    legumain_dir = cohort_dir / "_step1_legumain"

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=2) as ex:
        fut_D = ex.submit(
            run_prediction_batched,
            fasta_file, _CATHEPSIN_D, cathD_dir,
            batch_size=batch_size, batch_dir=batch_dir,
            process_num=pn, verbose=True,
        )
        fut_L = ex.submit(
            run_prediction_batched,
            fasta_file, _LEGUMAIN, legumain_dir,
            batch_size=batch_size, batch_dir=batch_dir,
            process_num=pn, verbose=True,
        )
        df_D = fut_D.result()
        df_L = fut_L.result()

    elapsed = time.time() - t0
    print(f"  Predictions done in {elapsed/60:.1f} min  "
          f"({len(df_D):,} CathD rows, {len(df_L):,} Legumain rows)")
    print(f"  Splitting results per protein…")

    # Write per-protein CSVs (skip proteins that already have them)
    written = 0
    missing = []
    for protein_id, _ in all_proteins:
        safe_id  = protein_id.replace("/", "_").replace("\\", "_")
        pred_dir = cohort_dir / safe_id / "step1_predictions"
        pred_dir.mkdir(parents=True, exist_ok=True)

        cathD_csv    = pred_dir / "cathepsinD_cleavage.csv"
        legumain_csv = pred_dir / "legumain_cleavage.csv"

        prot_D = df_D[df_D["sequence_id"] == protein_id]
        prot_L = df_L[df_L["sequence_id"] == protein_id]

        if prot_D.empty and prot_L.empty:
            missing.append(protein_id)
            continue

        if not (cathD_csv.exists() and cathD_csv.stat().st_size > 0):
            prot_D.to_csv(cathD_csv, index=False)
            written += 1

        if not (legumain_csv.exists() and legumain_csv.stat().st_size > 0):
            prot_L.to_csv(legumain_csv, index=False)

    print(f"  {written} protein(s) written, "
          f"{len(all_proteins) - written - len(missing)} already existed")
    if missing:
        print(f"  WARNING: {len(missing)} protein(s) had no predictions "
              f"(IDs may not match FASTA headers): {missing[:5]}")


# ── Worker ────────────────────────────────────────────────────────────────────

def _run_protein(args: tuple) -> tuple:
    """Worker: write a temp FASTA, run run_pipeline.py, delete the temp file.

    Returns (protein_id, success, elapsed_seconds).
    """
    protein_id, sequence, output_path, cores_per_job, pipeline_script = args
    t0 = time.time()

    # run_pipeline.py prepends "results/" to --output; strip that prefix
    results_root = Path("results")
    try:
        rel_output = output_path.relative_to(results_root)
    except ValueError:
        rel_output = output_path

    # Write a temporary single-protein FASTA — deleted after subprocess exits
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".fasta")
    try:
        with os.fdopen(tmp_fd, "w") as f:
            f.write(f">{protein_id}\n{sequence}\n")

        cmd = [
            sys.executable, str(pipeline_script),
            "--fasta",       tmp_path,
            "--output",      str(rel_output),
            "--cores",       str(cores_per_job),
            "--batch-size",  "1",          # single protein → single batch
            "--skip-step1",                # step1 CSVs written by Phase 1
        ]

        log_path = output_path / "pipeline.log"
        output_path.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as log:
            result = subprocess.run(cmd, stdout=log, stderr=log, text=True)

    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    return protein_id, result.returncode == 0, time.time() - t0


# ── Merge ─────────────────────────────────────────────────────────────────────

_SEQ_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|(\S+?)_(\d+-\d+)$")


def _parse_seq_id(seq_id: str):
    m = _SEQ_RE.match(seq_id)
    return (m.group(1), m.group(2), m.group(3)) if m else (seq_id, seq_id, "")


def _merge_results(cohort_dir: Path, total_initial: int = 0, n_with_regions: int = 0):
    """Collect all per-protein biomarker_candidates.csv and write master files."""
    candidate_files = sorted(
        p for p in cohort_dir.rglob("step5_scores/biomarker_candidates.csv")
    )

    if not candidate_files:
        print("  [merge] No biomarker_candidates.csv files found.")
        return pd.DataFrame(), pd.DataFrame()

    MUTATION_THRESHOLD = 0.475   # mirrors run_pipeline.py constant

    frames = []
    migrated = 0
    for csv_path in candidate_files:
        df = pd.read_csv(csv_path)
        # Auto-migrate old window-level format (has mutation_sum, lacks contrib_total)
        if "mutation_sum" in df.columns and "contrib_total" not in df.columns:
            residue_csv = csv_path.parent / "per_residue_contributions.csv"
            if residue_csv.exists():
                residues = pd.read_csv(residue_csv)
                df = residues[residues["contrib_total"] >= MUTATION_THRESHOLD].copy()
                df = df.rename(columns={"region_id": "sequence_id", "abs_position": "position"})
                df = df.sort_values("contrib_total", ascending=False).reset_index(drop=True)
                migrated += 1
            else:
                continue  # can't migrate without residue data — skip
        if not df.empty:
            frames.append(df)

    if migrated:
        print(f"  [merge] Auto-migrated {migrated} protein(s) from old window-level format")

    if not frames:
        print("  [merge] All candidate files are empty.")
        return pd.DataFrame(), pd.DataFrame()

    master = pd.concat(frames, ignore_index=True)

    parsed = master["sequence_id"].apply(
        lambda s: pd.Series(_parse_seq_id(s), index=["uniprot_id", "protein_name", "region"])
    )
    master = pd.concat([master, parsed], axis=1)

    col_order = [
        "uniprot_id", "protein_name", "region",
        "sequence_id", "position", "aa",
        "contrib_L", "contrib_B", "contrib_total",
    ]
    master = master[[c for c in col_order if c in master.columns]]
    master = master.sort_values("contrib_total", ascending=False).reset_index(drop=True)

    summary_rows = []
    for (uid, pname), grp in master.groupby(["uniprot_id", "protein_name"], sort=False):
        top = grp.iloc[0]
        summary_rows.append({
            "uniprot_id":        uid,
            "protein_name":      pname,
            "n_candidates":      len(grp),
            "n_regions":         grp["region"].nunique(),
            "top_position":      top["position"],
            "top_aa":            top["aa"],
            "top_region":        top["region"],
            "top_contrib_total": round(top["contrib_total"], 4),
        })

    summary = (
        pd.DataFrame(summary_rows)
        .sort_values("top_contrib_total", ascending=False)
        .reset_index(drop=True)
    )

    master_path  = cohort_dir / "master_biomarker_candidates.csv"
    summary_path = cohort_dir / "summary.csv"
    master.to_csv(master_path,  index=False)
    summary.to_csv(summary_path, index=False)

    # ── Funnel counts ─────────────────────────────────────────────────────────
    # Count proteins that had at least one region detected (non-empty step2 CSV)
    n_with_regions = 0
    for region_csv in cohort_dir.glob("*/step2_regions/low_cleavage_regions.csv"):
        try:
            if region_csv.stat().st_size > 50:   # header-only file is ~40 bytes
                n_with_regions += 1
        except OSError:
            pass

    n_with_candidates = len(summary)
    n_total_regions   = int(summary["n_regions"].sum())
    n_total_residues  = len(master)

    funnel_rows = [
        {"stage": "Initial proteins",               "count": total_initial,      "pct_of_initial": 100.0 if total_initial else None},
        {"stage": "Proteins with regions (step 2)", "count": n_with_regions,     "pct_of_initial": round(n_with_regions / total_initial * 100, 1) if total_initial else None},
        {"stage": "Proteins with candidates",       "count": n_with_candidates,  "pct_of_initial": round(n_with_candidates / total_initial * 100, 1) if total_initial else None},
        {"stage": "Total candidate regions",        "count": n_total_regions,    "pct_of_initial": None},
        {"stage": "Total candidate residues",       "count": n_total_residues,   "pct_of_initial": None},
    ]
    funnel_path = cohort_dir / "funnel.csv"
    pd.DataFrame(funnel_rows).to_csv(funnel_path, index=False)

    # ── Print full funnel ─────────────────────────────────────────────────────
    def _pct(n):
        return f"  ({n / total_initial * 100:.1f} %)" if total_initial else ""

    print(f"\n  {'='*45}")
    print(f"  Cohort funnel")
    print(f"  {'='*45}")
    if total_initial:
        print(f"  Initial proteins         : {total_initial:>6}")
    if n_with_regions:
        print(f"  With regions  (step 2)   : {n_with_regions:>6}{_pct(n_with_regions)}")
    print(f"  With candidates (step 5) : {n_with_candidates:>6}{_pct(n_with_candidates)}")
    print(f"  Candidate regions        : {n_total_regions:>6}")
    print(f"  Candidate residues       : {n_total_residues:>6}")
    print(f"  {'='*45}")

    print(f"\n  Top 10 residue candidates:")
    print(f"  {'#':>3}  {'UniProt':>10}  {'Protein':<16}  {'Pos':>6}  {'AA':>3}  {'Total':>8}")
    print(f"  {'-'*55}")
    for i, row in master.head(10).iterrows():
        print(f"  {i+1:>3}  {row['uniprot_id']:>10}  {row['protein_name']:<16}  "
              f"{row['position']:>6}  {row['aa']:>3}  {row['contrib_total']:>8.4f}")
    print(f"\n  Saved: {master_path}")
    print(f"  Saved: {summary_path}")
    print(f"  Saved: {funnel_path}")

    return master, summary


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run the biomarker pipeline on every protein in a FASTA, in parallel."
    )
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument(
        "--parallel", type=int, default=5,
        help="How many proteins to run simultaneously (default: 5)",
    )
    parser.add_argument(
        "--cores", type=int, default=20,
        help="Total CPU cores, divided equally across --parallel jobs (default: 20)",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Only process the first N proteins (useful for testing)",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output directory under results/ (default: cohort_<stem>_YYYY-MM-DD_HHMM)",
    )
    parser.add_argument(
        "--merge-only", action="store_true", default=False,
        help="Skip all pipeline runs; only (re-)merge existing protein outputs",
    )
    parser.add_argument(
        "--step1-only", action="store_true", default=False,
        help="Run Phase 1 (Step 1 cohort predictions) and exit — "
             "useful for pre-computing step1 CSVs before running Steps 2-5",
    )
    args = parser.parse_args()

    fasta_file = Path(args.fasta)
    if not fasta_file.exists():
        sys.exit(f"ERROR: FASTA not found: {fasta_file}")

    timestamp   = datetime.datetime.now().strftime("%Y-%m-%d_%H%M")
    cohort_name = args.output or f"cohort_{fasta_file.stem}_{timestamp}"
    cohort_dir  = Path("results") / cohort_name
    cohort_dir.mkdir(parents=True, exist_ok=True)

    parallel      = args.parallel
    total_cores   = args.cores
    cores_per_job = max(1, total_cores // parallel)
    pipeline_script = Path(__file__).parent / "run_pipeline.py"

    # Parse input FASTA — just id + sequence strings (no BioPython objects across processes)
    all_proteins = [
        (record.id, str(record.seq))
        for record in SeqIO.parse(fasta_file, "fasta")
    ]
    if args.limit:
        all_proteins = all_proteins[:args.limit]
    total_proteins = len(all_proteins)

    print(f"\n{'='*65}")
    print(f" Cohort run : {cohort_name}")
    print(f" Input      : {fasta_file}  ({total_proteins} protein(s))")
    print(f" Step 1     : batches of {STEP1_BATCH_SIZE}  |  CathD+Legumain parallel  "
          f"|  {max(1, total_cores // 2)} cores each")
    print(f" Steps 2-5  : {parallel} protein(s) at a time  |  {cores_per_job} cores each")
    print(f" Output     : {cohort_dir}/")
    if args.limit:
        print(f" Limit      : first {args.limit} proteins")
    print(f"{'='*65}")

    # ── Phase 1: Cohort Step 1 (CathD + Legumain on all proteins) ────────────
    if not args.merge_only:
        print(f"\n[Phase 1]  Step 1 — CathD + Legumain on {total_proteins} protein(s) "
              f"in batches of {STEP1_BATCH_SIZE}")
        _run_step1_cohort(
            fasta_file, cohort_dir, all_proteins,
            batch_size=STEP1_BATCH_SIZE,
            total_cores=total_cores,
        )

    if args.step1_only:
        print(f"\n[--step1-only]  Phase 1 complete — exiting before Steps 2-5.")
        print(f"  To continue, run without --step1-only (same --output dir):\n"
              f"    python run_cohort.py --fasta {fasta_file} "
              f"--output {cohort_name} --parallel {parallel} --cores {total_cores}")
        return

    # ── Build job list (skip already-complete proteins) ───────────────────────
    job_args   = []
    skipped    = []
    proto_dirs = []   # for progress monitor

    for protein_id, sequence in all_proteins:
        safe_id    = protein_id.replace("/", "_").replace("\\", "_")
        out_path   = cohort_dir / safe_id
        candidates = out_path / "step5_scores" / "biomarker_candidates.csv"
        short_name = protein_id.split("|")[-1] if "|" in protein_id else protein_id

        proto_dirs.append((short_name, out_path))

        if candidates.exists() and candidates.stat().st_size > 0:
            skipped.append(protein_id)
        elif not args.merge_only:
            job_args.append((protein_id, sequence, out_path, cores_per_job, pipeline_script))

    if skipped:
        print(f"\n[Resume] {len(skipped)} protein(s) already complete — skipping")

    # ── Run ───────────────────────────────────────────────────────────────────
    completed_counter = [len(skipped)]

    if job_args and not args.merge_only:
        todo = len(job_args)
        print(f"\n[Phase 2]  Steps 2-5 — {todo} protein(s)  "
              f"[{parallel} parallel  |  {cores_per_job} cores each]")
        print(f"      Logs  →  <cohort_dir>/<protein_id>/pipeline.log")
        t_run = time.time()

        stop_event = threading.Event()

        timing_window: deque = deque(maxlen=20)
        n_ok = n_fail = 0

        with ProcessPoolExecutor(max_workers=parallel) as ex:
            futures = {ex.submit(_run_protein, a): a[0] for a in job_args}
            for fut in as_completed(futures):
                protein_id, success, elapsed = fut.result()
                completed_counter[0] += 1
                timing_window.append(elapsed)
                done  = completed_counter[0]
                avg   = sum(timing_window) / len(timing_window)
                eta   = _eta(done, total_proteins, t_run)
                short = protein_id.split("|")[-1] if "|" in protein_id else protein_id
                status = "done" if success else "FAILED"
                print(
                    f"  [{done:>3}/{total_proteins}]  {short:<32}  "
                    f"{status}  ({elapsed/60:.1f}min)  "
                    f"avg {avg/60:.1f}min/prot  ETA {eta}",
                    flush=True,
                )
                if success:
                    n_ok += 1
                else:
                    n_fail += 1

        stop_event.set()

        print(f"\n  Complete: {n_ok + len(skipped)}/{total_proteins}"
              + (f"  ({n_fail} FAILED — check <protein>/pipeline.log)" if n_fail else ""))
        print(f"  Total elapsed: {_elapsed(t_run)}")

    elif args.merge_only:
        print("\n[Phase 2] Skipping pipeline runs (--merge-only)")
    else:
        print("\n[Phase 2] All proteins already complete — skipping to merge")

    # ── Merge ─────────────────────────────────────────────────────────────────
    print(f"\n[Merge] Collecting results…")
    t_merge = time.time()
    master, summary = _merge_results(cohort_dir, total_initial=total_proteins)
    print(f"  → done in {_elapsed(t_merge)}")

    print(f"\n{'='*65}")
    print(f" Done!  {len(summary)} protein(s) with candidates")
    print(f" Master : {cohort_dir}/master_biomarker_candidates.csv")
    print(f" Summary: {cohort_dir}/summary.csv")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
