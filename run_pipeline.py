#!/usr/bin/env python3
"""
Biomarker Discovery Pipeline
==============================
Criteria defined in pipeline.MD (last updated 12/23/2025).

Steps
-----
1. Run protease cleavage predictions (Cathepsin D + Legumain) on the input FASTA.
2. Find protein regions that are resistant to cleavage by BOTH proteases.
3. Create mutant FASTAs: S,T→P (structural mimetic) and S,T→E (charge mimetic),
   each extended by EXTENSION residues on each side of the detected region.
4. Run Cathepsin L (C01.032) and Cathepsin B (C01.060) predictions on original
   and mutant sequences.
5. Calculate mutation sums:
       sum = (CathL_original − CathL_ST→P) + (CathB_original − CathB_ST→E)
   and report positions where sum ≥ MUTATION_THRESHOLD as biomarker candidates.

Usage
-----
    conda activate prosperousplus
    cd /shared/macdata/groups/ppc/projects/ProsperousPlus
    python run_pipeline.py                          # runs on tau (default)
    python run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta
    python run_pipeline.py --fasta data/original_files/tau_only.fasta --batch-size 1

All results are written to:
    results/<input_stem>_YYYY-MM-DD_HHMM/
        step1_predictions/
            cathepsinD_cleavage.csv
            legumain_cleavage.csv
        step2_regions/
            low_cleavage_regions.csv
        step3_fastas/
            original.fasta
            mutant_ST_to_P.fasta
            mutant_ST_to_E.fasta
        step4_mutant_predictions/
            cathepsinL_original.csv
            cathepsinL_ST_to_P.csv
            cathepsinB_original.csv
            cathepsinB_ST_to_E.csv
        step5_scores/
            mutation_scores.csv
            biomarker_candidates.csv

Completed steps are detected automatically so interrupted runs resume from
where they left off (just re-run the same command).
"""

import argparse
import datetime
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

import pandas as pd  # noqa: E402 (requires prosperousplus conda env)

sys.path.insert(0, str(Path(__file__).parent))

from lib import (
    run_prediction, run_prediction_batched,
    find_low_cleavage_regions,
    load_sequences, create_mutation_fastas,
    calculate_mutation_sum, find_candidates,
)

# ── Protease model IDs ────────────────────────────────────────────────────────
CATHEPSIN_D = "A01.009"   # region detection
LEGUMAIN    = "C13.004"   # region detection
CATHEPSIN_L = "C01.032"   # mutation analysis (S,T→P)
CATHEPSIN_B = "C01.060"   # mutation analysis (S,T→E)

# ── Pipeline criteria (pipeline.MD) ──────────────────────────────────────────
MIN_LENGTH         = 25     # minimum consecutive low-cleavage residues
MAX_EXCEPTIONS     = 2      # total high-scoring positions allowed across both proteases
SCORE_THRESHOLD    = 0.3    # positions above this count as "exceptions"
EXTENSION          = 5      # residues to add on each side of detected region
MUTATION_THRESHOLD = 0.475  # minimum mutation sum to call a candidate (= pTau217 level)


def _csv_ready(path: Path) -> bool:
    """Return True if a CSV file already exists and is non-empty."""
    return path.exists() and path.stat().st_size > 0


def predict(fasta_file: Path, protease: str, out_csv: Path,
            out_dir: Path, batch_size: int, workers: Optional[int],
            batch_dir: Optional[Path] = None) -> pd.DataFrame:
    """Run prediction, skipping if out_csv already exists."""
    if _csv_ready(out_csv):
        print(f"  [resume] {out_csv.name} already exists — skipping")
        return pd.read_csv(out_csv)

    if batch_size > 1:
        df = run_prediction_batched(fasta_file, protease, out_dir,
                                    batch_size=batch_size, workers=workers,
                                    batch_dir=batch_dir)
    else:
        df = run_prediction(fasta_file, protease, out_dir)

    df.to_csv(out_csv, index=False)
    return df


def main():
    parser = argparse.ArgumentParser(description="ProsperousPlus Biomarker Pipeline")
    parser.add_argument(
        "--fasta", default="data/original_files/tau_only.fasta",
        help="Input FASTA file (default: tau)",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output directory name under results/ (default: <fasta_stem>_YYYY-MM-DD_HHMM)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=50,
        help="Sequences per batch for large inputs (use 1 to disable batching)",
    )
    parser.add_argument(
        "--workers", type=int, default=None,
        help="Parallel worker processes (default: cpu_count - 1)",
    )
    args = parser.parse_args()

    fasta_file = Path(args.fasta)
    if not fasta_file.exists():
        sys.exit(f"ERROR: FASTA file not found: {fasta_file}")

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M")  # e.g. 2026-04-23_1430
    run_name  = args.output or f"{fasta_file.stem}_{timestamp}"
    out_root = Path("results") / run_name
    batch_size = args.batch_size
    workers    = args.workers

    # ── Output paths ─────────────────────────────────────────────────────────
    pred_dir    = out_root / "step1_predictions"
    regions_dir = out_root / "step2_regions"
    fastas_dir  = out_root / "step3_fastas"
    mutpred_dir = out_root / "step4_mutant_predictions"
    scores_dir  = out_root / "step5_scores"

    for d in [pred_dir, regions_dir, fastas_dir, mutpred_dir, scores_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Shared batch-split dirs: one split per unique FASTA, reused across proteases
    # that run on the same input — avoids writing duplicate batch files.
    step1_batches   = pred_dir    / "_batches"          # input FASTA  → CathD + Legumain
    orig_batches    = mutpred_dir / "_batches_original"  # original.fasta → CathL + CathB
    mutp_batches    = mutpred_dir / "_batches_mutp"      # mutant_ST_to_P → CathL
    mute_batches    = mutpred_dir / "_batches_mute"      # mutant_ST_to_E → CathB

    cathD_csv    = pred_dir / "cathepsinD_cleavage.csv"
    legumain_csv = pred_dir / "legumain_cleavage.csv"
    regions_csv  = regions_dir / "low_cleavage_regions.csv"
    orig_fasta   = fastas_dir / "original.fasta"
    mutp_fasta   = fastas_dir / "mutant_ST_to_P.fasta"
    mute_fasta   = fastas_dir / "mutant_ST_to_E.fasta"
    cathL_orig_csv  = mutpred_dir / "cathepsinL_original.csv"
    cathL_mutp_csv  = mutpred_dir / "cathepsinL_ST_to_P.csv"
    cathB_orig_csv  = mutpred_dir / "cathepsinB_original.csv"
    cathB_mute_csv  = mutpred_dir / "cathepsinB_ST_to_E.csv"
    scores_csv     = scores_dir / "mutation_scores.csv"
    candidates_csv = scores_dir / "biomarker_candidates.csv"

    # ── Step 1: Cleavage predictions for region detection (parallel) ─────────
    print(f"\n{'='*60}")
    print(f" Input:      {fasta_file}")
    print(f" Output:     {out_root}")
    print(f" Batch size: {batch_size}")
    print(f"{'='*60}")
    print("\nStep 1: Running cleavage predictions (Cathepsin D + Legumain) [parallel]")

    with ThreadPoolExecutor(max_workers=2) as ex:
        fut_D = ex.submit(predict, fasta_file, CATHEPSIN_D, cathD_csv,
                          pred_dir / "cathepsinD", batch_size, workers, step1_batches)
        fut_L = ex.submit(predict, fasta_file, LEGUMAIN,    legumain_csv,
                          pred_dir / "legumain",   batch_size, workers, step1_batches)
        pred_D = fut_D.result()
        pred_L = fut_L.result()

    print(f"  Cathepsin D:  {len(pred_D):,} positions")
    print(f"  Legumain:     {len(pred_L):,} positions")

    # ── Step 2: Find low-cleavage regions ────────────────────────────────────
    print("\nStep 2: Finding low-cleavage regions")
    if _csv_ready(regions_csv):
        print("  [resume] low_cleavage_regions.csv already exists — skipping")
        regions = pd.read_csv(regions_csv)
    else:
        regions = find_low_cleavage_regions(
            pred_D, pred_L,
            CATHEPSIN_D, LEGUMAIN,
            min_length=MIN_LENGTH,
            max_exceptions=MAX_EXCEPTIONS,
            threshold=SCORE_THRESHOLD,
        )
        regions.to_csv(regions_csv, index=False)

    print(f"  Found {len(regions)} region(s):")
    for _, r in regions.iterrows():
        print(f"    {r['sequence_id']}: {r['start_position']}–{r['end_position']} ({r['length']} aa)")

    if len(regions) == 0:
        print("\nNo regions found — nothing to mutate. Done.")
        return

    # ── Step 3: Create mutant FASTAs ─────────────────────────────────────────
    print("\nStep 3: Creating mutant FASTAs")
    if orig_fasta.exists() and mutp_fasta.exists() and mute_fasta.exists():
        print("  [resume] FASTAs already exist — skipping")
    else:
        full_sequences = load_sequences(fasta_file)
        create_mutation_fastas(
            regions, full_sequences, fastas_dir, EXTENSION,
            orig_name="original.fasta",
            mutp_name="mutant_ST_to_P.fasta",
            mute_name="mutant_ST_to_E.fasta",
        )
    print(f"  {fastas_dir}/original.fasta")
    print(f"  {fastas_dir}/mutant_ST_to_P.fasta")
    print(f"  {fastas_dir}/mutant_ST_to_E.fasta")

    # ── Step 4: Mutant predictions (all 4 in parallel) ───────────────────────
    print("\nStep 4: Running mutant predictions [parallel]")

    with ThreadPoolExecutor(max_workers=4) as ex:
        fut_orig_L = ex.submit(predict, orig_fasta, CATHEPSIN_L, cathL_orig_csv,
                               mutpred_dir / "cathL_original", batch_size, workers, orig_batches)
        fut_mutp_L = ex.submit(predict, mutp_fasta, CATHEPSIN_L, cathL_mutp_csv,
                               mutpred_dir / "cathL_mutp",     batch_size, workers, mutp_batches)
        fut_orig_B = ex.submit(predict, orig_fasta, CATHEPSIN_B, cathB_orig_csv,
                               mutpred_dir / "cathB_original", batch_size, workers, orig_batches)
        fut_mute_B = ex.submit(predict, mute_fasta, CATHEPSIN_B, cathB_mute_csv,
                               mutpred_dir / "cathB_mute",     batch_size, workers, mute_batches)
        orig_L = fut_orig_L.result()
        mutp_L = fut_mutp_L.result()
        orig_B = fut_orig_B.result()
        mute_B = fut_mute_B.result()

    print(f"  CathL original:  {len(orig_L):,} positions")
    print(f"  CathL ST→P:      {len(mutp_L):,} positions")
    print(f"  CathB original:  {len(orig_B):,} positions")
    print(f"  CathB ST→E:      {len(mute_B):,} positions")

    # ── Step 5: Mutation sums and candidates ─────────────────────────────────
    print("\nStep 5: Calculating mutation sums")
    if _csv_ready(scores_csv) and _csv_ready(candidates_csv):
        print("  [resume] mutation_scores.csv already exists — skipping")
        scores     = pd.read_csv(scores_csv)
        candidates = pd.read_csv(candidates_csv)
    else:
        scores = calculate_mutation_sum(orig_L, mutp_L, orig_B, mute_B)
        scores.to_csv(scores_csv, index=False)

        candidates = find_candidates(scores, MUTATION_THRESHOLD)
        candidates.to_csv(candidates_csv, index=False)

    print(f"  Total positions scored: {len(scores):,}")
    print(f"  Candidates (sum ≥ {MUTATION_THRESHOLD}): {len(candidates)}")

    if len(candidates) > 0:
        print("\n  Top biomarker candidates:")
        print(f"  {'Position':>10}  {'Peptide':<10}  {'Mut Sum':>8}  {'CathL diff':>10}  {'CathB diff':>10}")
        print(f"  {'-'*55}")
        for _, c in candidates.head(15).iterrows():
            print(f"  {c['position']:>10}  {c['peptide']:<10}  {c['mutation_sum']:>8.4f}"
                  f"  {c['cathL_diff']:>10.4f}  {c['cathB_diff']:>10.4f}")

    print(f"\nDone. All results in: {out_root}/")
    print(f"  Key outputs:")
    print(f"    {candidates_csv}")
    print(f"    {scores_csv}")


if __name__ == "__main__":
    main()
