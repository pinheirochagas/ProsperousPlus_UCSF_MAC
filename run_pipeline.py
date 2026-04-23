#!/usr/bin/env python3
"""
Biomarker Discovery Pipeline
==============================
Criteria defined in pipeline.MD (last updated 12/23/2025).

Steps
-----
1. Run protease cleavage predictions (Cathepsin D + Legumain) on the input FASTA.
2. Find protein regions that are resistant to cleavage by BOTH proteases.
3. Create single-mutant FASTAs: for each S or T in each extended region, one
   sequence is produced with only that residue changed (S/T→P for CathL,
   S/T→E for CathB).  Each extended region is padded by EXTENSION residues.
4. Run Cathepsin L (C01.032) and Cathepsin B (C01.060) predictions on original
   and single-mutant sequences (4 jobs in parallel).
5. Calculate per-residue contributions and window-level mutation sums:
       contrib_k = (CathL_orig[k] − CathL_single_P_k[k])
                 + (CathB_orig[k] − CathB_single_E_k[k])
       mutation_sum[window] = Σ contrib_k  for S/T positions k within that window
   Report windows where mutation_sum ≥ MUTATION_THRESHOLD as biomarker candidates.

Usage
-----
    conda activate prosperousplus
    cd /shared/macdata/groups/ppc/projects/ProsperousPlus
    python run_pipeline.py                          # runs on tau (default)
    python run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta
    python run_pipeline.py --fasta data/original_files/tau_only.fasta --batch-size 1
    python run_pipeline.py --cores 20               # default; adjust to available CPUs

All results are written to:
    results/<input_stem>_YYYY-MM-DD_HHMM/
        step1_predictions/
            cathepsinD_cleavage.csv
            legumain_cleavage.csv
        step2_regions/
            low_cleavage_regions.csv
        step3_fastas/
            original.fasta
            single_ST_to_P.fasta   ← one sequence per S/T per region (→P)
            single_ST_to_E.fasta   ← one sequence per S/T per region (→E)
        step4_mutant_predictions/
            cathepsinL_original.csv
            cathepsinL_single_P.csv
            cathepsinB_original.csv
            cathepsinB_single_E.csv
        step5_scores/
            mutation_scores.csv
            per_residue_contributions.csv
            biomarker_candidates.csv

Completed steps are detected automatically so interrupted runs resume from
where they left off (just re-run the same command).
"""

import argparse
import datetime
import time
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

import pandas as pd  # noqa: E402 (requires prosperousplus conda env)

sys.path.insert(0, str(Path(__file__).parent))

from lib import (
    run_prediction, run_prediction_batched,
    find_low_cleavage_regions,
    load_sequences, create_single_mutant_fastas,
    calculate_single_residue_mutation_sum, find_candidates,
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
            out_dir: Path, batch_size: int, process_num: int,
            batch_dir: Optional[Path] = None) -> pd.DataFrame:
    """Run prediction, skipping if out_csv already exists."""
    if _csv_ready(out_csv):
        return pd.read_csv(out_csv)

    if batch_size > 1:
        df = run_prediction_batched(fasta_file, protease, out_dir,
                                    batch_size=batch_size,
                                    process_num=process_num,
                                    batch_dir=batch_dir,
                                    verbose=False)
    else:
        df = run_prediction(fasta_file, protease, out_dir,
                            process_num=process_num)

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
        "--cores", type=int, default=20,
        help="Total CPU cores to use, divided evenly across concurrent jobs (default: 20)",
    )
    args = parser.parse_args()

    fasta_file = Path(args.fasta)
    if not fasta_file.exists():
        sys.exit(f"ERROR: FASTA file not found: {fasta_file}")

    timestamp  = datetime.datetime.now().strftime("%Y-%m-%d_%H%M")  # e.g. 2026-04-23_1430
    run_name   = args.output or f"{fasta_file.stem}_{timestamp}"
    out_root   = Path("results") / run_name
    batch_size = args.batch_size
    cores      = args.cores

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
    step1_batches   = pred_dir    / "_batches"           # input FASTA   → CathD + Legumain
    orig_batches    = mutpred_dir / "_batches_original"  # original.fasta → CathL + CathB
    mutp_batches    = mutpred_dir / "_batches_single_P"  # single_ST_to_P → CathL
    mute_batches    = mutpred_dir / "_batches_single_E"  # single_ST_to_E → CathB

    cathD_csv        = pred_dir    / "cathepsinD_cleavage.csv"
    legumain_csv     = pred_dir    / "legumain_cleavage.csv"
    regions_csv      = regions_dir / "low_cleavage_regions.csv"
    orig_fasta       = fastas_dir  / "original.fasta"
    single_P_fasta   = fastas_dir  / "single_ST_to_P.fasta"
    single_E_fasta   = fastas_dir  / "single_ST_to_E.fasta"
    cathL_orig_csv   = mutpred_dir / "cathepsinL_original.csv"
    cathL_single_csv = mutpred_dir / "cathepsinL_single_P.csv"
    cathB_orig_csv   = mutpred_dir / "cathepsinB_original.csv"
    cathB_single_csv = mutpred_dir / "cathepsinB_single_E.csv"
    scores_csv       = scores_dir  / "mutation_scores.csv"
    residue_csv      = scores_dir  / "per_residue_contributions.csv"
    candidates_csv   = scores_dir  / "biomarker_candidates.csv"

    pn_step1 = max(1, cores // 2)   # 2 concurrent jobs in Step 1
    pn_step4 = max(1, cores // 4)   # 4 concurrent jobs in Step 4

    pipeline_start = time.time()

    def _elapsed(t0):
        s = time.time() - t0
        return f"{s:.1f}s" if s < 60 else f"{s/60:.1f}min"

    def _step(n, title):
        print(f"\n[Step {n}] {title}")

    def _ok(t0, **stats):
        parts = "  ".join(f"{k}: {v}" for k, v in stats.items())
        print(f"  {parts}")
        print(f"  done in {_elapsed(t0)}")

    # Count input sequences
    from Bio import SeqIO as _SeqIO
    n_input_seqs = sum(1 for _ in _SeqIO.parse(fasta_file, "fasta"))

    print(f"\n{'='*60}")
    print(f" Input:      {fasta_file}  ({n_input_seqs} sequence(s))")
    print(f" Output:     {out_root}")
    print(f" Batch size: {batch_size}")
    print(f" Cores:      {cores} total  ({pn_step1}/job in Step 1, {pn_step4}/job in Step 4)")
    print(f" Started:    {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")

    # ── Step 1 ───────────────────────────────────────────────────────────────
    t0 = time.time()
    print("\n[Step 1] Cleavage predictions (Cathepsin D + Legumain) [2 parallel]")

    with ThreadPoolExecutor(max_workers=2) as ex:
        fut_D = ex.submit(predict, fasta_file, CATHEPSIN_D, cathD_csv,
                          pred_dir / "cathepsinD", batch_size, pn_step1, step1_batches)
        fut_L = ex.submit(predict, fasta_file, LEGUMAIN,    legumain_csv,
                          pred_dir / "legumain",   batch_size, pn_step1, step1_batches)
        pred_D = fut_D.result()
        pred_L = fut_L.result()

    print(f"  Cathepsin D:  {len(pred_D):,} window scores  ({pred_D['sequence_id'].nunique()} sequences)")
    print(f"  Legumain:     {len(pred_L):,} window scores  ({pred_L['sequence_id'].nunique()} sequences)")
    print(f"  → done in {_elapsed(t0)}")

    # ── Step 2 ───────────────────────────────────────────────────────────────
    t0 = time.time()
    print("\n[Step 2] Finding low-cleavage regions")
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

    print(f"  Found {len(regions)} protected region(s)  "
          f"(min_length={MIN_LENGTH} aa, max_exceptions={MAX_EXCEPTIONS}, threshold={SCORE_THRESHOLD}):")
    for _, r in regions.iterrows():
        print(f"    {r['sequence_id']}: pos {r['start_position']}–{r['end_position']} "
              f"({r['length']} aa)")
    print(f"  → done in {_elapsed(t0)}")

    if len(regions) == 0:
        print("\nNo regions found — nothing to mutate. Done.")
        return

    # ── Step 3 ───────────────────────────────────────────────────────────────
    t0 = time.time()
    print("\n[Step 3] Creating single-mutant FASTAs")
    full_sequences = load_sequences(fasta_file)
    already_exist  = orig_fasta.exists() and single_P_fasta.exists() and single_E_fasta.exists()

    st_positions = create_single_mutant_fastas(
        regions, full_sequences, fastas_dir, EXTENSION,
        orig_name="original.fasta",
        mutp_name="single_ST_to_P.fasta",
        mute_name="single_ST_to_E.fasta",
    )
    total_st  = sum(len(v) for v in st_positions.values())
    n_regions = len(st_positions)

    if already_exist:
        print("  [resume] FASTAs already exist")
    per_region = {rid: len(v) for rid, v in st_positions.items()}
    for rid, n in per_region.items():
        print(f"    {rid}: {n} mutable S/T residues")
    print(f"  Total: {n_regions} region(s), {total_st} single-mutant sequences per protease")
    print(f"  single_ST_to_P.fasta  ({total_st} seqs × S/T→P)")
    print(f"  single_ST_to_E.fasta  ({total_st} seqs × S/T→E)")
    print(f"  → done in {_elapsed(t0)}")

    # ── Step 4 ───────────────────────────────────────────────────────────────
    t0 = time.time()
    print(f"\n[Step 4] Predictions on original + single mutants "
          f"[4 parallel, {pn_step4} cores each]")

    with ThreadPoolExecutor(max_workers=4) as ex:
        fut_orig_L   = ex.submit(predict, orig_fasta,     CATHEPSIN_L, cathL_orig_csv,
                                 mutpred_dir / "cathL_original", batch_size, pn_step4, orig_batches)
        fut_single_L = ex.submit(predict, single_P_fasta, CATHEPSIN_L, cathL_single_csv,
                                 mutpred_dir / "cathL_single_P", batch_size, pn_step4, mutp_batches)
        fut_orig_B   = ex.submit(predict, orig_fasta,     CATHEPSIN_B, cathB_orig_csv,
                                 mutpred_dir / "cathB_original", batch_size, pn_step4, orig_batches)
        fut_single_B = ex.submit(predict, single_E_fasta, CATHEPSIN_B, cathB_single_csv,
                                 mutpred_dir / "cathB_single_E", batch_size, pn_step4, mute_batches)
        orig_L   = fut_orig_L.result()
        single_L = fut_single_L.result()
        orig_B   = fut_orig_B.result()
        single_E = fut_single_B.result()

    orig_seqs    = orig_L['sequence_id'].nunique()
    single_seqs  = single_L['sequence_id'].nunique()
    print(f"  CathL original:   {len(orig_L):,} rows  ({orig_seqs} sequences)")
    print(f"  CathL single-P:   {len(single_L):,} rows  ({single_seqs} mutant sequences × windows)")
    print(f"  CathB original:   {len(orig_B):,} rows  ({orig_B['sequence_id'].nunique()} sequences)")
    print(f"  CathB single-E:   {len(single_E):,} rows  ({single_E['sequence_id'].nunique()} mutant sequences × windows)")
    print(f"  → done in {_elapsed(t0)}")

    # ── Step 5 ───────────────────────────────────────────────────────────────
    t0 = time.time()
    print("\n[Step 5] Per-residue contributions and mutation sums")
    if _csv_ready(scores_csv) and _csv_ready(candidates_csv):
        print("  [resume] mutation_scores.csv already exists — skipping")
        scores     = pd.read_csv(scores_csv)
        candidates = pd.read_csv(candidates_csv)
    else:
        scores, residues = calculate_single_residue_mutation_sum(
            orig_L, single_L, orig_B, single_E, st_positions
        )
        scores.to_csv(scores_csv, index=False)
        residues.to_csv(residue_csv, index=False)

        candidates = find_candidates(scores, MUTATION_THRESHOLD)
        candidates.to_csv(candidates_csv, index=False)

    n_above = (scores['mutation_sum'] > 0).sum()
    print(f"  Window positions scored:  {len(scores):,}")
    print(f"  S/T residues contributing: {total_st}")
    print(f"  Windows with sum > 0:     {n_above:,}")
    print(f"  Candidates (sum ≥ {MUTATION_THRESHOLD}): {len(candidates)}")
    print(f"  → done in {_elapsed(t0)}")

    if len(candidates) > 0:
        print("\n  Top biomarker candidates:")
        print(f"  {'#':>3}  {'Position':>8}  {'Peptide':<10}  {'Mut Sum':>8}  "
              f"{'CathL Δ':>8}  {'CathB Δ':>8}")
        print(f"  {'-'*58}")
        for i, (_, c) in enumerate(candidates.head(15).iterrows(), 1):
            print(f"  {i:>3}  {c['position']:>8}  {c['peptide']:<10}  "
                  f"{c['mutation_sum']:>8.4f}  {c['cathL_diff']:>8.4f}  {c['cathB_diff']:>8.4f}")

    total_elapsed = _elapsed(pipeline_start)
    print(f"\n{'='*60}")
    print(f" Done in {total_elapsed}  —  results in: {out_root}/")
    print(f"   {candidates_csv.name}  ({len(candidates)} candidates)")
    print(f"   {scores_csv.name}  ({len(scores):,} positions)")
    print(f"   {residue_csv.name}  ({total_st} S/T residues)")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
