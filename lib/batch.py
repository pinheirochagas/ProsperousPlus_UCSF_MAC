"""Batched parallel prediction for large FASTA inputs.

For inputs with many sequences (e.g. brain proteome), this module splits
the FASTA into chunks, runs ProsperousPlus on each chunk in parallel, and
combines the results.  Completed chunks are detected automatically so
interrupted runs resume from where they left off.
"""

import subprocess
import multiprocessing as mp
from typing import List, Optional, Tuple

import pandas as pd
from pathlib import Path
from Bio import SeqIO

from .predict import CODE_ROOT


# ── batch size tuned for ProsperousPlus memory usage ──────────────────────────
DEFAULT_BATCH_SIZE = 50


def split_fasta(fasta_file: Path, batch_dir: Path, batch_size: int = DEFAULT_BATCH_SIZE) -> List[Path]:
    """Split a FASTA file into numbered batch files.

    Returns the list of batch file paths (already-existing batches are skipped).
    """
    fasta_path = Path(fasta_file)
    batch_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(fasta_path, "fasta"))
    batches: list[Path] = []

    for i in range(0, len(records), batch_size):
        idx = i // batch_size + 1
        batch_path = batch_dir / f"batch_{idx:04d}.fasta"
        if not batch_path.exists():
            SeqIO.write(records[i : i + batch_size], batch_path, "fasta")
        batches.append(batch_path)

    return batches


def _run_batch(args: Tuple) -> dict:
    """Worker: run ProsperousPlus on a single batch file.  Returns status dict."""
    batch_file, protease, out_dir = args
    batch_path = Path(batch_file)
    out_path = Path(out_dir)
    results_csv = out_path / "results.csv"

    # Resume: skip if already done
    if results_csv.exists() and results_csv.stat().st_size > 0:
        return {"batch": batch_path.name, "protease": protease, "status": "skipped"}

    out_path.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python", str(CODE_ROOT / "Prosperousplus.py"),
        "--predictfile", str(batch_path.absolute()),
        "--outputpath",  str(out_path.absolute()),
        "--inputType",   "fasta",
        "--protease",    protease,
        "--mode",        "prediction",
        "--processNum",  "10",
        "--PLOT",        "No",
    ]

    result = subprocess.run(cmd, cwd=CODE_ROOT, capture_output=True, text=True)

    if result.returncode != 0:
        return {
            "batch": batch_path.name,
            "protease": protease,
            "status": "failed",
            "error": result.stderr[:400],
        }

    return {"batch": batch_path.name, "protease": protease, "status": "ok"}


def _combine_batches(batch_out_dirs: List[Path]) -> pd.DataFrame:
    """Concatenate results.csv from each batch directory into one DataFrame."""
    frames = []
    for d in sorted(batch_out_dirs):
        csv = d / "results.csv"
        if csv.exists():
            frames.append(pd.read_csv(csv))
    if not frames:
        raise RuntimeError("No batch results found to combine.")
    return pd.concat(frames, ignore_index=True)


def run_prediction_batched(
    fasta_file,
    protease: str,
    output_dir,
    batch_size: int = DEFAULT_BATCH_SIZE,
    workers: Optional[int] = None,
    batch_dir=None,
) -> pd.DataFrame:
    """Run prediction for a single protease on a (potentially large) FASTA.

    Splits *fasta_file* into batches, processes them in parallel, combines
    the results, and returns a single DataFrame identical in schema to what
    ``run_prediction`` returns for small inputs.

    Parameters
    ----------
    fasta_file:  path to input FASTA
    protease:    protease model ID (e.g. "A01.009")
    output_dir:  directory where per-batch results are cached
    batch_size:  sequences per batch (default 50)
    workers:     parallel processes (default: cpu_count)
    batch_dir:   where to write the split FASTA chunks (default: output_dir/_batches).
                 Pass a shared directory to avoid splitting the same FASTA multiple
                 times when several proteases run on the same input.
    """
    fasta_path = Path(fasta_file)
    out_root = Path(output_dir)
    batch_dir = Path(batch_dir) if batch_dir is not None else out_root / "_batches"

    if workers is None:
        workers = max(1, mp.cpu_count() - 1)

    print(f"  [{protease}] splitting into batches of {batch_size}…")
    batch_files = split_fasta(fasta_path, batch_dir, batch_size)
    print(f"  [{protease}] {len(batch_files)} batches, {workers} workers")

    # Build work items; worker skips already-completed batches automatically
    work = [
        (str(bf), protease, str(out_root / f"batch_{i+1:04d}"))
        for i, bf in enumerate(batch_files)
    ]

    with mp.Pool(workers) as pool:
        results = pool.map(_run_batch, work)

    failed = [r for r in results if r["status"] == "failed"]
    if failed:
        for f in failed:
            print(f"  [!] {f['batch']} failed: {f.get('error','')}")
        raise RuntimeError(f"{len(failed)} batch(es) failed for protease {protease}")

    done = sum(1 for r in results if r["status"] == "ok")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    print(f"  [{protease}] done={done}, resumed={skipped}")

    batch_out_dirs = [out_root / f"batch_{i+1:04d}" for i in range(len(batch_files))]
    combined = _combine_batches(batch_out_dirs)

    # Save combined CSV alongside batches for quick reuse
    combined.to_csv(out_root / "results.csv", index=False)
    return combined
