"""Run ProsperousPlus predictions."""

import subprocess
import pandas as pd
import shutil
from pathlib import Path

CODE_ROOT = Path('/shared/macdata/groups/ppc/code/ProsperousPlus')


def run_prediction(fasta_file, protease, output_dir, process_num: int = 10, plot: bool = False):
    """Run prediction for a single protease.

    Parameters
    ----------
    process_num : int
        Number of CPU cores passed to ProsperousPlus via ``--processNum``.
    """
    fasta_path = Path(fasta_file).absolute()
    output_path = Path(output_dir).absolute()
    
    # Remove existing output dir to avoid conflicts
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        'python', str(CODE_ROOT / 'Prosperousplus.py'),
        '--predictfile', str(fasta_path),
        '--outputpath', str(output_path),
        '--inputType', 'fasta',
        '--protease', protease,
        '--mode', 'prediction',
        '--processNum', str(process_num),
        '--PLOT', 'Yes' if plot else 'No'
    ]
    
    result = subprocess.run(cmd, cwd=CODE_ROOT, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr[:500]}")
        raise RuntimeError(f"Prediction failed for {protease}")
    
    return pd.read_csv(output_path / 'results.csv')


def load_predictions(results_path):
    """Load predictions from CSV."""
    return pd.read_csv(results_path)
