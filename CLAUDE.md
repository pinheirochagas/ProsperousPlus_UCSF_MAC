# ProsperousPlus — Biomarker Discovery Workspace

## Repository layout

| Path | Role |
|------|------|
| `/shared/macdata/groups/ppc/code/ProsperousPlus/` | ProsperousPlus prediction tool (do not edit) |
| `/shared/macdata/groups/ppc/projects/ProsperousPlus/` | **This workspace** — analysis code, data, results |
| `/shared/macdata/groups/ppc/projects/ProsperousPlus_older/` | Archive of superseded scripts and old result trees |

On macOS the shared path mounts as `/Volumes/macdata/groups/ppc/…`.

---

## Quick start

```bash
conda activate prosperousplus
cd /shared/macdata/groups/ppc/projects/ProsperousPlus

# Run on tau (single protein, fast — ~2 min)
python run_pipeline.py

# Run on brain-enriched proteins (large dataset, batched)
python run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta
python run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enhanced.fasta
python run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Group_Group_enriched.fasta

# Resume an interrupted run (just re-run the same command — completed steps are skipped)
python run_pipeline.py --fasta data/brain_elevated/...same file...
```

---

## Pipeline criteria (`pipeline.MD`)

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Region-detection proteases | A01.009, C13.004 | Cathepsin D AND Legumain |
| `min_length` | 25 | Minimum consecutive low-cleavage residues |
| `max_exceptions` | 2 | Max high-scoring positions allowed across both proteases |
| `score_threshold` | 0.3 | Positions above this count as "exceptions" |
| `extension` | 5 | Residues added on each side of detected region for mutation analysis |
| Mutation proteases | C01.032 (CathL), C01.060 (CathB) | Run on original + mutant sequences |
| `mutation_threshold` | 0.475 | Min mutation sum to call a candidate (= pTau217 level) |

**Mutation formula:**
```
mutation_sum = (CathL_original − CathL_ST→P) + (CathB_original − CathB_ST→E)
```

---

## Pipeline steps and output files

Each run writes to `results/<input_stem>_YYYY-MM-DD_HHMM/`:

```
results/tau_only_2026-04-23_1430/        ← one folder per input FASTA + date + time
  step1_predictions/
    cathepsinD_cleavage.csv              # CathD predictions on full input
    legumain_cleavage.csv                # Legumain predictions on full input
  step2_regions/
    low_cleavage_regions.csv             # regions resistant to both proteases
  step3_fastas/
    original.fasta                       # extended region sequences
    mutant_ST_to_P.fasta                 # S,T→P (structural mimetic)
    mutant_ST_to_E.fasta                 # S,T→E (charge mimetic)
  step4_mutant_predictions/
    cathepsinL_original.csv              # CathL on original
    cathepsinL_ST_to_P.csv              # CathL on S,T→P mutants
    cathepsinB_original.csv              # CathB on original
    cathepsinB_ST_to_E.csv              # CathB on S,T→E mutants
  step5_scores/
    mutation_scores.csv                  # all positions with mutation sums
    biomarker_candidates.csv             # positions above 0.475 threshold
```

---

## Batching (for large inputs)

For inputs with many sequences, predictions are split into batches and run in parallel.

```bash
# Default: 50 sequences/batch, cpu_count-1 workers
python run_pipeline.py --fasta data/brain_elevated/...fasta

# Custom batch size and workers
python run_pipeline.py --fasta data/brain_elevated/...fasta --batch-size 25 --workers 8

# Disable batching (single-sequence files like tau)
python run_pipeline.py --fasta data/original_files/tau_only.fasta --batch-size 1
```

Batch results are cached under `results/<run>/predictions/*_batches/` and `mutant_predictions/*_batches/`.
If a run is interrupted, re-running the same command skips finished batches automatically.

---

## Library (`lib/`)

| Module | Functions |
|--------|-----------|
| `predict.py` | `run_prediction()` — single FASTA, single protease |
| `batch.py` | `run_prediction_batched()` — splits large FASTA, parallel, resumable |
| `regions.py` | `find_low_cleavage_regions()` — detect resistant stretches |
| `mutate.py` | `load_sequences()`, `create_mutation_fastas()` — generate S,T→P/E FASTAs |
| `compare.py` | `calculate_mutation_sum()`, `find_candidates()` — score and rank positions |

---

## Input data

```
data/
  original_files/          # canonical reference sequences
    tau_only.fasta          # full-length human tau (P10636)
    brain_single.fasta      # small brain-enriched protein subset
    brain_some.fasta        # medium brain-enriched protein subset
    brain_many.fasta        # large brain-enriched protein subset
    UP000005640_9606.fasta  # full human proteome (UniProt)
    *.json                  # protein metadata (tissue category, RNA levels)
  brain_elevated/          # primary target dataset
    tissue_category_rna_brain_Group_Group_enriched.fasta   # 452 proteins
    tissue_category_rna_brain_Tissue_Tissue_enhanced.fasta # 1230 proteins
    tissue_category_rna_brain_Tissue_Tissue_enriched.fasta # 468 proteins
  proteases_list.txt        # default proteases: A01.009 C01.034 C01.060 C01.032 C13.004
  cathepsin_B_list.txt      # C01.060 only
  cathepsin_L_list.txt      # C01.032 only
  cathepsins_combined_list.txt  # C01.032 + C01.060
  website_tau_*.csv         # reference predictions from ProsperousPlus website (tau)
  sample_tau_mutation.csv   # sample mutation analysis output for tau
```

---

## ProsperousPlus tool reference

Tool location: `/shared/macdata/groups/ppc/code/ProsperousPlus/`

```bash
# Direct tool usage (prediction mode)
conda activate prosperousplus
cd /shared/macdata/groups/ppc/code/ProsperousPlus
python Prosperousplus.py \
  --predictfile  /path/to/input.fasta \
  --outputpath   /path/to/output_dir \
  --inputType    fasta \
  --protease     A01.009 \
  --mode         prediction \
  --processNum   10 \
  --PLOT         No
```

Output: `output_dir/results.csv` with columns `protease, sequence_id, position, seqs, prediction, pro`.

Available models: 110+ in `finalModels/` (C01.032 = CathL, C01.060 = CathB, A01.009 = CathD, C13.004 = Legumain, …).

---

## Environment setup

```bash
conda create -n prosperousplus python=3.7
conda activate prosperousplus
pip install -r /shared/macdata/groups/ppc/code/ProsperousPlus/requirements.txt
# Also requires Java (JDK 1.8) and R
```

---

## Archive

Old scripts and result trees are preserved in:
```
/shared/macdata/groups/ppc/projects/ProsperousPlus_older/
  older/            # old standalone scripts (superseded by lib/)
  data_stale/       # pre-generated batch files (regenerated automatically)
  results_stale/    # old result trees from exploratory runs
```
