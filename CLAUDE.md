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
python scripts/run_pipeline.py

# Run on brain-enriched proteins (large dataset, batched)
python scripts/run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta
python scripts/run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enhanced.fasta
python scripts/run_pipeline.py --fasta data/brain_elevated/tissue_category_rna_brain_Group_Group_enriched.fasta

# Resume an interrupted run (just re-run the same command — completed steps are skipped)
python scripts/run_pipeline.py --fasta data/brain_elevated/...same file...
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

**Per-residue contribution (Step 4–5):**

For every S or T at position *k* in an extended region, a **single-position mutant** is generated (only that one residue changed, everything else unchanged):

```
contrib_L(k) = CathL_original[window_k] − CathL_singleP_k[window_k]   # S/T→P
contrib_B(k) = CathB_original[window_k] − CathB_singleE_k[window_k]   # S/T→E
contrib_total(k) = contrib_L(k) + contrib_B(k)
```

`window_k` is the 8-mer window whose cleavage site falls on position *k*.
A **positive** contribution means that S/T is protecting the peptide from cleavage — mutating it away *reduces* the predicted cleavage score.

**Window-level mutation sum (Step 5 output):**

```
mutation_sum(window) = Σ contrib_L(k) + Σ contrib_B(k)   for all S/T k within the 8-mer
```

Windows with `mutation_sum ≥ 0.475` (= pTau217 threshold) are called **biomarker candidates**.

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
  step6_ptm/                            ← written by scripts/annotate_ptm.py (optional)
    annotated_candidates.csv             # long format: one row per (candidate × source) hit
    annotated_summary.csv               # one row per residue, ranked by n_sources then contrib_total
```

---

## PTM annotation (Step 6 — optional, standalone)

Annotates Step 5 biomarker candidates with known human phosphorylation evidence
from three open databases (EPSD 2.0, dbPTM, UniProt).

### One-time data download

```bash
# EPSD 2.0 — human phosphorylation sites
mkdir -p data/reference/ptm_sources/epsd
curl -L -o data/reference/ptm_sources/epsd/Homo_sapiens.zip \
  "http://epsd.biocuckoo.cn/Download/Homo%20sapiens.zip"
unzip -d data/reference/ptm_sources/epsd/ \
  data/reference/ptm_sources/epsd/Homo_sapiens.zip

# dbPTM — experimental phosphorylation sites (all species; filtered to human at runtime)
# Visit https://biomics.lab.nycu.edu.tw/dbPTM/download.php
# → "Experimental & Putative PTM Sites" → Phosphorylation → download MAC/Linux .tgz
mkdir -p data/reference/ptm_sources/dbptm
tar -xzf ~/Downloads/Phosphorylation.tgz -C data/reference/ptm_sources/dbptm/

# UniProt — Swiss-Prot curated modified-residue annotations (human, reviewed)
mkdir -p data/reference/ptm_sources/uniprot
curl -L --compressed \
  -o data/reference/ptm_sources/uniprot/human_reviewed_modres.tsv \
  "https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&format=tsv&fields=accession,id,ft_mod_res"
```

### Run annotation

```bash
conda activate prosperousplus
cd /shared/macdata/groups/ppc/projects/ProsperousPlus

# Annotate tau candidates (output written to results/<run>/step6_ptm/)
python scripts/annotate_ptm.py \
  --candidates results/tau_only_2026-04-23_1430/step5_scores/biomarker_candidates.csv

# Annotate ALL scored S/T residues (not just above-threshold candidates)
python scripts/annotate_ptm.py \
  --candidates results/tau_only_2026-04-23_1430/step5_scores/per_residue_contributions.csv

# Custom source file paths (defaults shown above are used when flags are omitted)
python scripts/annotate_ptm.py \
  --candidates results/<run>/step5_scores/biomarker_candidates.csv \
  --epsd    "data/reference/ptm_sources/epsd/Homo sapiens.txt" \
  --dbptm   data/reference/ptm_sources/dbptm/Phosphorylation \
  --uniprot data/reference/ptm_sources/uniprot/human_reviewed_modres.tsv \
  --out     results/<run>/step6_ptm/
```

### Output columns

`annotated_summary.csv` (one row per residue, ranked by evidence strength):

| Column | Meaning |
|--------|---------|
| `uniprot_id` | UniProt accession |
| `position` | 1-based residue position in full protein |
| `aa` | Amino acid (S or T) |
| `contrib_total` | Pipeline cleavage-protection score |
| `n_sources` | Number of databases confirming this site (max 3) |
| `sources` | Semicolon-separated database names (EPSD, dbPTM, UniProt) |
| `n_pubmed` | Count of unique PubMed IDs across all confirming sources |
| `pubmed_ids` | Semicolon-separated PMIDs (look up at pubmed.ncbi.nlm.nih.gov) |

---

## Batching (for large inputs)

For inputs with many sequences, predictions are split into batches and run in parallel.

```bash
# Default: 50 sequences/batch, cpu_count-1 workers
python scripts/run_pipeline.py --fasta data/brain_elevated/...fasta

# Custom batch size and workers
python scripts/run_pipeline.py --fasta data/brain_elevated/...fasta --batch-size 25 --workers 8

# Disable batching (single-sequence files like tau)
python scripts/run_pipeline.py --fasta data/original_files/tau_only.fasta --batch-size 1
```

Batch results are cached under `results/<run>/predictions/*_batches/` and `mutant_predictions/*_batches/`.
If a run is interrupted, re-running the same command skips finished batches automatically.

---

## Library (`lib/`)

| Module | Key functions |
|--------|---------------|
| `predict.py` | `run_prediction()` — single FASTA, single protease |
| `batch.py` | `run_prediction_batched()` — splits large FASTA, parallel, resumable |
| `regions.py` | `find_low_cleavage_regions()` — detect resistant stretches |
| `mutate.py` | `load_sequences()`, `create_single_mutant_fastas()` — one mutant per S/T position |
| `compare.py` | `calculate_single_residue_mutation_sum()`, `find_candidates()` — per-residue contributions & window sums |
| `ptm_annotate.py` | `load_epsd()`, `load_dbptm()`, `load_uniprot_modres()`, `annotate_candidates()`, `summarize_annotations()` — PTM annotation (Step 6) |

---

## Input data

```
data/
  reference/                   # PTM annotation reference databases (one-time download)
    ptm_sources/
      epsd/
        Homo sapiens.txt        # EPSD 2.0 human phospho sites (tab-delimited, has header)
        Homo sapiens.fasta      # EPSD protein FASTA sequences
      dbptm/
        Phosphorylation         # dbPTM experimental phosphorylation (all species, no header)
      uniprot/
        human_reviewed_modres.tsv  # UniProt Swiss-Prot human modified-residue features
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
