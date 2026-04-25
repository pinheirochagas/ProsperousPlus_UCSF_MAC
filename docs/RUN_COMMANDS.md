# Run Commands

## Setup

```bash
cd /shared/macdata/groups/ppc/projects/ProsperousPlus
conda activate prosperousplus
```

---

## Single protein (tau)

```bash
python run_pipeline.py
```

---

## Cohort — brain-enriched tissue (468 proteins)

The cohort runner uses a two-phase architecture:

- **Phase 1** — CathD + Legumain on all proteins in batches of 50 (20 model loads total, both proteases in parallel). Results are split per-protein into each protein's `step1_predictions/` folder.
- **Phase 2** — Steps 2–5 per protein, run `--parallel` at a time. Step 1 is skipped because the CSVs already exist.

---

**Option A — Run both phases in one command** (recommended for short runs or overnight):

```bash
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \
  --parallel 5 --cores 20
```

---

**Option B — Run Phase 1 first, then Phase 2 separately** (useful when you want to verify Step 1 before proceeding, or split across sessions):

```bash
# Step 1: Phase 1 only — CathD + Legumain on all 468 proteins
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \
  --cores 20 --step1-only
```

The command prints the exact output directory name (e.g. `cohort_tissue_category_rna_brain_Tissue_Tissue_enriched_2026-04-24_1030`). Use that name to continue:

```bash
# Step 2: Phase 2 — Steps 2-5 per protein (pass the same output dir)
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \
  --output cohort_tissue_category_rna_brain_Tissue_Tissue_enriched_2026-04-24_1030 \
  --parallel 5 --cores 20
```

Phase 1 will detect all step1 CSVs already exist and skip straight to Phase 2.

---

**Test run** (first 5 proteins only — benchmark timing before committing to full run):

```bash
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \
  --parallel 5 --cores 20 --limit 5
```

---

**Resume** interrupted run (both phases detect completed work and skip it automatically):

```bash
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \
  --output cohort_tissue_category_rna_brain_Tissue_Tissue_enriched_2026-04-24_1030 \
  --parallel 5 --cores 20
```

---

**Merge only** (regenerate master files from completed proteins without re-running anything):

```bash
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enriched.fasta \
  --output cohort_tissue_category_rna_brain_Tissue_Tissue_enriched_2026-04-24_1030 \
  --parallel 5 --cores 20 --merge-only
```

---

## Other brain-enriched FASTAs

```bash
# tissue enhanced (1230 proteins)
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Tissue_Tissue_enhanced.fasta \
  --parallel 5 --cores 20

# group enriched (452 proteins)
python run_cohort.py \
  --fasta data/brain_elevated/tissue_category_rna_brain_Group_Group_enriched.fasta \
  --parallel 5 --cores 20
```

---

## Visualise results (single run)

```bash
python plot_results.py --run-dir results/tau_only_2026-04-23_2204
```

---

## Output locations

| Run | Directory |
|-----|-----------|
| Tau (single) | `results/tau_only_YYYY-MM-DD_HHMM/` |
| Cohort | `results/cohort_<fasta_stem>_YYYY-MM-DD_HHMM/` |
| Phase 1 batch cache | `results/cohort_.../_step1_batches/` · `_step1_cathD/` · `_step1_legumain/` |
| Per-protein step1 | `results/cohort_.../sp\|XXXXX\|PROTEIN/step1_predictions/` |
| Master candidates | `results/cohort_.../master_biomarker_candidates.csv` |
| Per-protein summary | `results/cohort_.../summary.csv` |
