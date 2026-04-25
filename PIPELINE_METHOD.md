# ProsperousPlus — Pipeline Method

> Full algorithmic description of the biomarker discovery pipeline.
> For quick-start commands and file paths see `CLAUDE.md`.

---

## Overview

The pipeline starts from a protein FASTA and returns ranked 8-mer peptides that
are predicted to be **resistant to lysosomal protease cleavage** — candidate
biomarkers in CSF or blood that survive intact long enough to be detected.

```
Input FASTA
    │
    ▼
┌─────────────────────────────────────────────┐
│  Step 1  Cleavage prediction (CathD + Lgmn) │  ← "sentinel" proteases
└─────────────────────────────────────────────┘
    │  per-residue cleavage probability (0–1)
    ▼
┌─────────────────────────────────────────────┐
│  Step 2  Low-cleavage region detection      │  stretch resistant to both
└─────────────────────────────────────────────┘
    │  regions (start, end) per protein
    ▼
┌─────────────────────────────────────────────┐
│  Step 3  Single-mutant FASTA generation     │  one S/T mutated at a time
└─────────────────────────────────────────────┘
    │  original.fasta  +  single_ST_to_P.fasta  +  single_ST_to_E.fasta
    ▼
┌─────────────────────────────────────────────┐
│  Step 4  Cleavage prediction (CathL + CathB)│  ← "scoring" proteases
└─────────────────────────────────────────────┘
    │  4 prediction CSVs (orig/mutant × L/B)
    ▼
┌─────────────────────────────────────────────┐
│  Step 5  Per-residue contribution scoring   │  diff original vs. single mutant
└─────────────────────────────────────────────┘
    │  per_residue_contributions.csv
    │  mutation_scores.csv (window sums)
    │  biomarker_candidates.csv (≥ 0.475)
    │     [cohort: master_biomarker_candidates.csv]
    ▼
┌─────────────────────────────────────────────┐
│  Step 6  Phosphorylation site annotation    │  open databases (EPSD/dbPTM/UniProt)
└─────────────────────────────────────────────┘
    │  annotated_candidates.csv (all hits)
    │  annotated_summary.csv   (ranked by evidence strength)
    ▼
┌─────────────────────────────────────────────┐
│  Visualisation  pipeline_summary.png        │
│                 residue_detail.png          │
└─────────────────────────────────────────────┘
```

---

## Step 1 — Sentinel protease prediction

Two proteases are used to scan the entire input proteome in **parallel**:

| Protease    | MEROPS ID | Role in biology                        |
|-------------|-----------|----------------------------------------|
| Cathepsin D | A01.009   | Major lysosomal aspartyl protease      |
| Legumain    | C13.004   | Asparagine endopeptidase (lysosomal)   |

ProsperousPlus scores every 8-mer window along each protein sequence and outputs
a **cleavage probability** (`pro`, range 0–1) at each position.

```
Protein sequence (example — tau fragment):
  ...K Y V S S V T S R T G S S K T N I...
     ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑ ...   (each = one 8-mer window)

CathD score:  0.05 0.04 0.03 0.02 0.02 0.02 0.03 0.04 ...
Legumain:     0.07 0.06 0.04 0.03 0.02 0.03 0.04 0.06 ...
```

Scores are written to `step1_predictions/cathepsinD_cleavage.csv` and
`legumain_cleavage.csv`.

---

## Step 2 — Low-cleavage region detection

A sliding window scans for stretches where **both** sentinel proteases predict
low cleavage, using these criteria:

| Parameter         | Value | Meaning                                          |
|-------------------|-------|--------------------------------------------------|
| `min_length`      | 25    | Minimum consecutive positions in the region      |
| `max_exceptions`  | 2     | Max positions above threshold across both scores |
| `score_threshold` | 0.3   | Probability above which a position is "cleavable"|

```
Position:  ...438 439 440 441 442 443 444 445 446 447 448 449 450 451 452...
CathD:        0.1  0.1  0.1  0.1  0.0  0.0  0.0  0.0  0.1  0.1  0.1  0.0 ...
Legumain:     0.1  0.0  0.1  0.0  0.1  0.0  0.0  0.1  0.0  0.0  0.1  0.1 ...
                                                                            
                ←──────────── RESISTANT REGION (≤2 exceptions) ──────────→
                    written to step2_regions/low_cleavage_regions.csv
```

Overlapping candidate regions are deduplicated, keeping the longest.  
Each region is then **extended by 5 residues** on each side to provide flanking
context for the mutation analysis.

```
Detected region:  |───── 39 residues ─────|
Extended region:  |──5──|───── 39 ─────|──5──|   (total ≈ 49 aa)
```

---

## Step 3 — Single-mutant FASTA generation

For every S or T in the extended region, **one mutant sequence** is created
where only that single residue is changed — all other positions are unchanged.

Two mutation types are generated in parallel (one per scoring protease):

| Mutation | Protease scored | Rationale                                      |
|----------|-----------------|------------------------------------------------|
| S/T → P  | Cathepsin L     | Proline is rigid, blocks S/T-dependent cleavage|
| S/T → E  | Cathepsin B     | Glutamate mimics phosphorylation (charge mimic)|

```
Extended region:  V  S  S  V  T  S  R  T  G  S  S  ...
                     ↑           ↑     ↑     ↑  ↑
                    st1          st5   st7   st9 st10   (1-based within region)

For st1  (S → P):  V [P] S  V  T  S  R  T  G  S  S  ...  id = region_id_st1
For st1  (S → E):  V [E] S  V  T  S  R  T  G  S  S  ...  id = region_id_st1

For st5  (T → P):  V  S  S  V [P] S  R  T  G  S  S  ...  id = region_id_st5
For st5  (T → E):  V  S  S  V [E] S  R  T  G  S  S  ...  id = region_id_st5

...one pair of sequences per S/T position
```

Sequences are written to:
- `step3_fastas/original.fasta`
- `step3_fastas/single_ST_to_P.fasta`
- `step3_fastas/single_ST_to_E.fasta`

---

## Step 4 — Scoring protease prediction

CathL and CathB are run on the original and all single-mutant sequences:

```
                ┌────────────────┐
original.fasta  │  CathL predict │ → cathepsinL_original.csv
                └────────────────┘
                ┌────────────────┐
single_P.fasta  │  CathL predict │ → cathepsinL_single_P.csv   (S/T→P mutants)
                └────────────────┘

                ┌────────────────┐
original.fasta  │  CathB predict │ → cathepsinB_original.csv
                └────────────────┘
                ┌────────────────┐
single_E.fasta  │  CathB predict │ → cathepsinB_single_E.csv   (S/T→E mutants)
                └────────────────┘
```

CathL and CathB predictions run **in parallel** using half the available cores each.

---

## Step 5 — Per-residue contribution scoring

### What is "contribution"?

For each S or T at position *k* in the extended region, we compare the original
cleavage score to the score of the mutant where **only that residue** is changed,
at the specific 8-mer window (`window_k`) where position *k* acts as the
cleavage site.

```
                     window_k (8 aa, position k = cleavage site)
                     ┌─────────────┐
original sequence:   │  V S[k] V T │  CathL score = 0.78
single-P mutant:     │  V P[k] V T │  CathL score = 0.08
                     └─────────────┘
                              ↓
              contrib_L(k) = 0.78 − 0.08 = +0.70   ← large positive → S is protecting!
```

More formally:

```
contrib_L(k) = CathL_original[window_k] − CathL_singleP_k[window_k]
contrib_B(k) = CathB_original[window_k] − CathB_singleE_k[window_k]
contrib_total(k) = contrib_L(k) + contrib_B(k)
```

**Positive contribution** → this S/T is protecting the peptide from cleavage.  
**Negative contribution** → mutating this residue actually *increases* predicted
cleavage (the wild-type residue was not contributing to resistance).

### Window-level mutation sum

All per-residue contributions are aggregated to a **mutation sum** for each
sliding 8-mer window:

```
Window [i+1 … i+8]:

  mutation_sum = Σ contrib_L(k) + Σ contrib_B(k)
                   for all S/T k that fall within this 8-mer window
```

```
8-mer window (positions 443–450 in tau):

  K  Y  V [S] [S]  V [T] [S]
            443  444    446  447

  contrib_L:   0.00  0.55   0.70  0.20
  contrib_B:   0.00  0.08   0.20  0.05
                ─────────────────────────
  mutation_sum = (0.55+0.70+0.20) + (0.08+0.20+0.05) = 1.78
```

### Threshold

Windows with `mutation_sum ≥ 0.475` are called **biomarker candidates**.
The threshold is set to match the pTau217 reference level.

### Output files

| File | Content |
|------|---------|
| `step5_scores/per_residue_contributions.csv` | One row per S/T: `region_id, abs_position, aa, contrib_L, contrib_B, contrib_total` |
| `step5_scores/mutation_scores.csv` | One row per 8-mer window across all regions |
| `step5_scores/biomarker_candidates.csv` | Subset of windows above the 0.475 threshold |

---

## Step 6 — Phosphorylation site annotation

### Biological rationale

Steps 1–5 identify specific S/T residues whose presence **protects a peptide
from lysosomal protease cleavage**.  Step 6 asks: are those same residues
**known to be phosphorylated in humans**?

A residue that is both protease-resistant and a confirmed human phosphosite
is a strong biomarker candidate — it survives intact long enough to be
detected in CSF or blood, and its phosphorylation status can be measured by
standard phosphoproteomics assays.

### Input

The starting point is **`master_biomarker_candidates.csv`** from a completed
cohort run (or the equivalent `biomarker_candidates.csv` from a single-protein
run).  Each row represents one S/T residue that individually crossed the
`contrib_total ≥ 0.475` threshold when mutated.

```
uniprot_id | protein_name | region | sequence_id | position | aa | contrib_L | contrib_B | contrib_total
Q13009     | TIAM1_HUMAN  | 1421-1462 | sp|Q13009|... | 1459 | S | 0.6604 | 0.3732 | 1.0336
P23471     | PTPRZ_HUMAN  | 1360-1396 | sp|P23471|... | 1388 | T | 0.6024 | 0.4056 | 1.008
...        (916 rows across 264 proteins in the enriched-brain cohort)
```

### Databases queried

For each `(uniprot_id, position)` pair, the script queries three open
databases in order of breadth:

| # | Database | Source file | Coverage | Value |
|---|----------|-------------|----------|-------|
| 1 | **EPSD 2.0** | `epsd/Homo sapiens.txt` | 1M+ human phosphosites | Broadest bulk human phospho data |
| 2 | **dbPTM** | `dbptm/Phosphorylation` | 1.6M sites all species → filtered to human | Catches sites not in EPSD |
| 3 | **UniProt** | `uniprot/human_reviewed_modres.tsv` | Swiss-Prot curated features | Highest confidence, manually reviewed |

All three are loaded from local files (one-time download — see `CLAUDE.md`).
The join key is exact `(uniprot_id, position)` — 1-based residue number in
the canonical UniProt sequence.

### Method

#### 1. Database loading and parsing

Each database is loaded once into memory as a normalised table with a common
schema: `uniprot_id, position, aa, pubmed_ids, source`.

**EPSD 2.0** (`Homo sapiens.txt`) is tab-delimited with a header row.  The
`UniProt ID` and `Position` columns are used directly.  Only rows with a
valid integer position are kept.  The file covers human only — no species
filtering is needed.

**dbPTM** (`Phosphorylation`) is tab-delimited with no header.  The six
columns are: `entry_name, uniprot_id, position, ptm_type, pubmed_ids,
sequence_context`.  Because the file covers all eukaryotic species, it is
filtered at load time to rows whose `uniprot_id` appears in the set of human
UniProt accessions derived from the UniProt TSV (`Entry` column).  The amino
acid at the modified position is inferred from column 11 of the 21-character
sequence context (the central residue of the ±10 window).

**UniProt** (`human_reviewed_modres.tsv`) is a REST-streamed TSV of all
Swiss-Prot reviewed human proteins with their `Modified residue` feature
annotations.  The `Modified residue` field is free text and may contain
multiple entries per protein, e.g.:

```
MOD_RES 519; /note="Phosphoserine; by CK1 and PDPK1"; /evidence="ECO:0000269|PubMed:14761950;..."
```

Each `MOD_RES` token is parsed with a regular expression.  Only entries
whose `/note` contains the word `phospho` (case-insensitive) are kept.
The residue type (S, T, or Y) is inferred from the note text
(`Phosphoserine`, `Phosphothreonine`, `Phosphotyrosine`).  PubMed IDs are
extracted from the `/evidence` field.

#### 2. Annotation join

For each of the three normalised site tables, an **inner join** is performed
against the candidates DataFrame on the composite key `(uniprot_id, position)`.
An inner join means only exact matches are retained — a candidate that has
no entry at that exact position in a given database produces no row for that
source.

```
candidates (uniprot_id, position, ...)
    ⟕ EPSD   (uniprot_id, position)  → hits labelled source = "EPSD"
    ⟕ dbPTM  (uniprot_id, position)  → hits labelled source = "dbPTM"
    ⟕ UniProt(uniprot_id, position)  → hits labelled source = "UniProt"

All three hit tables are concatenated → annotated_candidates.csv
```

The result (`annotated_candidates.csv`) is in long format: one row per
`(candidate residue × database)` hit.

#### 3. Summary collapse

The long-format hits are collapsed to one row per candidate residue.  For
each residue:

- `n_sources` — count of distinct databases where the residue was found (0–3)
- `sources` — semicolon-separated sorted list of database names
- `n_pubmed` — count of unique PubMed IDs across all sources combined
- `pubmed_ids` — semicolon-separated sorted list of those PMIDs

PubMed IDs are de-duplicated across databases (the same paper may appear in
EPSD and dbPTM; it is counted once).

### Prioritization

```
For each candidate residue:

  n_sources = count of databases where (uniprot_id, position) appears

  Primary rank  : n_sources   (3 > 2 > 1 > 0)
  Secondary rank: contrib_total descending  (stronger protease resistance first)
```

A residue supported by all three databases at a position that strongly
confers protease resistance is the most actionable biomarker candidate.
Residues with `n_sources = 0` (no phosphorylation evidence at that exact
position) do not appear in the output but remain in
`master_biomarker_candidates.csv`.

### Output files

Written to `step6_ptm/` inside the run or cohort folder:

| File | Content |
|------|---------|
| `annotated_candidates.csv` | Long format — one row per (residue × database) hit. Includes `ptm_source`, `ptm_pubmed_ids`, `ptm_aa`. |
| `annotated_summary.csv` | One row per residue. Columns: `uniprot_id`, `position`, `aa`, `contrib_total`, `n_sources`, `sources`, `n_pubmed`, `pubmed_ids`. Sorted by `n_sources` then `contrib_total`. |

Residues with no phosphorylation evidence in any database are absent from
both outputs — they remain in `master_biomarker_candidates.csv` for reference.

### Command

```bash
# Recommended: point to summary.csv — auto-resolves to master_biomarker_candidates.csv
python annotate_ptm.py \
    --candidates results/cohort_tissue_category_rna_brain_Tissue_Tissue_enriched_2026-04-24_0525/summary.csv

# Or point directly to master file
python annotate_ptm.py \
    --candidates results/cohort_.../master_biomarker_candidates.csv
```

---

## Visualisation

Two figures are produced by `plot_results.py`:

### `pipeline_summary.png` — genome browser view

Four horizontal tracks sharing a common protein-position x-axis:

```
Track 1: CathD cleavage score   (blue fill)
Track 2: Legumain score         (purple fill)
         ║ shaded bands = low-cleavage regions detected in Step 2
Track 3: CathL + CathB original scores (green / orange)
Track 4: mutation_sum per window (red bars)
         ▶ top candidate labelled with peptide + position
```

### `residue_detail.png` — per-residue contribution grid

Top-10 candidates arranged in a 2 × 5 grid (left→right, top→bottom):

```
┌──────────────────────┐  ┌──────────────────────┐
│  #1  VSSVTSRT  1.58  │  │  #2  VTSRTGSS  1.47  │
│  ████                │  │  ████                 │
│  ▒▒                  │  │         ▒▒            │
│  S   S   V   T ...   │  │  T   S   R   T ...    │
│  443 444 445 446     │  │  446 447 448 449       │
└──────────────────────┘  └──────────────────────┘
   Green = CathL contribution (S/T→P effect)
   Orange = CathB contribution (S/T→E effect)
   Bold coloured x-labels = S or T residues
```

---

## Parameters reference

### Steps 1–5 (`run_pipeline.py` / `run_cohort.py`)

| Parameter          | Default | Flag            | Meaning                             |
|--------------------|---------|-----------------|-------------------------------------|
| `--fasta`          | tau     | required        | Input protein FASTA                 |
| `--cores`          | 20      | `--cores N`     | Total CPU cores; split across jobs  |
| `min_length`       | 25      | hardcoded       | Minimum resistant region length     |
| `max_exceptions`   | 2       | hardcoded       | Max cleavable positions per region  |
| `score_threshold`  | 0.3     | hardcoded       | Cleavage probability cutoff         |
| `extension`        | 5       | hardcoded       | Flanking residues added per side    |
| `mutation_threshold`| 0.475  | hardcoded       | Min contrib_total for candidacy     |
| Sentinel proteases | CathD + Lgmn | hardcoded  | Used in Step 1 region detection     |
| Scoring proteases  | CathL + CathB | hardcoded | Used in Steps 4–5 mutation scoring  |

### Step 6 (`annotate_ptm.py`)

| Parameter    | Default | Meaning |
|--------------|---------|---------|
| `--candidates` | required | Path to `master_biomarker_candidates.csv`, `summary.csv` (auto-resolved), or single-run `biomarker_candidates.csv` |
| `--epsd`     | `data/reference/ptm_sources/epsd/Homo sapiens.txt` | EPSD 2.0 human phospho sites |
| `--dbptm`    | `data/reference/ptm_sources/dbptm/Phosphorylation` | dbPTM experimental sites |
| `--uniprot`  | `data/reference/ptm_sources/uniprot/human_reviewed_modres.tsv` | UniProt curated modified residues |
| `--out`      | `<run_root>/step6_ptm/` (inferred) | Output directory for annotation CSVs |
