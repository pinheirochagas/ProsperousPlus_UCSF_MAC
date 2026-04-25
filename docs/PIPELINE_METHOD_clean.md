# ProsperousPlus — Pipeline Method

> Full algorithmic description of the biomarker discovery pipeline.
> For quick-start commands and file paths see `CLAUDE.md`.

---

## Overview

The pipeline starts from a protein FASTA and returns ranked 8-mer peptides that are predicted to be **resistant to lysosomal protease cleavage** — candidate biomarkers in CSF or blood that survive intact long enough to be detected.

| Step | Name | Key output |
|:---:|---|---|
| **Input** | Protein FASTA | — |
| ↓ | | |
| **1** | Cleavage prediction — *CathD + Legumain (sentinel proteases)* | Per-residue cleavage probability (0–1) |
| ↓ | | |
| **2** | Low-cleavage region detection — *stretches resistant to both proteases* | Regions (start, end) per protein |
| ↓ | | |
| **3** | Single-mutant FASTA generation — *one S/T mutated at a time* | `original.fasta` · `single_ST_to_P.fasta` · `single_ST_to_E.fasta` |
| ↓ | | |
| **4** | Cleavage prediction — *CathL + CathB (scoring proteases)* | 4 prediction CSVs (orig/mutant × L/B) |
| ↓ | | |
| **5** | Per-residue contribution scoring — *diff original vs. single mutant* | `per_residue_contributions.csv` · `mutation_scores.csv` · `biomarker_candidates.csv` |
| ↓ | | |
| **6** | Phosphorylation site annotation — *EPSD / dbPTM / UniProt* | `annotated_candidates.csv` · `annotated_summary.csv` |
| ↓ | | |
| **Viz** | Visualisation | `pipeline_summary.png` · `residue_detail.png` |

---

## Step 1 — Sentinel protease prediction

Two proteases scan the entire input proteome **in parallel**:

| Protease | MEROPS ID | Role in biology |
|---|---|---|
| Cathepsin D | A01.009 | Major lysosomal aspartyl protease |
| Legumain | C13.004 | Asparagine endopeptidase (lysosomal) |

ProsperousPlus scores every 8-mer window along each protein sequence and outputs a **cleavage probability** (`pro`, range 0–1) at each position.

**Example — tau fragment:**

| Position | K | Y | V | S | S | V | T | S | R | … |
|---|---|---|---|---|---|---|---|---|---|---|
| CathD score | 0.05 | 0.04 | 0.03 | 0.02 | 0.02 | 0.02 | 0.03 | 0.04 | … | … |
| Legumain score | 0.07 | 0.06 | 0.04 | 0.03 | 0.02 | 0.03 | 0.04 | 0.06 | … | … |

Scores are written to `step1_predictions/cathepsinD_cleavage.csv` and `legumain_cleavage.csv`.

---

## Step 2 — Low-cleavage region detection

A sliding window scans for stretches where **both** sentinel proteases predict low cleavage:

| Parameter | Value | Meaning |
|---|---|---|
| `min_length` | 25 | Minimum consecutive positions in the region |
| `max_exceptions` | 2 | Max positions above threshold across both scores |
| `score_threshold` | 0.3 | Probability above which a position is "cleavable" |

**Example — positions 438–452:**

| Position | 438 | 439 | 440 | 441 | 442 | 443 | 444 | 445 | 446 | 447 | 448 | 449 | 450 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CathD | 0.1 | 0.1 | 0.1 | 0.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 | 0.1 | 0.1 | 0.0 | … |
| Legumain | 0.1 | 0.0 | 0.1 | 0.0 | 0.1 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.1 | 0.1 | … |

The entire span qualifies as a **resistant region** (≤ 2 exceptions). Results are written to `step2_regions/low_cleavage_regions.csv`.

Overlapping candidate regions are deduplicated, keeping the longest. Each region is then **extended by 5 residues on each side**:

```
Detected region:  |───── 39 residues ─────|
Extended region:  |──5──|───── 39 ─────|──5──|   (total ≈ 49 aa)
```

---

## Step 3 — Single-mutant FASTA generation

For every S or T in the extended region, **one mutant sequence** is created where only that single residue is changed.

Two mutation types are generated in parallel (one per scoring protease):

| Mutation | Protease scored | Rationale |
|---|---|---|
| S/T → P | Cathepsin L | Proline is rigid, blocks S/T-dependent cleavage |
| S/T → E | Cathepsin B | Glutamate mimics phosphorylation (charge mimic) |

**Example — extended region `V S S V T S R T G S S …`:**

| S/T site | Original | S/T → P mutant | S/T → E mutant |
|---|---|---|---|
| st1 (S) | V **S** S V T S R … | V **P** S V T S R … | V **E** S V T S R … |
| st5 (T) | V S S V **T** S R … | V S S V **P** S R … | V S S V **E** S R … |
| st7 (S) | … R **S** T G S S … | … R **P** T G S S … | … R **E** T G S S … |

One pair of sequences is written per S/T position, to:

- `step3_fastas/original.fasta`
- `step3_fastas/single_ST_to_P.fasta`
- `step3_fastas/single_ST_to_E.fasta`

---

## Step 4 — Scoring protease prediction

CathL and CathB are run on the original and all single-mutant sequences:

| Input | Protease | Output |
|---|---|---|
| `original.fasta` | CathL | `cathepsinL_original.csv` |
| `single_P.fasta` | CathL | `cathepsinL_single_P.csv` |
| `original.fasta` | CathB | `cathepsinB_original.csv` |
| `single_E.fasta` | CathB | `cathepsinB_single_E.csv` |

CathL and CathB predictions run **in parallel** using half the available cores each.

---

## Step 5 — Per-residue contribution scoring

### What is "contribution"?

For each S or T at position *k*, we compare the original cleavage score to the score of the mutant where **only that residue** is changed, evaluated at the 8-mer window where position *k* acts as the cleavage site.

**Example — position k = S:**

| | Sequence window | CathL score |
|---|---|---|
| Original | V **S**[k] V T | 0.78 |
| S → P mutant | V **P**[k] V T | 0.08 |
| **contrib_L(k)** | 0.78 − 0.08 | **+0.70** ← S is protecting! |

More formally:

```
contrib_L(k)     = CathL_original[window_k] − CathL_singleP_k[window_k]
contrib_B(k)     = CathB_original[window_k] − CathB_singleE_k[window_k]
contrib_total(k) = contrib_L(k) + contrib_B(k)
```

- **Positive contribution** → this S/T protects the peptide from cleavage.
- **Negative contribution** → the wild-type residue was not contributing to resistance.

### Window-level mutation sum

Per-residue contributions are summed across every 8-mer window:

```
mutation_sum = Σ contrib_L(k) + Σ contrib_B(k)
               for all S/T k within the window
```

**Example — 8-mer window at positions 443–450 in tau (`K Y V S S V T S`):**

| Position | 443 | 444 | 446 | 447 |
|---|---|---|---|---|
| Residue | S | S | T | S |
| contrib_L | 0.00 | 0.55 | 0.70 | 0.20 |
| contrib_B | 0.00 | 0.08 | 0.20 | 0.05 |

**mutation_sum = (0.55 + 0.70 + 0.20) + (0.08 + 0.20 + 0.05) = 1.78**

### Threshold

Windows with `mutation_sum ≥ 0.475` are called **biomarker candidates**. The threshold is set to match the pTau217 reference level.

### Output files

| File | Content |
|---|---|
| `step5_scores/per_residue_contributions.csv` | One row per S/T: `region_id, abs_position, aa, contrib_L, contrib_B, contrib_total` |
| `step5_scores/mutation_scores.csv` | One row per 8-mer window across all regions |
| `step5_scores/biomarker_candidates.csv` | Subset of windows above the 0.475 threshold |

---

## Step 6 — Phosphorylation site annotation

### Biological rationale

Steps 1–5 identify S/T residues whose presence **protects a peptide from lysosomal protease cleavage**. Step 6 asks: are those same residues **known to be phosphorylated in humans**?

A residue that is both protease-resistant and a confirmed human phosphosite is a strong biomarker candidate — it survives intact long enough to be detected in CSF or blood, and its phosphorylation status can be measured by standard phosphoproteomics assays.

### Input

The starting point is **`master_biomarker_candidates.csv`** from a completed cohort run. Each row represents one S/T residue that crossed the `contrib_total ≥ 0.475` threshold.

| uniprot_id | protein_name | region | position | aa | contrib_L | contrib_B | contrib_total |
|---|---|---|---|---|---|---|---|
| Q13009 | TIAM1_HUMAN | 1421–1462 | 1459 | S | 0.6604 | 0.3732 | 1.0336 |
| P23471 | PTPRZ_HUMAN | 1360–1396 | 1388 | T | 0.6024 | 0.4056 | 1.008 |
| … | … | … | … | … | … | … | … |

*(916 rows across 264 proteins in the enriched-brain cohort)*

### Databases queried

For each `(uniprot_id, position)` pair, three open databases are queried in order of breadth:

| # | Database | Source file | Coverage | Value |
|---|---|---|---|---|
| 1 | **EPSD 2.0** | `epsd/Homo sapiens.txt` | 1M+ human phosphosites | Broadest bulk human phospho data |
| 2 | **dbPTM** | `dbptm/Phosphorylation` | 1.6M sites (filtered to human) | Catches sites not in EPSD |
| 3 | **UniProt** | `uniprot/human_reviewed_modres.tsv` | Swiss-Prot curated features | Highest confidence, manually reviewed |

All three are loaded from local files (one-time download — see `CLAUDE.md`). The join key is exact `(uniprot_id, position)` — 1-based residue number in the canonical UniProt sequence.

### Method

#### 1. Database loading and parsing

Each database is loaded once into memory as a normalised table with schema: `uniprot_id, position, aa, pubmed_ids, source`.

**EPSD 2.0** — Tab-delimited with a header. `UniProt ID` and `Position` columns used directly. Only rows with a valid integer position are kept. Human-only file; no species filtering needed.

**dbPTM** — Tab-delimited, no header. Columns: `entry_name, uniprot_id, position, ptm_type, pubmed_ids, sequence_context`. Filtered at load time to rows whose `uniprot_id` appears in the human UniProt accession set. The amino acid is inferred from the central residue of the 21-character sequence context.

**UniProt** — REST-streamed TSV of all Swiss-Prot reviewed human proteins with `Modified residue` feature annotations. Each `MOD_RES` token is parsed with a regular expression; only entries whose `/note` contains `phospho` (case-insensitive) are kept. Example entry:

```
MOD_RES 519; /note="Phosphoserine; by CK1 and PDPK1"; /evidence="ECO:0000269|PubMed:14761950;..."
```

#### 2. Annotation join

An **inner join** is performed against the candidates on the composite key `(uniprot_id, position)` for each database. Only exact matches are retained.

```
candidates  ⟕ EPSD    → hits labelled source = "EPSD"
            ⟕ dbPTM   → hits labelled source = "dbPTM"
            ⟕ UniProt → hits labelled source = "UniProt"

All three hit tables concatenated → annotated_candidates.csv  (long format)
```

#### 3. Summary collapse

Long-format hits are collapsed to one row per candidate residue:

| Column | Description |
|---|---|
| `n_sources` | Count of distinct databases where the residue was found (0–3) |
| `sources` | Semicolon-separated sorted list of database names |
| `n_pubmed` | Count of unique PubMed IDs across all sources |
| `pubmed_ids` | Semicolon-separated sorted list of PMIDs (de-duplicated) |

### Prioritization

Candidates are ranked by:

1. **Primary:** `n_sources` descending (3 > 2 > 1 > 0)
2. **Secondary:** `contrib_total` descending (stronger protease resistance first)

A residue supported by all three databases at a position that strongly confers protease resistance is the most actionable biomarker candidate.

### Output files

Written to `step6_ptm/` inside the run or cohort folder:

| File | Content |
|---|---|
| `annotated_candidates.csv` | Long format — one row per (residue × database) hit. Includes `ptm_source`, `ptm_pubmed_ids`, `ptm_aa`. |
| `annotated_summary.csv` | One row per residue, sorted by `n_sources` then `contrib_total`. |

Residues with no phosphorylation evidence in any database are absent from both outputs but remain in `master_biomarker_candidates.csv`.

### Command

```bash
# Recommended: point to summary.csv — auto-resolves to master_biomarker_candidates.csv
python annotate_ptm.py \
    --candidates results/cohort_.../summary.csv

# Or point directly to the master file
python annotate_ptm.py \
    --candidates results/cohort_.../master_biomarker_candidates.csv
```

---

## Visualisation

Two figures are produced by `plot_results.py`:

### `pipeline_summary.png` — genome browser view

Four horizontal tracks sharing a common protein-position x-axis:

| Track | Signal | Colour |
|---|---|---|
| 1 | CathD cleavage score | Blue fill |
| 2 | Legumain score (shaded bands = low-cleavage regions from Step 2) | Purple fill |
| 3 | CathL + CathB original scores | Green / Orange |
| 4 | mutation_sum per window (top candidate labelled with peptide + position) | Red bars |

### `residue_detail.png` — per-residue contribution grid

Top-10 candidates in a 2 × 5 grid. Each panel shows:

- **Green bars** — CathL contribution (S/T → P effect)
- **Orange bars** — CathB contribution (S/T → E effect)
- **Bold coloured x-labels** — S or T residues

---

## Parameters reference

### Steps 1–5 (`run_pipeline.py` / `run_cohort.py`)

| Parameter | Default | Flag | Meaning |
|---|---|---|---|
| `--fasta` | tau | required | Input protein FASTA |
| `--cores` | 20 | `--cores N` | Total CPU cores; split across jobs |
| `min_length` | 25 | hardcoded | Minimum resistant region length |
| `max_exceptions` | 2 | hardcoded | Max cleavable positions per region |
| `score_threshold` | 0.3 | hardcoded | Cleavage probability cutoff |
| `extension` | 5 | hardcoded | Flanking residues added per side |
| `mutation_threshold` | 0.475 | hardcoded | Min contrib_total for candidacy |
| Sentinel proteases | CathD + Lgmn | hardcoded | Used in Step 1 region detection |
| Scoring proteases | CathL + CathB | hardcoded | Used in Steps 4–5 mutation scoring |

### Step 6 (`annotate_ptm.py`)

| Parameter | Default | Meaning |
|---|---|---|
| `--candidates` | required | Path to `master_biomarker_candidates.csv`, `summary.csv` (auto-resolved), or single-run `biomarker_candidates.csv` |
| `--epsd` | `data/reference/ptm_sources/epsd/Homo sapiens.txt` | EPSD 2.0 human phospho sites |
| `--dbptm` | `data/reference/ptm_sources/dbptm/Phosphorylation` | dbPTM experimental sites |
| `--uniprot` | `data/reference/ptm_sources/uniprot/human_reviewed_modres.tsv` | UniProt curated modified residues |
| `--out` | `<run_root>/step6_ptm/` (inferred) | Output directory for annotation CSVs |
