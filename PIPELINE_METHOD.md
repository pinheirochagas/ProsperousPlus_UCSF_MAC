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

| Parameter          | Default | Flag            | Meaning                             |
|--------------------|---------|-----------------|-------------------------------------|
| `--fasta`          | tau     | required        | Input protein FASTA                 |
| `--cores`          | 20      | `--cores N`     | Total CPU cores; split across jobs  |
| `min_length`       | 25      | hardcoded       | Minimum resistant region length     |
| `max_exceptions`   | 2       | hardcoded       | Max cleavable positions per region  |
| `score_threshold`  | 0.3     | hardcoded       | Cleavage probability cutoff         |
| `extension`        | 5       | hardcoded       | Flanking residues added per side    |
| `mutation_threshold`| 0.475  | hardcoded       | Min mutation sum for candidacy      |
| Sentinel proteases | CathD + Lgmn | hardcoded  | Used in Step 1 region detection     |
| Scoring proteases  | CathL + CathB | hardcoded | Used in Steps 4–5 mutation scoring  |
