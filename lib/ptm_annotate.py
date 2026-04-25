"""Phosphorylation site annotation from EPSD, dbPTM, and UniProt.

This module loads bulk phosphorylation site data from three open databases
and annotates candidate protease-resistant residues (from Step 5 of the
ProsperousPlus pipeline) with known human phosphorylation evidence.

Sources
-------
EPSD 2.0 : Eukaryotic Phosphorylation Site Database
    http://epsd.biocuckoo.cn/Download.php  →  "Homo sapiens.zip"
    File: data/reference/ptm_sources/epsd/Homo sapiens.txt
    Columns (tab-delimited, has header):
        EPSD ID, UniProt ID, AA, Position, Source, Reference

dbPTM    : post-translational modification database
    https://biomics.lab.nycu.edu.tw/dbPTM/download.php  →  Phosphorylation
    File: data/reference/ptm_sources/dbptm/Phosphorylation
    Columns (tab-delimited, no header):
        entry_name, uniprot_id, position, ptm_type, pubmed_ids, context
    Note: covers all species — filtered to human using UniProt accessions.

UniProt  : Swiss-Prot curated modified-residue annotations
    REST stream: organism_id:9606, reviewed:true, fields: accession,id,ft_mod_res
    File: data/reference/ptm_sources/uniprot/human_reviewed_modres.tsv
    Columns (tab-delimited, has header):
        Entry, Entry Name, Modified residue
    Note: Modified residue field is free text; phospho sites are parsed
    from entries of the form  MOD_RES N; /note="Phospho..."; /evidence="..."
"""

import re
import warnings
from pathlib import Path
from typing import Optional, Set

import pandas as pd


# ── sequence_id helpers (mirrors lib/merge.py) ────────────────────────────────
_SEQ_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|(\S+?)_(\d+-\d+)$")


def _parse_sequence_id(seq_id: str) -> tuple:
    """Return (uniprot_id, protein_name, region) from a pipeline sequence_id."""
    m = _SEQ_RE.match(seq_id)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return seq_id, seq_id, ""


# ── UniProt MOD_RES field parsers ─────────────────────────────────────────────
_MODRES_RE = re.compile(
    r'MOD_RES (\d+); /note="([^"]+)"(?:; /evidence="([^"]*)")?'
)
_PUBMED_RE = re.compile(r"PubMed:(\d+)")


def _parse_modres_field(accession: str, modres_str: str) -> list:
    """Parse a UniProt MOD_RES feature string, keeping only phospho sites.

    Parameters
    ----------
    accession : str
        UniProt accession (Entry column).
    modres_str : str
        Full content of the 'Modified residue' TSV column for one protein.

    Returns
    -------
    list of dict with keys: uniprot_id, position, aa, pubmed_ids, source
    """
    rows = []
    for m in _MODRES_RE.finditer(modres_str):
        note = m.group(2)
        if "phospho" not in note.lower():
            continue
        position = int(m.group(1))
        evidence = m.group(3) or ""

        note_lower = note.lower()
        if "serine" in note_lower:
            aa = "S"
        elif "threonine" in note_lower:
            aa = "T"
        elif "tyrosine" in note_lower:
            aa = "Y"
        else:
            aa = ""

        pubmed_ids = ";".join(sorted(set(_PUBMED_RE.findall(evidence))))
        rows.append({
            "uniprot_id": accession,
            "position":   position,
            "aa":         aa,
            "pubmed_ids": pubmed_ids,
            "source":     "UniProt",
        })
    return rows


# ── Loaders ───────────────────────────────────────────────────────────────────

def load_epsd(path) -> pd.DataFrame:
    """Load EPSD 2.0 human phosphorylation sites.

    Parameters
    ----------
    path : str or Path
        Path to 'Homo sapiens.txt' (tab-delimited, has header).

    Returns
    -------
    DataFrame with columns: uniprot_id, position, aa, pubmed_ids, source
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df = df.rename(columns={
        "UniProt ID": "uniprot_id",
        "AA":         "aa",
        "Position":   "position",
        "Reference":  "pubmed_ids",
    })
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df = df.dropna(subset=["position"])
    df["position"] = df["position"].astype(int)
    df["source"] = "EPSD"
    return df[["uniprot_id", "position", "aa", "pubmed_ids", "source"]].copy()


def load_dbptm(path, human_accessions: Set[str]) -> pd.DataFrame:
    """Load dbPTM phosphorylation sites, filtered to human proteins.

    Parameters
    ----------
    path : str or Path
        Path to the dbPTM 'Phosphorylation' file (tab-delimited, no header).
    human_accessions : set of str
        UniProt accessions known to be human; rows whose uniprot_id is not in
        this set are dropped.

    Returns
    -------
    DataFrame with columns: uniprot_id, position, aa, pubmed_ids, source
    """
    path = Path(path)
    col_names = [
        "entry_name", "uniprot_id", "position",
        "ptm_type", "pubmed_ids", "context",
    ]
    df = pd.read_csv(
        path, sep="\t", header=None, names=col_names,
        dtype=str, low_memory=False,
    )
    df = df[df["uniprot_id"].isin(human_accessions)].copy()
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df = df.dropna(subset=["position"])
    df["position"] = df["position"].astype(int)

    # Sequence context is 21 chars: 10 upstream + modified residue + 10 downstream.
    # The modified residue is the character at index 10 (0-based).
    def _center_aa(ctx):
        if isinstance(ctx, str) and len(ctx) >= 11:
            return ctx[10]
        return ""

    df["aa"] = df["context"].apply(_center_aa)
    df["source"] = "dbPTM"
    return df[["uniprot_id", "position", "aa", "pubmed_ids", "source"]].copy()


def load_uniprot_modres(path) -> tuple:
    """Load UniProt curated phosphorylation sites and human accession set.

    Parameters
    ----------
    path : str or Path
        Path to human_reviewed_modres.tsv (tab-delimited, has header).
        Columns: Entry, Entry Name, Modified residue

    Returns
    -------
    sites : DataFrame
        Columns: uniprot_id, position, aa, pubmed_ids, source
    human_accessions : set of str
        All UniProt accessions in the file.  Pass this to load_dbptm() to
        filter the all-species dbPTM file down to human proteins.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    human_accessions = set(df["Entry"].dropna().str.strip())

    rows = []
    for _, row in df.iterrows():
        acc = row["Entry"]
        modres = row.get("Modified residue", "")
        if not isinstance(modres, str) or not modres.strip():
            continue
        rows.extend(_parse_modres_field(acc, modres))

    sites = (
        pd.DataFrame(rows)
        if rows
        else pd.DataFrame(
            columns=["uniprot_id", "position", "aa", "pubmed_ids", "source"]
        )
    )
    return sites, human_accessions


# ── Annotation ────────────────────────────────────────────────────────────────

def annotate_candidates(
    candidates: pd.DataFrame,
    epsd: Optional[pd.DataFrame] = None,
    dbptm: Optional[pd.DataFrame] = None,
    uniprot: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Join pipeline biomarker candidates with known phosphorylation sites.

    For each candidate S/T residue, look up whether its
    (uniprot_id, position) pair appears in any of the supplied databases.

    Parameters
    ----------
    candidates : DataFrame
        Accepts two formats:

        - Per-protein ``biomarker_candidates.csv``: columns include
          ``sequence_id``, ``position``, ``aa``, ``contrib_total``.
          ``uniprot_id`` is parsed from ``sequence_id``.

        - Cohort ``master_biomarker_candidates.csv``: columns include
          ``uniprot_id``, ``protein_name``, ``region``, ``sequence_id``,
          ``position``, ``aa``, ``contrib_L``, ``contrib_B``, ``contrib_total``.
          ``uniprot_id`` is already present and used directly.

    epsd, dbptm, uniprot : DataFrame or None
        Site DataFrames returned by load_epsd / load_dbptm /
        load_uniprot_modres.  Pass None to skip a source.

    Returns
    -------
    DataFrame with one row per (candidate × source) match.
    Columns: sequence_id, uniprot_id, protein_name, region,
             position, aa, contrib_L, contrib_B, contrib_total,
             ptm_source, ptm_pubmed_ids, ptm_aa
    Sorted by contrib_total descending.
    """
    if candidates.empty:
        return pd.DataFrame()

    cand = candidates.copy()

    # If uniprot_id is not already a column, parse it from sequence_id.
    # master_biomarker_candidates.csv already has uniprot_id, protein_name,
    # and region; per-protein biomarker_candidates.csv does not.
    if "uniprot_id" not in cand.columns:
        parsed = cand["sequence_id"].apply(
            lambda s: pd.Series(
                _parse_sequence_id(s),
                index=["uniprot_id", "protein_name", "region"],
            )
        )
        cand = pd.concat([cand, parsed], axis=1)

    sources = {}
    if epsd is not None and not epsd.empty:
        sources["EPSD"] = epsd
    if dbptm is not None and not dbptm.empty:
        sources["dbPTM"] = dbptm
    if uniprot is not None and not uniprot.empty:
        sources["UniProt"] = uniprot

    if not sources:
        warnings.warn(
            "No PTM source databases provided — returning empty annotation.",
            stacklevel=2,
        )
        return pd.DataFrame()

    hit_frames = []
    for source_name, src_df in sources.items():
        merged = cand.merge(
            src_df[["uniprot_id", "position", "aa", "pubmed_ids"]].rename(
                columns={
                    "aa":        "ptm_aa",
                    "pubmed_ids": "ptm_pubmed_ids",
                }
            ),
            on=["uniprot_id", "position"],
            how="inner",
        )
        merged["ptm_source"] = source_name
        hit_frames.append(merged)

    if not hit_frames:
        return pd.DataFrame()

    hits = pd.concat(hit_frames, ignore_index=True)
    col_order = [
        "sequence_id", "uniprot_id", "protein_name", "region",
        "position", "aa", "contrib_L", "contrib_B", "contrib_total",
        "ptm_source", "ptm_pubmed_ids", "ptm_aa",
    ]
    hits = hits[[c for c in col_order if c in hits.columns]]
    return (
        hits.sort_values(["contrib_total", "position"], ascending=[False, True])
        .reset_index(drop=True)
    )


def summarize_annotations(hits: pd.DataFrame) -> pd.DataFrame:
    """Collapse per-source hits to one summary row per candidate residue.

    Parameters
    ----------
    hits : DataFrame returned by annotate_candidates.

    Returns
    -------
    DataFrame with columns:
        sequence_id, uniprot_id, protein_name, region,
        position, aa, contrib_L, contrib_B, contrib_total,
        n_sources, sources, n_pubmed, pubmed_ids
    Sorted by n_sources (desc) then contrib_total (desc).
    ``n_pubmed`` is the count of unique PubMed IDs across all supporting
    databases.  ``pubmed_ids`` lists them semicolon-separated.
    """
    if hits.empty:
        return pd.DataFrame()

    key_cols = [
        c for c in [
            "sequence_id", "uniprot_id", "protein_name", "region",
            "position", "aa", "contrib_L", "contrib_B", "contrib_total",
        ]
        if c in hits.columns
    ]

    rows = []
    for key, grp in hits.groupby(key_cols, sort=False):
        key_vals = key if isinstance(key, tuple) else (key,)
        sources_list = sorted(grp["ptm_source"].unique())
        all_pmids: set = set()
        for pmid_str in grp["ptm_pubmed_ids"].dropna():
            all_pmids.update(pmid_str.split(";"))
        all_pmids.discard("")

        sorted_pmids = sorted(all_pmids)
        row = dict(zip(key_cols, key_vals))
        row["n_sources"]  = len(sources_list)
        row["sources"]    = ";".join(sources_list)
        row["n_pubmed"]   = len(sorted_pmids)
        row["pubmed_ids"] = ";".join(sorted_pmids)
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary = summary.sort_values(
        ["n_sources", "contrib_total"], ascending=[False, False]
    ).reset_index(drop=True)
    return summary
