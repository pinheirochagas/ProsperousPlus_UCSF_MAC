"""Create mutant sequences and save FASTAs."""

import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_sequences(fasta_file):
    """Load sequences from FASTA into dict."""
    return {r.id: str(r.seq) for r in SeqIO.parse(fasta_file, "fasta")}


def extend_region(full_seq, start, end, extension=5):
    """Extend region by N amino acids on each side."""
    new_start = max(0, start - 1 - extension)
    new_end = min(len(full_seq), end + extension)
    return full_seq[new_start:new_end], new_start + 1, new_end


def mutate_ST_to_P(sequence):
    """Replace S and T with P."""
    return sequence.replace('S', 'P').replace('T', 'P')


def mutate_ST_to_E(sequence):
    """Replace S and T with E."""
    return sequence.replace('S', 'E').replace('T', 'E')


def create_single_mutant_fastas(
    regions_df,
    full_sequences,
    output_dir,
    extension: int = 5,
    orig_name: str = "original.fasta",
    mutp_name: str = "single_ST_to_P.fasta",
    mute_name: str = "single_ST_to_E.fasta",
):
    """Create FASTA files for single-position S/T mutations.

    For each S or T residue in every extended region, one mutant sequence is
    written: the residue is changed to P (for CathL) or E (for CathB) while
    every other position is left unchanged.  Record IDs use the scheme
    ``{region_id}_st{k}`` where *k* is the 1-based position within the
    extended region.

    FASTAs are only written if they do not already exist, so resume works
    without re-running this function.

    Returns
    -------
    st_positions : dict
        ``{region_id: [(k, aa, abs_pos), ...]}`` — for each region, the list
        of mutable positions where *k* is 1-based within the extended region,
        *aa* is the original amino acid (S or T), and *abs_pos* is the
        1-based position in the full input protein.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    orig_out = output_path / orig_name
    mutp_out = output_path / mutp_name
    mute_out = output_path / mute_name
    need_write = not (orig_out.exists() and mutp_out.exists() and mute_out.exists())

    original_records = []
    single_P_records = []
    single_E_records = []
    st_positions = {}

    for _, row in regions_df.iterrows():
        seq_id   = row['sequence_id']
        full_seq = full_sequences[seq_id]

        extended, new_start, new_end = extend_region(
            full_seq, row['start_position'], row['end_position'], extension
        )

        region_id = f"{seq_id}_{new_start}-{new_end}"
        original_records.append(SeqRecord(Seq(extended), id=region_id, description=""))

        st_list = []
        for k, aa in enumerate(extended, 1):
            if aa in ('S', 'T'):
                mutant_id = f"{region_id}_st{k}"
                single_P_records.append(
                    SeqRecord(Seq(extended[:k-1] + 'P' + extended[k:]), id=mutant_id, description="")
                )
                single_E_records.append(
                    SeqRecord(Seq(extended[:k-1] + 'E' + extended[k:]), id=mutant_id, description="")
                )
                st_list.append((k, aa, new_start + k - 1))

        st_positions[region_id] = st_list

    if need_write:
        SeqIO.write(original_records, orig_out, "fasta")
        SeqIO.write(single_P_records, mutp_out, "fasta")
        SeqIO.write(single_E_records, mute_out, "fasta")

    return st_positions


def create_mutation_fastas(
    regions_df,
    full_sequences,
    output_dir,
    extension: int = 5,
    orig_name: str = "original.fasta",
    mutp_name: str = "mutant_ST_to_P.fasta",
    mute_name: str = "mutant_ST_to_E.fasta",
):
    """Create FASTA files for original and mutant sequences.

    Parameters
    ----------
    orig_name / mutp_name / mute_name:
        Output filenames for the three FASTA files (default names match the
        descriptive convention used by run_pipeline.py).
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    original_records = []
    mutant_P_records = []
    mutant_E_records = []

    for _, row in regions_df.iterrows():
        seq_id   = row['sequence_id']
        full_seq = full_sequences[seq_id]

        extended, new_start, new_end = extend_region(
            full_seq, row['start_position'], row['end_position'], extension
        )

        record_id = f"{seq_id}_{new_start}-{new_end}"

        original_records.append(SeqRecord(Seq(extended), id=record_id, description=""))
        mutant_P_records.append(SeqRecord(Seq(mutate_ST_to_P(extended)), id=f"{record_id}_P", description=""))
        mutant_E_records.append(SeqRecord(Seq(mutate_ST_to_E(extended)), id=f"{record_id}_E", description=""))

    SeqIO.write(original_records, output_path / orig_name, "fasta")
    SeqIO.write(mutant_P_records, output_path / mutp_name, "fasta")
    SeqIO.write(mutant_E_records, output_path / mute_name, "fasta")

    return output_path

