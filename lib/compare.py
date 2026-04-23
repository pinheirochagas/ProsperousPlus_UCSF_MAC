"""Compare predictions and calculate mutation sums."""

import pandas as pd


def calculate_single_residue_mutation_sum(orig_L, single_P, orig_B, single_E, st_positions):
    """Calculate per-residue contributions and window-level mutation sums.

    For each S/T position *k* in each extended region the contribution is
    measured as the drop in cleavage score when only that single residue is
    mutated (S/T→P for CathL, S/T→E for CathB):

        contrib_k_L = cathL_orig[window_k] − cathL_single_P_k[window_k]
        contrib_k_B = cathB_orig[window_k] − cathB_single_E_k[window_k]

    where ``window_k`` is the 8-mer window whose cleavage site falls on
    position *k*.  Window-level mutation scores aggregate the per-residue
    contributions of every S/T within that window.

    Parameters
    ----------
    orig_L, orig_B   : DataFrames from CathL/CathB predictions on original sequences
    single_P, single_E : DataFrames from CathL/CathB on single-mutant sequences
    st_positions     : dict returned by ``create_single_mutant_fastas``
                       ``{region_id: [(k, aa, abs_pos), ...]}``

    Returns
    -------
    mutation_scores : DataFrame  (one row per 8-mer window; same schema as before)
    residue_contribs : DataFrame (one row per S/T position)
    """
    mutation_score_rows = []
    residue_rows = []

    for region_id, st_list in st_positions.items():
        oL = orig_L[orig_L['sequence_id'] == region_id].sort_values('position').reset_index(drop=True)
        oB = orig_B[orig_B['sequence_id'] == region_id].sort_values('position').reset_index(drop=True)

        if len(oL) == 0:
            continue

        n_windows = len(oL)
        L = n_windows + 7  # extended region length

        parts = region_id.split('_')
        start_offset = int(parts[-1].split('-')[0])

        # Per-residue contributions: contrib_by_k[k] = (contrib_L, contrib_B)
        contrib_by_k = {}
        for k, aa, abs_pos in st_list:
            mutant_id = f"{region_id}_st{k}"
            sP = single_P[single_P['sequence_id'] == mutant_id].sort_values('position').reset_index(drop=True)
            sE = single_E[single_E['sequence_id'] == mutant_id].sort_values('position').reset_index(drop=True)

            # Window where position k is the cleavage site (0-based index)
            window_idx = max(0, min(k - 4, n_windows - 1))

            if window_idx < len(sP) and window_idx < len(sE):
                contrib_L = float(oL.iloc[window_idx]['pro']) - float(sP.iloc[window_idx]['pro'])
                contrib_B = float(oB.iloc[window_idx]['pro']) - float(sE.iloc[window_idx]['pro'])
            else:
                contrib_L = contrib_B = 0.0

            contrib_by_k[k] = (contrib_L, contrib_B)

            residue_rows.append({
                'region_id':   region_id,
                'abs_position': abs_pos,
                'aa':           aa,
                'contrib_L':    round(contrib_L, 4),
                'contrib_B':    round(contrib_B, 4),
                'contrib_total': round(contrib_L + contrib_B, 4),
            })

        # Aggregate to per-window scores
        for i in range(n_windows):
            actual_position = start_offset + i + 3  # 8-mer center (matches original convention)
            # Positions within this window are i+1 … i+8 (1-indexed in extended region)
            w_start, w_end = i + 1, i + 8
            cathL_diff = sum(
                contrib_by_k[k][0] for k, aa, abs_pos in st_list if w_start <= k <= w_end
            )
            cathB_diff = sum(
                contrib_by_k[k][1] for k, aa, abs_pos in st_list if w_start <= k <= w_end
            )
            mutation_score_rows.append({
                'sequence_id':    region_id,
                'position':       actual_position,
                'peptide':        oL.iloc[i]['seqs'],
                'cathL_original': round(float(oL.iloc[i]['pro']), 4),
                'cathB_original': round(float(oB.iloc[i]['pro']), 4),
                'cathL_diff':     round(cathL_diff, 4),
                'cathB_diff':     round(cathB_diff, 4),
                'mutation_sum':   round(cathL_diff + cathB_diff, 4),
            })

    return pd.DataFrame(mutation_score_rows), pd.DataFrame(residue_rows)


def calculate_mutation_sum(orig_L, mut_P, orig_B, mut_E):
    """Calculate mutation sum for each position."""
    results = []
    
    for seq_id in orig_L['sequence_id'].unique():
        oL = orig_L[orig_L['sequence_id'] == seq_id].sort_values('position').reset_index(drop=True)
        oB = orig_B[orig_B['sequence_id'] == seq_id].sort_values('position').reset_index(drop=True)
        
        # Match mutant by removing suffix
        mut_seq_id = seq_id + "_P"
        mP = mut_P[mut_P['sequence_id'] == mut_seq_id].sort_values('position').reset_index(drop=True)
        
        mut_seq_id_E = seq_id + "_E"
        mE = mut_E[mut_E['sequence_id'] == mut_seq_id_E].sort_values('position').reset_index(drop=True)
        
        # Extract start position from seq_id (format: sp|...|NAME_START-END)
        parts = seq_id.split('_')
        pos_range = parts[-1]  # e.g., "419-461"
        start_offset = int(pos_range.split('-')[0])
        
        for i in range(min(len(oL), len(mP), len(oB), len(mE))):
            actual_position = start_offset + i + 3  # +3 for 8-mer center offset
            
            cathL_diff = oL.iloc[i]['pro'] - mP.iloc[i]['pro']
            cathB_diff = oB.iloc[i]['pro'] - mE.iloc[i]['pro']
            
            results.append({
                'sequence_id': seq_id,
                'position': actual_position,
                'peptide': oL.iloc[i]['seqs'],
                'cathL_original': round(oL.iloc[i]['pro'], 4),
                'cathL_P_mutant': round(mP.iloc[i]['pro'], 4),
                'cathL_diff': round(cathL_diff, 4),
                'cathB_original': round(oB.iloc[i]['pro'], 4),
                'cathB_E_mutant': round(mE.iloc[i]['pro'], 4),
                'cathB_diff': round(cathB_diff, 4),
                'mutation_sum': round(cathL_diff + cathB_diff, 4)
            })
    
    return pd.DataFrame(results)


def find_candidates(df, threshold=0.45):
    """Find positions with mutation_sum >= threshold."""
    return df[df['mutation_sum'] >= threshold].sort_values('mutation_sum', ascending=False)
