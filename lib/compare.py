"""Compare predictions and calculate mutation sums."""

import pandas as pd


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
