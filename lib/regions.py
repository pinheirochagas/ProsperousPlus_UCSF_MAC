"""Find low-cleavage regions - adapted from older/scripts/analyze_low_score_regions.py"""

import pandas as pd
from collections import defaultdict


def find_low_cleavage_regions(df1, df2, protease1, protease2, min_length=25, max_exceptions=2, threshold=0.3):
    """Find regions where both proteases have low scores."""
    merged = pd.merge(
        df1[['sequence_id', 'position', 'pro', 'seqs']], 
        df2[['sequence_id', 'position', 'pro']], 
        on=['sequence_id', 'position'], 
        suffixes=(f'_{protease1}', f'_{protease2}')
    )
    
    regions = []
    for sequence_id, group in merged.groupby('sequence_id'):
        group = group.sort_values('position').reset_index(drop=True)
        sequence_regions = _find_regions_in_sequence(
            group, sequence_id, protease1, protease2, min_length, max_exceptions, threshold
        )
        regions.extend(sequence_regions)
    
    regions = _remove_overlapping_regions(regions)
    return pd.DataFrame(regions)


def _find_regions_in_sequence(group, sequence_id, protease1, protease2, min_length, max_exceptions, threshold):
    """Find qualifying regions within a single sequence."""
    regions = []
    n = len(group)
    
    for start_idx in range(n - min_length + 1):
        for end_idx in range(start_idx + min_length - 1, n):
            window = group.iloc[start_idx:end_idx + 1]
            
            p1_exceptions = sum(window[f'pro_{protease1}'] >= threshold)
            p2_exceptions = sum(window[f'pro_{protease2}'] >= threshold)
            total_exceptions = p1_exceptions + p2_exceptions
            
            if total_exceptions <= max_exceptions:
                continue
            else:
                if end_idx > start_idx + min_length - 1:
                    prev_window = group.iloc[start_idx:end_idx]
                    prev_p1 = sum(prev_window[f'pro_{protease1}'] >= threshold)
                    prev_p2 = sum(prev_window[f'pro_{protease2}'] >= threshold)
                    
                    if prev_p1 + prev_p2 <= max_exceptions:
                        regions.append(_create_region_info(prev_window, sequence_id))
                break
        else:
            window = group.iloc[start_idx:]
            p1_exc = sum(window[f'pro_{protease1}'] >= threshold)
            p2_exc = sum(window[f'pro_{protease2}'] >= threshold)
            
            if p1_exc + p2_exc <= max_exceptions and len(window) >= min_length:
                regions.append(_create_region_info(window, sequence_id))
    
    return regions


def _create_region_info(window, sequence_id):
    """Create region info dictionary."""
    return {
        'sequence_id': sequence_id,
        'start_position': int(window['position'].iloc[0]),
        'end_position': int(window['position'].iloc[-1]),
        'length': len(window),
        'sequence_region': ''.join(window['seqs'].str[4])
    }


def _remove_overlapping_regions(regions):
    """Remove overlapping regions, keeping the longest ones."""
    if not regions:
        return regions
    
    sequence_groups = defaultdict(list)
    for region in regions:
        sequence_groups[region['sequence_id']].append(region)
    
    filtered_regions = []
    
    for sequence_id, seq_regions in sequence_groups.items():
        seq_regions.sort(key=lambda x: x['start_position'])
        
        non_overlapping = []
        for region in seq_regions:
            overlaps = False
            for existing in non_overlapping:
                if (region['start_position'] <= existing['end_position'] and 
                    region['end_position'] >= existing['start_position']):
                    if region['length'] > existing['length']:
                        non_overlapping.remove(existing)
                        non_overlapping.append(region)
                    overlaps = True
                    break
            
            if not overlaps:
                non_overlapping.append(region)
        
        filtered_regions.extend(non_overlapping)
    
    return filtered_regions
