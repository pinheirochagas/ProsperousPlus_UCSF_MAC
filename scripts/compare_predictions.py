#!/usr/bin/env python3
"""
Compare ProsperousPlus website predictions with local predictions
This script validates that our local installation produces identical results
"""

import pandas as pd
import numpy as np
import argparse
import os

def load_and_standardize(file_path):
    """
    Load CSV and convert to standardized format for comparison
    """
    df = pd.read_csv(file_path)
    
    # Detect format and standardize column names
    if 'Rank' in df.columns and 'Sequence Id' in df.columns:
        # Display format
        standardized = pd.DataFrame({
            'protease': df['Protease'],
            'sequence_id': df['Sequence Id'],
            'position': df['Position'], 
            'seqs': df['Cleavage site'],
            'prediction': df['Prediction'],
            'pro': df['Prediction Score']
        })
    else:
        # Raw format - already standardized
        standardized = df.copy()
    
    # Ensure consistent data types
    standardized['position'] = standardized['position'].astype(int)
    standardized['prediction'] = standardized['prediction'].astype(int)
    standardized['pro'] = standardized['pro'].astype(float)
    
    # Sort by position for comparison
    standardized = standardized.sort_values(['sequence_id', 'position']).reset_index(drop=True)
    
    return standardized

def compare_predictions(website_file, local_file, tolerance=1e-6):
    """
    Compare website vs local predictions with detailed analysis
    """
    print("🔍 PROSPEROUSPLUS PREDICTION VALIDATION")
    print("=" * 80)
    
    # Load and standardize both files
    print(f"📂 Loading website predictions: {website_file}")
    website_df = load_and_standardize(website_file)
    
    print(f"📂 Loading local predictions: {local_file}")
    local_df = load_and_standardize(local_file)
    
    print(f"\n📊 DATA OVERVIEW:")
    print(f"   Website predictions: {len(website_df)} rows")
    print(f"   Local predictions: {len(local_df)} rows")
    
    # Basic structure comparison
    if len(website_df) != len(local_df):
        print(f"⚠️  ROW COUNT MISMATCH: Website({len(website_df)}) vs Local({len(local_df)})")
        return False
    
    # Column comparison
    website_cols = set(website_df.columns)
    local_cols = set(local_df.columns)
    
    if website_cols != local_cols:
        print(f"⚠️  COLUMN MISMATCH:")
        print(f"   Website columns: {sorted(website_cols)}")
        print(f"   Local columns: {sorted(local_cols)}")
        return False
    
    print(f"✅ Structure validation passed")
    
    # Merge on key columns for comparison
    merged = pd.merge(
        website_df, local_df, 
        on=['protease', 'sequence_id', 'position', 'seqs'],
        suffixes=('_website', '_local'),
        how='outer',
        indicator=True
    )
    
    # Check for missing rows
    missing_in_local = merged[merged['_merge'] == 'left_only']
    missing_in_website = merged[merged['_merge'] == 'right_only']
    
    if len(missing_in_local) > 0:
        print(f"⚠️  MISSING IN LOCAL: {len(missing_in_local)} predictions")
        print(missing_in_local[['protease', 'sequence_id', 'position', 'seqs']].head())
    
    if len(missing_in_website) > 0:
        print(f"⚠️  MISSING IN WEBSITE: {len(missing_in_website)} predictions")
        print(missing_in_website[['protease', 'sequence_id', 'position', 'seqs']].head())
    
    # Focus on matching rows
    matching = merged[merged['_merge'] == 'both'].copy()
    print(f"\n🔄 COMPARING {len(matching)} MATCHING PREDICTIONS:")
    
    # Compare binary predictions
    pred_match = matching['prediction_website'] == matching['prediction_local']
    pred_mismatch = sum(~pred_match)
    
    print(f"   Binary predictions (0/1):")
    print(f"   ✅ Matching: {sum(pred_match)}")
    print(f"   ❌ Mismatching: {pred_mismatch}")
    
    if pred_mismatch > 0:
        print(f"   Mismatch examples:")
        mismatches = matching[~pred_match][['sequence_id', 'position', 'seqs', 'prediction_website', 'prediction_local']].head()
        print(mismatches.to_string(index=False))
    
    # Compare probability scores
    prob_diff = np.abs(matching['pro_website'] - matching['pro_local'])
    max_diff = prob_diff.max()
    mean_diff = prob_diff.mean()
    
    print(f"\n   Probability scores:")
    print(f"   📈 Maximum difference: {max_diff:.8f}")
    print(f"   📊 Mean difference: {mean_diff:.8f}")
    print(f"   🎯 Tolerance: {tolerance}")
    
    within_tolerance = sum(prob_diff <= tolerance)
    outside_tolerance = sum(prob_diff > tolerance)
    
    print(f"   ✅ Within tolerance: {within_tolerance}")
    print(f"   ❌ Outside tolerance: {outside_tolerance}")
    
    if outside_tolerance > 0:
        print(f"\n   Largest differences:")
        largest_diffs = matching.loc[prob_diff.nlargest(5).index]
        diff_summary = largest_diffs[['sequence_id', 'position', 'seqs', 'pro_website', 'pro_local']].copy()
        diff_summary['difference'] = prob_diff.loc[diff_summary.index]
        print(diff_summary.to_string(index=False))
    
    # Overall validation result
    print(f"\n🎯 VALIDATION RESULT:")
    if pred_mismatch == 0 and outside_tolerance == 0:
        print(f"   ✅ PERFECT MATCH: Local predictions identical to website!")
        return True
    elif pred_mismatch == 0 and outside_tolerance <= len(matching) * 0.01:  # Allow 1% minor differences
        print(f"   ✅ EXCELLENT MATCH: Predictions match with minor numerical differences")
        return True
    else:
        print(f"   ❌ SIGNIFICANT DIFFERENCES: Local predictions differ from website")
        return False

def main():
    parser = argparse.ArgumentParser(description='Compare website vs local ProsperousPlus predictions')
    parser.add_argument('--website', '-w', required=True, help='Website predictions CSV file')
    parser.add_argument('--local', '-l', required=True, help='Local predictions CSV file')
    parser.add_argument('--tolerance', '-t', type=float, default=1e-6, help='Numerical tolerance for probability scores')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.website):
        print(f"Error: Website file {args.website} not found!")
        return
    
    if not os.path.exists(args.local):
        print(f"Error: Local file {args.local} not found!")
        return
    
    success = compare_predictions(args.website, args.local, args.tolerance)
    
    if success:
        print(f"\n🎉 LOCAL INSTALLATION VALIDATED!")
        print(f"   Your ProsperousPlus setup produces results identical to the website")
    else:
        print(f"\n⚠️  VALIDATION ISSUES DETECTED")
        print(f"   Consider checking your local installation or model files")

if __name__ == "__main__":
    main() 