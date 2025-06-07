#!/usr/bin/env python3
"""
Convert ProsperousPlus results to website-standard column format
"""

import pandas as pd
import argparse
import os

def convert_results_format(input_file, output_file=None):
    """
    Convert ProsperousPlus results CSV to website-standard format
    
    Website format columns:
    - Rank
    - Protease  
    - Sequence Id
    - Position
    - Cleavage site
    - Prediction Score
    """
    
    # Load original results
    df = pd.read_csv(input_file)
    
    # Create new DataFrame with website column names
    converted_df = pd.DataFrame()
    
    # Map columns to website format
    if 'protease' in df.columns:
        converted_df['Protease'] = df['protease']
    elif 'Protease' in df.columns:
        converted_df['Protease'] = df['Protease']
    
    if 'sequence_id' in df.columns:
        converted_df['Sequence Id'] = df['sequence_id']
    elif 'Sequence Id' in df.columns:
        converted_df['Sequence Id'] = df['Sequence Id']
    
    if 'position' in df.columns:
        converted_df['Position'] = df['position']
    elif 'Position' in df.columns:
        converted_df['Position'] = df['Position']
    
    if 'seqs' in df.columns:
        converted_df['Cleavage site'] = df['seqs']
    elif 'Cleavage site' in df.columns:
        converted_df['Cleavage site'] = df['Cleavage site']
    
    if 'pro' in df.columns:
        converted_df['Prediction Score'] = df['pro']
    elif 'Prediction Score' in df.columns:
        converted_df['Prediction Score'] = df['Prediction Score']
    
    if 'prediction' in df.columns:
        converted_df['Prediction'] = df['prediction']
    elif 'Prediction' in df.columns:
        converted_df['Prediction'] = df['Prediction']
    
    # Sort by prediction score (descending) and add rank
    converted_df = converted_df.sort_values('Prediction Score', ascending=False).reset_index(drop=True)
    converted_df.insert(0, 'Rank', range(1, len(converted_df) + 1))
    
    # Output file
    if output_file is None:
        output_file = input_file.replace('.csv', '_website_format.csv')
    
    # Save converted results
    converted_df.to_csv(output_file, index=False)
    
    print(f"✅ Converted results saved to: {output_file}")
    print(f"📊 {len(converted_df)} predictions converted")
    
    # Show preview
    print(f"\n📋 Preview of converted format:")
    print(converted_df.head().to_string(index=False))
    
    return output_file

def main():
    parser = argparse.ArgumentParser(description='Convert ProsperousPlus results to website format')
    parser.add_argument('--input', '-i', required=True, help='Input results CSV file')
    parser.add_argument('--output', '-o', help='Output file (optional, will auto-generate if not provided)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found!")
        return
    
    convert_results_format(args.input, args.output)

if __name__ == "__main__":
    main() 