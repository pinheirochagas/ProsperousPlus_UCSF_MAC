#!/usr/bin/env python3
"""
Simple script to plot protease cleavage prediction results (no GUI display)
Automatically detects and handles both raw and display formats
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

def detect_format(df):
    """
    Detect whether the CSV is in raw format or display format
    """
    if 'Rank' in df.columns and 'Sequence Id' in df.columns and 'Prediction Score' in df.columns:
        return 'display'
    elif 'protease' in df.columns and 'sequence_id' in df.columns and 'pro' in df.columns:
        return 'raw'
    else:
        raise ValueError("Unknown CSV format. Expected either raw format (protease, sequence_id, pro) or display format (Rank, Sequence Id, Prediction Score)")

def standardize_columns(df, format_type):
    """
    Standardize column names for plotting regardless of input format
    """
    if format_type == 'display':
        # Display format -> standard names
        return {
            'protease_col': 'Protease',
            'sequence_col': 'Sequence Id', 
            'position_col': 'Position',
            'score_col': 'Prediction Score',
            'prediction_col': 'Prediction',
            'cleavage_col': 'Cleavage site'
        }
    else:  # raw format
        # Raw format -> standard names
        return {
            'protease_col': 'protease',
            'sequence_col': 'sequence_id',
            'position_col': 'position', 
            'score_col': 'pro',
            'prediction_col': 'prediction',
            'cleavage_col': 'seqs'
        }

def plot_protein_predictions(results_file, output_dir="plots"):
    """
    Create plots with website styling (auto-detects format)
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load results and detect format
    df = pd.read_csv(results_file)
    format_type = detect_format(df)
    cols = standardize_columns(df, format_type)
    
    print(f"📋 Detected format: {format_type}")
    
    # Get unique proteins and proteases using detected column names
    proteins = df[cols['sequence_col']].unique()
    proteases = df[cols['protease_col']].unique()
    
    print(f"Processing {len(proteins)} proteins with {len(proteases)} protease(s)")
    
    # Create individual plots for each protein (website style)
    for i, protein in enumerate(proteins):
        protein_data = df[df[cols['sequence_col']] == protein].copy()
        
        # More landscape aspect ratio to match website
        fig, ax = plt.subplots(figsize=(16, 5))
        
        for j, protease in enumerate(proteases):
            protease_data = protein_data[protein_data[cols['protease_col']] == protease].copy()
            
            if len(protease_data) > 0:
                # Sort by position for proper line plotting
                protease_data = protease_data.sort_values(cols['position_col'])
                
                positions = protease_data[cols['position_col']].astype(int)
                scores = protease_data[cols['score_col']].astype(float)
                
                # Plot line only (no markers) with website styling - black line
                ax.plot(positions, scores, '-', color='black', 
                       linewidth=1, alpha=1)
        
        # Customize plot to match website exactly
        ax.set_xlabel('Residue Position', fontsize=20)
        ax.set_ylabel('Prediction Score', fontsize=20) 
        ax.set_title(f'ProsperousPlus Prediction for {protein}', fontsize=14, fontweight='bold')
        
        # Set larger fontsize for tick labels
        ax.tick_params(axis='both', which='major', labelsize=16)
        
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0.0, 1.0)
        
        # Add threshold line at 0.5 (dashed green line like website)
        ax.axhline(y=0.5, color='green', linestyle='--', alpha=0.6, linewidth=1)
        
        # Set background to white and make more landscape
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(output_dir, f'{protein}_cleavage_prediction.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"  ✓ Saved plot: {output_file}")
    
    # Create comparison plot if multiple proteins
    # if len(proteins) > 1:
    #     # Even more landscape for comparison
    #     fig, ax = plt.subplots(figsize=(16, 4))
        
    #     colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
        
    #     for i, protein in enumerate(proteins):
    #         protein_data = df[df[cols['sequence_col']] == protein].copy()
            
    #         # Use first protease for comparison (or combine if multiple)
    #         if len(proteases) == 1:
    #             protease_data = protein_data.copy()
    #         else:
    #             # Average scores across proteases for comparison
    #             protease_data = protein_data.groupby(cols['position_col']).agg({
    #                 cols['score_col']: 'mean',
    #                 cols['sequence_col']: 'first'
    #             }).reset_index()
            
    #         if len(protease_data) > 0:
    #             protease_data = protease_data.sort_values(cols['position_col'])
    #             positions = protease_data[cols['position_col']].astype(int)
    #             scores = protease_data[cols['score_col']].astype(float)
                
    #             color = colors[i % len(colors)]
    #             # Line only, no markers for comparison too
    #             ax.plot(positions, scores, '-', label=f'{protein}', color=color, 
    #                    linewidth=2, alpha=0.7)
        
    #     ax.set_xlabel('Residue Position', fontsize=20)
    #     ax.set_ylabel('Prediction Score', fontsize=20)
    #     ax.set_title('ProsperousPlus Prediction - All Proteins Comparison', fontsize=14, fontweight='bold')
    #     ax.grid(True, alpha=0.3)
    #     ax.set_ylim(0.0, 1.0)
    #     ax.legend(fontsize=10)
    #     ax.axhline(y=0.5, color='green', linestyle='--', alpha=0.6, linewidth=1)
        
    #     # Set background to white
    #     ax.set_facecolor('white')
    #     fig.patch.set_facecolor('white')
        
    #     plt.tight_layout()
        
    #     comparison_file = os.path.join(output_dir, 'all_proteins_comparison.png')
    #     plt.savefig(comparison_file, dpi=300, bbox_inches='tight', facecolor='white')
    #     plt.close()
        
    #     print(f"  ✓ Saved comparison plot: {comparison_file}")
    
    print(f"\n📊 All plots saved to: {output_dir}")

def print_top_sites(results_file, top_n=10):
    """
    Print summary of top cleavage sites (auto-detects format)
    """
    df = pd.read_csv(results_file)
    format_type = detect_format(df)
    cols = standardize_columns(df, format_type)
    
    print(f"\n📋 TOP {top_n} CLEAVAGE SITES (sorted by probability) - {format_type.upper()} FORMAT:")
    print("=" * 100)
    
    if format_type == 'display':
        # Display format already has Rank
        print(f"{'Rank':<6} {'Protease':<10} {'Sequence Id':<15} {'Position':<10} {'Cleavage site':<15} {'Prediction Score':<15}")
        print("=" * 100)
        
        # Sort by prediction score and get top sites
        top_sites = df.nlargest(top_n, cols['score_col'])
        
        for _, row in top_sites.iterrows():
            print(f"{row['Rank']:<6} {row[cols['protease_col']]:<10} {row[cols['sequence_col']]:<15} {row[cols['position_col']]:<10} "
                  f"{row[cols['cleavage_col']]:<15} {row[cols['score_col']]:<15.3f}")
    else:
        # Raw format - create ranking
        print(f"{'Rank':<6} {'Protease':<10} {'Sequence Id':<15} {'Position':<10} {'Cleavage site':<15} {'Prediction Score':<15}")
        print("=" * 100)
        
        # Sort by prediction score and get top sites
        top_sites = df.nlargest(top_n, cols['score_col'])
        
        for idx, (_, row) in enumerate(top_sites.iterrows(), 1):
            print(f"{idx:<6} {row[cols['protease_col']]:<10} {row[cols['sequence_col']]:<15} {row[cols['position_col']]:<10} "
                  f"{row[cols['cleavage_col']]:<15} {row[cols['score_col']]:<15.3f}")
    
    # Summary statistics
    print(f"\n📊 SUMMARY STATISTICS:")
    print(f"   Total predictions: {len(df)}")
    print(f"   Predicted cleavages: {sum(df[cols['prediction_col']] == 1)}")
    print(f"   Average prediction score: {df[cols['score_col']].mean():.3f}")
    print(f"   Max prediction score: {df[cols['score_col']].max():.3f}")
    print(f"   Sites above 0.5: {sum(df[cols['score_col']] > 0.5)}")
    print(f"   Sites above 0.8: {sum(df[cols['score_col']] > 0.8)}")

def main():
    parser = argparse.ArgumentParser(description='Plot protease cleavage prediction results (auto-detects format)')
    parser.add_argument('--results', '-r', required=True, help='Path to results.csv file')
    parser.add_argument('--output', '-o', default='plots', help='Output directory for plots')
    parser.add_argument('--top-sites', '-t', type=int, default=10, help='Number of top sites to display')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results):
        print(f"Error: Results file {args.results} not found!")
        return
    
    print(f"📈 Generating plots from: {args.results}")
    plot_protein_predictions(args.results, args.output)
    print_top_sites(args.results, args.top_sites)

if __name__ == "__main__":
    main() 