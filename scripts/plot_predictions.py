#!/usr/bin/env python3
"""
Plot protease cleavage prediction results
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import argparse

def load_and_plot_predictions(results_file, output_dir="plots", threshold=0.5):
    """
    Load prediction results and create plots for each protein
    
    Args:
        results_file: Path to the results.csv file
        output_dir: Directory to save plots
        threshold: Threshold for highlighting high-probability cleavage sites
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the results
    df = pd.read_csv(results_file)
    
    # Get unique proteins
    proteins = df['sequence_id'].unique()
    proteases = df['protease'].unique()
    
    print(f"Found {len(proteins)} proteins and {len(proteases)} proteases")
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create plots for each protein
    for protein in proteins:
        protein_data = df[df['sequence_id'] == protein].copy()
        protein_data = protein_data.sort_values('position')
        
        # Create figure with subplots for each protease
        n_proteases = len(proteases)
        fig, axes = plt.subplots(n_proteases, 1, figsize=(12, 4*n_proteases), squeeze=False)
        
        for i, protease in enumerate(proteases):
            ax = axes[i, 0]
            protease_data = protein_data[protein_data['protease'] == protease]
            
            if len(protease_data) == 0:
                continue
                
            positions = protease_data['position']
            scores = protease_data['pro']
            predictions = protease_data['prediction']
            
            # Plot the prediction scores as a line
            ax.plot(positions, scores, 'o-', linewidth=2, markersize=4, alpha=0.7, label='Prediction Score')
            
            # Highlight cleavage sites (prediction = 1)
            cleavage_sites = protease_data[protease_data['prediction'] == 1]
            if len(cleavage_sites) > 0:
                ax.scatter(cleavage_sites['position'], cleavage_sites['pro'], 
                          color='red', s=50, alpha=0.8, label='Predicted Cleavage', zorder=5)
            
            # Add threshold line
            ax.axhline(y=threshold, color='gray', linestyle='--', alpha=0.5, label=f'Threshold ({threshold})')
            
            # Customize the plot
            ax.set_xlabel('Residue Position')
            ax.set_ylabel('Cleavage Probability')
            ax.set_title(f'{protein} - {protease} Cleavage Predictions')
            ax.set_ylim(0, 1)
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Add text annotations for high-probability sites
            high_prob_sites = protease_data[protease_data['pro'] > threshold]
            for _, site in high_prob_sites.iterrows():
                ax.annotate(f'{site["seqs"]}\n{site["pro"]:.3f}', 
                           xy=(site['position'], site['pro']),
                           xytext=(10, 10), textcoords='offset points',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                           arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
                           fontsize=8)
        
        plt.tight_layout()
        
        # Save the plot
        output_file = os.path.join(output_dir, f'{protein}_cleavage_predictions.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot for {protein}: {output_file}")
        
        plt.show()
    
    # Create a summary plot with all proteins
    create_summary_plot(df, output_dir, threshold)

def create_summary_plot(df, output_dir, threshold):
    """Create a summary plot showing all proteins together"""
    
    proteins = df['sequence_id'].unique()
    proteases = df['protease'].unique()
    
    fig, axes = plt.subplots(len(proteases), 1, figsize=(15, 4*len(proteases)), squeeze=False)
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(proteins)))
    
    for i, protease in enumerate(proteases):
        ax = axes[i, 0]
        protease_data = df[df['protease'] == protease]
        
        for j, protein in enumerate(proteins):
            protein_protease_data = protease_data[protease_data['sequence_id'] == protein]
            
            if len(protein_protease_data) == 0:
                continue
                
            protein_protease_data = protein_protease_data.sort_values('position')
            
            # Adjust positions to avoid overlap between proteins
            adjusted_positions = protein_protease_data['position'] + j * 100
            
            ax.plot(adjusted_positions, protein_protease_data['pro'], 
                   'o-', linewidth=2, markersize=3, alpha=0.7, 
                   color=colors[j], label=f'{protein}')
            
            # Highlight cleavage sites
            cleavage_sites = protein_protease_data[protein_protease_data['prediction'] == 1]
            if len(cleavage_sites) > 0:
                adjusted_cleavage_pos = cleavage_sites['position'] + j * 100
                ax.scatter(adjusted_cleavage_pos, cleavage_sites['pro'], 
                          color=colors[j], s=50, alpha=1.0, 
                          marker='s', edgecolor='black', linewidth=1, zorder=5)
        
        # Add threshold line
        ax.axhline(y=threshold, color='gray', linestyle='--', alpha=0.5, label=f'Threshold ({threshold})')
        
        ax.set_xlabel('Adjusted Residue Position (protein separated by 100)')
        ax.set_ylabel('Cleavage Probability')
        ax.set_title(f'All Proteins - {protease} Cleavage Predictions')
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    # Save the summary plot
    output_file = os.path.join(output_dir, 'all_proteins_summary.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved summary plot: {output_file}")
    
    plt.show()

def print_summary_stats(results_file):
    """Print summary statistics of the predictions"""
    
    df = pd.read_csv(results_file)
    
    print("\n" + "="*60)
    print("PREDICTION SUMMARY STATISTICS")
    print("="*60)
    
    for protease in df['protease'].unique():
        protease_data = df[df['protease'] == protease]
        print(f"\nProtease: {protease}")
        print("-" * 40)
        
        for protein in protease_data['sequence_id'].unique():
            protein_data = protease_data[protease_data['sequence_id'] == protein]
            
            total_sites = len(protein_data)
            cleavage_sites = len(protein_data[protein_data['prediction'] == 1])
            high_conf_sites = len(protein_data[protein_data['pro'] > 0.7])
            avg_score = protein_data['pro'].mean()
            max_score = protein_data['pro'].max()
            
            print(f"  {protein}:")
            print(f"    Total sites analyzed: {total_sites}")
            print(f"    Predicted cleavage sites: {cleavage_sites} ({cleavage_sites/total_sites*100:.1f}%)")
            print(f"    High confidence sites (>0.7): {high_conf_sites}")
            print(f"    Average prediction score: {avg_score:.3f}")
            print(f"    Maximum prediction score: {max_score:.3f}")
            
            # Show top 3 cleavage sites
            top_sites = protein_data.nlargest(3, 'pro')[['position', 'seqs', 'pro', 'prediction']]
            print(f"    Top 3 sites:")
            for _, site in top_sites.iterrows():
                print(f"      Pos {site['position']}: {site['seqs']} (score: {site['pro']:.3f}, pred: {site['prediction']})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot protease cleavage predictions')
    parser.add_argument('--results', '-r', default='example_results/results.csv', 
                       help='Path to results CSV file')
    parser.add_argument('--output', '-o', default='prediction_plots', 
                       help='Output directory for plots')
    parser.add_argument('--threshold', '-t', type=float, default=0.5, 
                       help='Threshold for highlighting cleavage sites')
    parser.add_argument('--stats-only', action='store_true',
                       help='Only print statistics, do not create plots')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results):
        print(f"Error: Results file {args.results} not found!")
        exit(1)
    
    # Print summary statistics
    print_summary_stats(args.results)
    
    if not args.stats_only:
        # Create plots
        load_and_plot_predictions(args.results, args.output, args.threshold)
        print(f"\nPlots saved to: {args.output}/") 