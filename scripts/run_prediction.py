#!/usr/bin/env python3
"""
Organized ProsperousPlus prediction wrapper script
This script uses the organized project directory structure for clean data management
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
from pathlib import Path

# Define paths relative to the project structure
PROJECT_ROOT = Path("/shared/macdata/groups/ppc/projects/ProsperousPlus")
CODE_ROOT = Path("/shared/macdata/groups/ppc/code/ProsperousPlus")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" 
PLOTS_DIR = PROJECT_ROOT / "plots"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

def list_available_proteases():
    """List available pre-trained protease models"""
    models_dir = CODE_ROOT / "finalModels"
    if models_dir.exists():
        proteases = [d.name for d in models_dir.iterdir() if d.is_dir()]
        print(f"📋 Available proteases ({len(proteases)} total):")
        for i, protease in enumerate(sorted(proteases), 1):
            print(f"  {i:3d}. {protease}")
    else:
        print("❌ Models directory not found!")

def list_data_files():
    """List available FASTA files"""
    if DATA_DIR.exists():
        fasta_files = list(DATA_DIR.glob("*.fasta")) + list(DATA_DIR.glob("*.fa"))
        print(f"📁 Available FASTA files ({len(fasta_files)} total):")
        for i, file in enumerate(fasta_files, 1):
            print(f"  {i:3d}. {file.name}")
    else:
        print("❌ Data directory not found!")

def read_proteases_file(proteases_file):
    """Read proteases from a file (one per line)"""
    file_path = DATA_DIR / proteases_file if not os.path.isabs(proteases_file) else Path(proteases_file)
    
    if not file_path.exists():
        print(f"❌ Proteases file not found: {file_path}")
        return []
    
    proteases = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip empty lines and comments
                proteases.append(line)
    
    print(f"📋 Loaded {len(proteases)} proteases from {file_path.name}")
    for i, protease in enumerate(proteases, 1):
        print(f"  {i:2d}. {protease}")
    
    return proteases

def create_display_format(results_file):
    """
    Create website-display format (what you see on the website UI)
    """
    df = pd.read_csv(results_file)
    
    # Create display format DataFrame
    display_df = pd.DataFrame()
    
    # Map to display column names (UI format)
    display_df['Protease'] = df['protease'] if 'protease' in df.columns else df.get('Protease', '')
    display_df['Sequence Id'] = df['sequence_id'] if 'sequence_id' in df.columns else df.get('Sequence Id', '')
    display_df['Position'] = df['position'] if 'position' in df.columns else df.get('Position', '')
    display_df['Cleavage site'] = df['seqs'] if 'seqs' in df.columns else df.get('Cleavage site', '')
    display_df['Prediction Score'] = df['pro'] if 'pro' in df.columns else df.get('Prediction Score', '')
    
    # Add ranking by prediction score (descending)
    display_df = display_df.sort_values('Prediction Score', ascending=False).reset_index(drop=True)
    display_df.insert(0, 'Rank', range(1, len(display_df) + 1))
    
    return display_df

def run_batch_prediction(fasta_file, proteases_list, output_base, output_format):
    """
    Run predictions for multiple proteases and organize results
    """
    # Validate inputs
    fasta_path = DATA_DIR / fasta_file if not os.path.isabs(fasta_file) else Path(fasta_file)
    if not fasta_path.exists():
        print(f"❌ FASTA file not found: {fasta_path}")
        return False
    
    # Create batch output directory
    batch_dir = RESULTS_DIR / output_base
    batch_dir.mkdir(parents=True, exist_ok=True)
    
    # Create plots directory
    batch_plots_dir = PLOTS_DIR / output_base
    batch_plots_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"🚀 Starting batch prediction for {len(proteases_list)} proteases...")
    print(f"📁 Results will be saved to: {batch_dir}")
    print(f"📊 Plots will be saved to: {batch_plots_dir}")
    
    all_results = []
    successful_predictions = 0
    
    for i, protease in enumerate(proteases_list, 1):
        print(f"\n{'='*60}")
        print(f"🔬 Processing {i}/{len(proteases_list)}: {protease}")
        print(f"{'='*60}")
        
        try:
            # Create individual protease directory
            protease_dir = batch_dir / protease
            protease_dir.mkdir(exist_ok=True)
            
            # Build command
            cmd = [
                'python', str(CODE_ROOT / 'Prosperousplus.py'),
                '--predictfile', str(fasta_path),
                '--outputpath', str(protease_dir),
                '--inputType', 'fasta',
                '--protease', protease,
                '--mode', 'prediction',
                '--processNum', '2',
                '--PLOT', 'Yes'
            ]
            
            print(f"Running: {' '.join(cmd)}")
            
            # Run prediction
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=CODE_ROOT)
            
            if result.returncode == 0:
                print(f"✅ {protease} prediction completed successfully!")
                
                # Load and process results
                results_file = protease_dir / 'results.csv'
                if results_file.exists():
                    df = pd.read_csv(results_file)
                    all_results.append(df)
                    
                    # Create individual plots using our plotting script
                    individual_plot_dir = batch_plots_dir / protease
                    individual_plot_dir.mkdir(exist_ok=True)
                    
                    plot_cmd = [
                        'python', str(SCRIPTS_DIR / 'plot_predictions_simple.py'),
                        '--results', str(results_file),
                        '--output', str(individual_plot_dir)
                    ]
                    
                    plot_result = subprocess.run(plot_cmd, capture_output=True, text=True)
                    if plot_result.returncode == 0:
                        print(f"📊 Plot generated for {protease}")
                    else:
                        print(f"⚠️ Plot generation failed for {protease}")
                    
                    successful_predictions += 1
                else:
                    print(f"⚠️ Results file not found for {protease}")
            else:
                print(f"❌ {protease} prediction failed!")
                print(f"Error: {result.stderr}")
                
        except Exception as e:
            print(f"❌ Error processing {protease}: {str(e)}")
            continue
    
    # Combine all results
    if all_results:
        print(f"\n{'='*60}")
        print(f"📋 Combining results from {successful_predictions} successful predictions...")
        
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Save combined results in requested format
        if output_format == 'display':
            display_df = create_display_format_from_combined(combined_df)
            combined_file = batch_dir / 'combined_results_display.csv'
            display_df.to_csv(combined_file, index=False)
            print(f"✅ Combined results saved: {combined_file}")
        else:
            combined_file = batch_dir / 'combined_results.csv'
            combined_df.to_csv(combined_file, index=False)
            print(f"✅ Combined results saved: {combined_file}")
        
        # Generate summary statistics
        create_batch_summary(combined_df, batch_dir, proteases_list)
        
        # Create comparison plots
        create_comparison_plots(combined_df, batch_plots_dir)
        
        print(f"\n🎉 Batch analysis complete!")
        print(f"📊 {successful_predictions}/{len(proteases_list)} proteases processed successfully")
        print(f"📁 Results: {batch_dir}")
        print(f"📈 Plots: {batch_plots_dir}")
        
        return True
    else:
        print(f"\n❌ No successful predictions to combine!")
        return False

def create_display_format_from_combined(combined_df):
    """Create display format from combined results"""
    display_df = pd.DataFrame()
    display_df['Protease'] = combined_df['protease']
    display_df['Sequence Id'] = combined_df['sequence_id']
    display_df['Position'] = combined_df['position']
    display_df['Cleavage site'] = combined_df['seqs']
    display_df['Prediction Score'] = combined_df['pro']
    
    # Add ranking by prediction score (descending)
    display_df = display_df.sort_values('Prediction Score', ascending=False).reset_index(drop=True)
    display_df.insert(0, 'Rank', range(1, len(display_df) + 1))
    
    return display_df

def create_batch_summary(combined_df, output_dir, proteases_list):
    """Create summary statistics for batch analysis"""
    summary_file = output_dir / 'batch_summary.txt'
    
    with open(summary_file, 'w') as f:
        f.write("PROSPEROUSPLUS BATCH ANALYSIS SUMMARY\n")
        f.write("="*50 + "\n\n")
        
        f.write(f"Proteases analyzed: {len(proteases_list)}\n")
        f.write(f"Proteases: {', '.join(proteases_list)}\n\n")
        
        f.write(f"Total predictions: {len(combined_df)}\n")
        f.write(f"Predicted cleavages: {sum(combined_df['prediction'] == 1)}\n")
        f.write(f"Average prediction score: {combined_df['pro'].mean():.3f}\n")
        f.write(f"Max prediction score: {combined_df['pro'].max():.3f}\n")
        f.write(f"Sites above 0.5: {sum(combined_df['pro'] > 0.5)}\n")
        f.write(f"Sites above 0.8: {sum(combined_df['pro'] > 0.8)}\n\n")
        
        # Top sites overall
        f.write("TOP 20 CLEAVAGE SITES (ALL PROTEASES):\n")
        f.write("-" * 40 + "\n")
        top_sites = combined_df.nlargest(20, 'pro')
        for i, (_, row) in enumerate(top_sites.iterrows(), 1):
            f.write(f"{i:2d}. {row['protease']} pos {row['position']} ({row['seqs']}) = {row['pro']:.3f}\n")
        
        # Per protease summary
        f.write("\nPER PROTEASE SUMMARY:\n")
        f.write("-" * 40 + "\n")
        for protease in proteases_list:
            protease_data = combined_df[combined_df['protease'] == protease]
            if len(protease_data) > 0:
                f.write(f"{protease}:\n")
                f.write(f"  Predictions: {len(protease_data)}\n")
                f.write(f"  Cleavages: {sum(protease_data['prediction'] == 1)}\n")
                f.write(f"  Avg score: {protease_data['pro'].mean():.3f}\n")
                f.write(f"  Max score: {protease_data['pro'].max():.3f}\n")
                f.write(f"  Top site: pos {protease_data.loc[protease_data['pro'].idxmax(), 'position']} = {protease_data['pro'].max():.3f}\n\n")
    
    print(f"📋 Summary saved: {summary_file}")

def create_comparison_plots(combined_df, plots_dir):
    """Create comparison plots across proteases"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    # Plot comparing max scores per protease
    protease_stats = combined_df.groupby('protease')['pro'].agg(['max', 'mean', 'count']).reset_index()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Max scores comparison
    ax1.bar(protease_stats['protease'], protease_stats['max'], color='skyblue', alpha=0.7)
    ax1.set_title('Maximum Prediction Score by Protease', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Protease', fontsize=12)
    ax1.set_ylabel('Max Prediction Score', fontsize=12)
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3)
    
    # Average scores comparison
    ax2.bar(protease_stats['protease'], protease_stats['mean'], color='lightcoral', alpha=0.7)
    ax2.set_title('Average Prediction Score by Protease', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Protease', fontsize=12)
    ax2.set_ylabel('Average Prediction Score', fontsize=12)
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    comparison_file = plots_dir / 'protease_comparison.png'
    plt.savefig(comparison_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"📊 Comparison plot saved: {comparison_file}")

def run_single_prediction(fasta_file, protease, output_name, output_format):
    """Run prediction for a single protease (existing functionality)"""
    # Validate inputs
    fasta_path = DATA_DIR / fasta_file if not os.path.isabs(fasta_file) else Path(fasta_file)
    if not fasta_path.exists():
        print(f"❌ FASTA file not found: {fasta_path}")
        return False
    
    # Create output directory
    output_dir = RESULTS_DIR / output_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Running prediction...")
    print(f"FASTA file: {fasta_path}")
    print(f"Protease: {protease}")
    print(f"Output directory: {output_dir}")
    
    # Build command
    cmd = [
        'python', str(CODE_ROOT / 'Prosperousplus.py'),
        '--predictfile', str(fasta_path),
        '--outputpath', str(output_dir),
        '--inputType', 'fasta',
        '--protease', protease,
        '--mode', 'prediction',
        '--processNum', '2',
        '--PLOT', 'Yes'
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    # Run prediction
    result = subprocess.run(cmd, cwd=CODE_ROOT)
    
    if result.returncode == 0:
        results_file = output_dir / 'results.csv'
        
        # Convert format if requested
        if output_format == 'display':
            display_df = create_display_format(results_file)
            display_file = output_dir / 'results_display_format.csv'
            display_df.to_csv(display_file, index=False)
            print(f"✅ Display format saved: {display_file}")
        
        print("✅ Prediction completed successfully!")
        print(f"📊 Results saved to: {results_file}")
        print(f"📋 Results in {output_format} format {'(matches website export)' if output_format == 'raw' else '(matches website display)'}")
        
        # Generate plots
        print(f"\n📈 Generating plots...")
        plot_output_dir = PLOTS_DIR / output_name
        plot_cmd = [
            'python', str(SCRIPTS_DIR / 'plot_predictions_simple.py'),
            '--results', str(results_file),
            '--output', str(plot_output_dir)
        ]
        
        plot_result = subprocess.run(plot_cmd, cwd=SCRIPTS_DIR)
        
        if plot_result.returncode == 0:
            print(f"📊 Plots saved to: {plot_output_dir}")
            
            print(f"\n🎉 Analysis complete!")
            print(f"📁 Project structure:")
            print(f"   Data: {DATA_DIR}")
            print(f"   Results: {RESULTS_DIR}")
            print(f"   Plots: {PLOTS_DIR}")
            
            print(f"\n📋 Format Info:")
            if output_format == 'raw':
                print(f"   ✅ Raw format (matches website export/download)")
                print(f"   📁 Columns: protease, sequence_id, position, seqs, prediction, pro")
            else:
                print(f"   ✅ Display format (matches website UI)")
                print(f"   📁 Columns: Rank, Protease, Sequence Id, Position, Cleavage site, Prediction Score")
            
            return True
        else:
            print(f"⚠️ Plot generation failed")
            return False
    else:
        print("❌ Prediction failed!")
        return False

def main():
    parser = argparse.ArgumentParser(description='ProsperousPlus Prediction Wrapper')
    parser.add_argument('--fasta', required=True, help='FASTA file name (in data directory)')
    parser.add_argument('--protease', help='Single protease to test (e.g., A01.001)')
    parser.add_argument('--proteases-file', help='File containing list of proteases (one per line)')
    parser.add_argument('--output', required=True, help='Output directory name')
    parser.add_argument('--format', choices=['raw', 'display'], default='raw',
                       help='Output format: raw (matches website export) or display (matches website UI)')
    parser.add_argument('--list-data', action='store_true', help='List available FASTA files')
    parser.add_argument('--list-proteases', action='store_true', help='List available proteases')
    
    args = parser.parse_args()
    
    # Handle list commands
    if args.list_data:
        list_data_files()
        return
    
    if args.list_proteases:
        list_available_proteases()
        return
    
    # Validate required arguments
    if not args.protease and not args.proteases_file:
        print("❌ Error: Must specify either --protease or --proteases-file")
        return
    
    if args.protease and args.proteases_file:
        print("❌ Error: Cannot specify both --protease and --proteases-file")
        return
    
    # Run batch or single prediction
    if args.proteases_file:
        # Batch processing
        proteases_list = read_proteases_file(args.proteases_file)
        if proteases_list:
            run_batch_prediction(args.fasta, proteases_list, args.output, args.format)
        else:
            print("❌ No proteases loaded from file")
    else:
        # Single protease
        run_single_prediction(args.fasta, args.protease, args.output, args.format)

if __name__ == "__main__":
    main() 