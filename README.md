# ProsperousPlus Project Organization

This directory contains an organized workspace for ProsperousPlus protease cleavage predictions, separate from the main git repository. All results and plots are automatically generated in the official ProsperousPlus website format.

## Directory Structure

```
/shared/macdata/groups/ppc/projects/ProsperousPlus/
├── data/           # Input FASTA files and protease lists
├── results/        # Prediction results (CSV files in website format)
├── plots/          # Generated visualization plots (website style)
├── scripts/        # Analysis and plotting scripts
└── README.md       # This file
```

## Quick Start

### 1. List Available Resources
```bash
cd /shared/macdata/groups/ppc/projects/ProsperousPlus/scripts

# List available FASTA files
python run_prediction.py --list-data

# List available protease models
python run_prediction.py --list-proteases
```

### 2. Single Protease Analysis
```bash
# Basic prediction with automatic plotting (website format)
python run_prediction.py --fasta your_proteins.fasta --protease A01.001 --output analysis_name

# With display format (matches website UI)
python run_prediction.py --fasta proteins.fasta --protease A01.001 --output analysis_name --format display
```

### 3. Batch Analysis (Multiple Proteases)
```bash
# Create a proteases list file
echo "C01.060
C01.032  
C01.034
C14.003
S01.001" > ../data/my_proteases.txt

# Run batch analysis for all proteases
python run_prediction.py --fasta proteins.fasta --proteases-file my_proteases.txt --output batch_analysis

# Batch with display format
python run_prediction.py --fasta proteins.fasta --proteases-file proteases_list.txt --output batch_analysis --format display
```

### 4. Generate Custom Plots
```bash
# Generate plots from existing results (website style)
python plot_predictions_simple.py --results ../results/analysis/results.csv --output ../plots/custom_plots
```

## Batch Processing Features

### Multi-Protease Analysis
The batch processing system allows you to test multiple proteases against your protein sequences efficiently:

**Benefits:**
- **Automated processing**: Run multiple proteases with a single command
- **Organized results**: Each protease gets its own directory
- **Combined analysis**: All results merged into comprehensive summaries
- **Comparison plots**: Automatic generation of protease comparison charts
- **Individual plots**: Separate visualization for each protease

### Batch Directory Structure
```
results/batch_analysis/
├── C01.060/
│   ├── results.csv
│   └── results/        # Original ProsperousPlus output
├── C01.032/
│   ├── results.csv
│   └── results/
├── combined_results.csv        # All proteases combined
└── batch_summary.txt          # Statistical summary

plots/batch_analysis/
├── C01.060/
│   └── protein_cleavage_prediction.png
├── C01.032/
│   └── protein_cleavage_prediction.png
└── protease_comparison.png    # Compare all proteases
```

### Proteases List File Format
Create a simple text file with one protease per line:
```
# My proteases of interest
C01.060
C01.032
C01.034
C01.036
C14.003
S01.001
# Lines starting with # are ignored
```

### Batch Analysis Output

**Combined Results CSV**: Contains predictions from all proteases, ranked by prediction score across the entire dataset.

**Batch Summary Report**: Includes:
- Total predictions across all proteases
- Per-protease statistics (max scores, cleavage counts, etc.)
- Top 20 cleavage sites overall
- Performance comparison between proteases

**Comparison Plots**: Visual comparison showing:
- Maximum prediction scores by protease
- Average prediction scores by protease
- Side-by-side performance metrics

## Understanding the Results

### Results CSV Format (Website Standard)
The main results file (`results.csv`) contains these columns exactly as shown on the ProsperousPlus website:

| Column | Description |
|--------|-------------|
| `protease` | Protease model used (e.g., A01.001) |
| `sequence_id` | Protein identifier |
| `position` | Residue position in the protein |
| `seqs` | 8-amino acid window around cleavage site |
| `prediction` | Binary prediction (1=cleave, 0=no cleave) |
| `pro` | Cleavage probability score (0-1) |

### Display Format (Website UI Style)
When using `--format display`, results include ranking and formatted column names:

| Column | Description |
|--------|-------------|
| `Rank` | Ranking by prediction score (1 = highest) |
| `Protease` | Protease model used |
| `Sequence Id` | Protein identifier |
| `Position` | Residue position in the protein |
| `Cleavage site` | 8-amino acid window around cleavage site |
| `Prediction Score` | Cleavage probability score (0-1) |

### Generated Plots (Website Style)
All plots match the official ProsperousPlus website styling:

1. **Individual protein plots**: `protein_cleavage_prediction.png`
   - **Title**: "ProsperousPlus Prediction for [SequenceID]"
   - **Y-axis**: "Prediction Score" (0.0 to 1.0)
   - **X-axis**: "Residue Position"
   - **Style**: Black line, green dashed threshold line at 0.5, white background
   
2. **Comparison plots**: Visual comparison of protease performance

## Input File Requirements

### FASTA Format
Your input files should be in standard FASTA format and can contain multiple sequences:
```
>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEE...
>protein2  
MGSSHHHHHHSSGLVPRGSHM...
>protein3
MASNDYTQQATQSYGAYPTQPGQGYSQQSS...
```

### Multi-Sequence Processing
- **Single proteases**: Each sequence gets analyzed separately
- **Batch analysis**: All proteases test all sequences
- **Combined results**: Full matrix of protease × sequence × position predictions
- **Individual plots**: Separate visualization for each sequence

### Supported Sequences
- Any protein sequence length
- Standard 20 amino acids
- Non-standard amino acids will be replaced with '-'

## Available Protease Models

ProsperousPlus includes 110+ pre-trained protease models including:
- **A01.001**: Pepsin A
- **C01.060**: Caspase-3
- **C01.032**: Caspase-1
- **M24.026**: Methionine aminopeptidase
- **S01.001**: Chymotrypsin
- **S01.140**: Trypsin
- And many more...

Use `--list-proteases` to see the complete list.

## Analysis Examples

### Single Protease Analysis
```bash
# Analyze TDP-43 with caspase-3
python run_prediction.py --fasta TDP43.fasta --protease C01.060 --output TDP43_caspase3
```

### Multi-Protease Batch Analysis
```bash
# Test multiple caspases
echo "C01.060
C01.032
C01.034
C14.003" > ../data/caspases.txt

python run_prediction.py --fasta proteins.fasta --proteases-file caspases.txt --output caspase_screen
```

### Large-Scale Screening
```bash
# Test many proteases against multiple proteins
python run_prediction.py --fasta multi_proteins.fasta --proteases-file proteases_list.txt --output comprehensive_screen
```

## Analysis Tips

### High-Confidence Predictions
- **Prediction scores > 0.8** are typically high-confidence
- **Sites above 0.5** are considered positive predictions (threshold line)
- Consider the biological context of predicted cleavage sites
- Multiple predictions at similar positions may indicate true cleavage sites

### Batch Analysis Strategy
1. **Start small**: Test a few proteases first to verify results
2. **Use biological relevance**: Select proteases relevant to your system
3. **Check individual results**: Verify individual protease predictions before interpreting combined results
4. **Compare systematically**: Use the comparison plots to identify most active proteases

### Performance Considerations
- **Large batches**: 9 proteases × 2 proteins = ~18 individual predictions
- **Processing time**: Each protease takes a few minutes per protein
- **Storage**: Results scale with number of proteases × sequences × sequence length

## Example Output Summary

### Single Analysis
```
TOP 10 CLEAVAGE SITES (sorted by probability):
====================================================================================================
Rank   Protease   Sequence Id     Position   Cleavage site   Prediction Score
====================================================================================================
1      C01.060    protein1        32         KGFGFVRF        0.982          
2      C01.060    protein1        114        RAFAFVTF        0.976          

SUMMARY STATISTICS:
   Total predictions: 291
   Predicted cleavages: 92
   Average prediction score: 0.368
   Max prediction score: 0.982
   Sites above 0.5: 92
   Sites above 0.8: 27
```

### Batch Analysis Summary
```
PROSPEROUSPLUS BATCH ANALYSIS SUMMARY
==================================================

Proteases analyzed: 9
Proteases: C01.060, C01.032, C01.034, C01.036, C01.009, C14.003, C14.004, C02.001, C02.002

Total predictions: 2619
Predicted cleavages: 876
Average prediction score: 0.334
Max prediction score: 0.982
Sites above 0.5: 876
Sites above 0.8: 156

TOP 20 CLEAVAGE SITES (ALL PROTEASES):
----------------------------------------
 1. C01.060 pos 32 (KGFGFVRF) = 0.982
 2. C01.032 pos 114 (RAFAFVTF) = 0.976
...

PER PROTEASE SUMMARY:
----------------------------------------
C01.060:
  Predictions: 291
  Cleavages: 92
  Avg score: 0.368
  Max score: 0.982
  Top site: pos 32 = 0.982
```

## Troubleshooting

### Common Issues
1. **"FASTA file not found"**: Make sure your file is in the `data/` directory
2. **"Protease not found"**: Use `--list-proteases` to see available models
3. **"Proteases file not found"**: Check the path to your proteases list file
4. **Python environment**: Make sure you're using the `prosperousplus` conda environment

### Batch Processing Issues
1. **Partial completion**: Check individual protease directories for error logs
2. **Memory usage**: Large batches may require more system resources
3. **Disk space**: Ensure sufficient space for all results and plots

### Performance Tips
- Start with small test runs to verify setup
- Monitor system resources during large batch jobs
- Consider splitting very large analyses into smaller batches

## Advanced Usage

### Custom Protease Lists
Create targeted lists for specific biological questions:
```bash
# Apoptosis-related proteases
echo "C01.060
C01.032
C14.003" > ../data/apoptosis_proteases.txt

# Digestive proteases  
echo "A01.001
S01.001
S01.140" > ../data/digestive_proteases.txt
```

### Format Conversion
Convert between raw and display formats:
```bash
# Raw format (default) - matches website downloads
python run_prediction.py --fasta proteins.fasta --protease A01.001 --output analysis --format raw

# Display format - matches website UI
python run_prediction.py --fasta proteins.fasta --protease A01.001 --output analysis --format display
```

## Website Format Features

### Automatic Standardization
- **All results** are automatically in website-standard format
- **All plots** use the exact styling from the ProsperousPlus website
- **Consistent formatting** across all analyses

### Visual Consistency
- Title format: "ProsperousPlus Prediction for [SequenceID]"
- Axis labels: "Prediction Score" (Y) and "Residue Position" (X)
- Green dashed threshold line at 0.5
- Clean white background matching website

## Related Files

- **Source code**: `/shared/macdata/groups/ppc/code/ProsperousPlus/`
- **Original README**: See the main repository for detailed algorithm information
- **Model files**: Pre-trained models are in the source code directory

## Getting Help

For questions about:
- **Algorithm details**: Check the original ProsperousPlus paper and repository
- **This organized setup**: Ask the person who set this up for you
- **Biological interpretation**: Consult domain experts
- **Website format**: Compare with official ProsperousPlus website results

---

*This organized workspace was created to keep your data and results separate from the git repository for better project management. All outputs match the official ProsperousPlus website format for consistency and standardization.* 