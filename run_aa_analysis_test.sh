#!/bin/bash

# Test version of AA analysis script with downsampled data
# Run with: bash run_aa_analysis_test.sh

# Activate conda environment
source activate loopseq || conda activate loopseq || echo "Warning: Could not activate loopseq environment"

# Install required packages if not present
pip install biopython --quiet

# Parameters - using downsampled data
INPUT_FASTA="L1e8e291d_sample_A1_downsampled_10x.fa"
REFERENCE_FASTA="sp23_dna.fa"
MIN_LENGTH=1500
GENE_START="ATGCCCAAAATAAATACATTTAATT"
GENE_LENGTH=1233
FREQUENCY_CUTOFF=0.01
OUTPUT_DIR="aa_analysis_results_A1_test_10x_v2"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting test AA analysis with 10x downsampled data..."
echo "Input FASTA: $INPUT_FASTA"
echo "Reference FASTA: $REFERENCE_FASTA"
echo "Min length: $MIN_LENGTH"
echo "Gene start: $GENE_START"
echo "Gene length: $GENE_LENGTH"
echo "Frequency cutoff: $FREQUENCY_CUTOFF"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run analysis
python analyze_aa_mutations.py \
    --input_fasta "$INPUT_FASTA" \
    --reference_fasta "$REFERENCE_FASTA" \
    --min_length "$MIN_LENGTH" \
    --gene_start "$GENE_START" \
    --gene_length "$GENE_LENGTH" \
    --frequency_cutoff "$FREQUENCY_CUTOFF" \
    --output_dir "$OUTPUT_DIR"

echo "Test analysis complete. Results saved to $OUTPUT_DIR"
