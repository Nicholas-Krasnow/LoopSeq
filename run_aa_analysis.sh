#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
#$ -m ea
#$ -M nkrasnow@broadinstitute.org
#$ -j y

# AA analysis script for cluster submission
# Submit with: qsub run_aa_analysis.sh

# Load modules
useuse Uger
use Anaconda3

# Activate conda environment
source activate loopseq

# Install required packages if not present
pip install biopython --quiet
# EVcouplings removed - requires pre-computed model files

# Parameters - EDIT THESE PATHS
INPUT_FASTA="L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_A1/L1e8e291d_sample_A1_contig_list_trimmed.fa"
REFERENCE_FASTA="sp23_dna.fa"
MIN_LENGTH=1500
GENE_START="ATGCCCAAAATAAATACATTTAATT"
GENE_LENGTH=1233
FREQUENCY_CUTOFF=0.01
OUTPUT_DIR="aa_analysis_results_A1_cluster_1based"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting cluster AA analysis..."
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

echo "Analysis complete. Results saved to $OUTPUT_DIR"

