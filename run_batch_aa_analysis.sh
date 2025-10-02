#!/bin/bash

# Batch AA analysis script
# Runs run_aa_analysis_local.sh for multiple samples with different configurations
# Usage: bash run_batch_aa_analysis.sh

# Configuration arrays - EDIT THESE FOR YOUR SAMPLES
# Each array should have the same number of elements, with corresponding indices

# Sample names (used for output directory naming)
SAMPLE_NAMES=(
    "A1"
    "A2" 
    "B3"
    "C3"
    "D3"
    "E3"
    "F3"
)

# Input FASTA file paths
INPUT_FASTAS=(
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_A1/L1e8e291d_sample_A1_contig_list_trimmed.fa"
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_A2/L1e8e291d_sample_A2_contig_list_trimmed.fa"
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_B3/L1e8e291d_sample_B3_contig_list_trimmed.fa"
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_C3/L1e8e291d_sample_C3_contig_list_trimmed.fa"
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_D3/L1e8e291d_sample_D3_contig_list_trimmed.fa"
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_E3/L1e8e291d_sample_E3_contig_list_trimmed.fa"
    "L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_F3/L1e8e291d_sample_F3_contig_list_trimmed.fa"
)

# Reference FASTA files
REFERENCE_FASTAS=(
    "sp23_dna.fa"
    "sp23_dna.fa"
    "sp24_dna.fa"
    "sp23_dna.fa"
    "sp24_dna.fa"
    "sp23_dna.fa"
    "sp24_dna.fa"
)

# Gene start sequences
GENE_STARTS=(
    "ATGCCCAAAATAAATACATTTAATT"
    "ATGCCCAAAATAAATACATTTAATT"
    "ATGCCCAAAATAAATACATTTAATT"
    "ATGCCCAAAATAAATACATTTAATT"
    "ATGCCCAAAATAAATACATTTAATT"
    "ATGCCCAAAATAAATACATTTAATT"
    "ATGCCCAAAATAAATACATTTAATT"
)

# Gene lengths
GENE_LENGTHS=(
    1233
    1233
    1233
    1233
    1233
    1233
    1233
)

# Common parameters
MIN_LENGTH=1500
FREQUENCY_CUTOFF=0.01

# Other sample directories for comparison (space-separated, optional)
# Leave empty to disable comparison
OTHER_SAMPLE_DIRS=""

# Check that all arrays have the same length
check_array_lengths() {
    local names_len=${#SAMPLE_NAMES[@]}
    local fastas_len=${#INPUT_FASTAS[@]}
    local refs_len=${#REFERENCE_FASTAS[@]}
    local starts_len=${#GENE_STARTS[@]}
    local lengths_len=${#GENE_LENGTHS[@]}
    
    if [ $names_len -ne $fastas_len ] || [ $names_len -ne $refs_len ] || [ $names_len -ne $starts_len ] || [ $names_len -ne $lengths_len ]; then
        echo "Error: All configuration arrays must have the same length"
        echo "SAMPLE_NAMES: $names_len"
        echo "INPUT_FASTAS: $fastas_len"
        echo "REFERENCE_FASTAS: $refs_len"
        echo "GENE_STARTS: $starts_len"
        echo "GENE_LENGTHS: $lengths_len"
        exit 1
    fi
}

# Function to run analysis for a single sample
run_single_analysis() {
    local sample_name="$1"
    local input_fasta="$2"
    local reference_fasta="$3"
    local gene_start="$4"
    local gene_length="$5"
    local output_dir="aa_analysis_results_${sample_name}_batch"
    
    echo "=========================================="
    echo "Processing sample: $sample_name"
    echo "Input FASTA: $input_fasta"
    echo "Reference FASTA: $reference_fasta"
    echo "Gene start: $gene_start"
    echo "Gene length: $gene_length"
    echo "Output directory: $output_dir"
    echo "=========================================="
    
    # Check if input file exists
    if [ ! -f "$input_fasta" ]; then
        echo "Warning: Input FASTA file not found: $input_fasta"
        echo "Skipping sample $sample_name"
        return 1
    fi
    
    # Check if reference file exists
    if [ ! -f "$reference_fasta" ]; then
        echo "Warning: Reference FASTA file not found: $reference_fasta"
        echo "Skipping sample $sample_name"
        return 1
    fi
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Run the analysis
    python analyze_aa_mutations.py \
        --input_fasta "$input_fasta" \
        --reference_fasta "$reference_fasta" \
        --min_length "$MIN_LENGTH" \
        --gene_start "$gene_start" \
        --gene_length "$gene_length" \
        --frequency_cutoff "$FREQUENCY_CUTOFF" \
        --output_dir "$output_dir" \
        --other_sample_dirs $OTHER_SAMPLE_DIRS
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "✓ Successfully completed analysis for $sample_name"
    else
        echo "✗ Analysis failed for $sample_name (exit code: $exit_code)"
    fi
    
    echo ""
    return $exit_code
}

# Main execution
main() {
    echo "Starting batch AA analysis..."
    echo "Number of samples: ${#SAMPLE_NAMES[@]}"
    echo "Other sample directories: ${OTHER_SAMPLE_DIRS:-'None'}"
    echo ""
    
    # Check array lengths
    check_array_lengths
    
    # Track results
    local success_count=0
    local total_count=${#SAMPLE_NAMES[@]}
    local failed_samples=()
    
    # Process each sample
    for i in "${!SAMPLE_NAMES[@]}"; do
        local sample_name="${SAMPLE_NAMES[$i]}"
        local input_fasta="${INPUT_FASTAS[$i]}"
        local reference_fasta="${REFERENCE_FASTAS[$i]}"
        local gene_start="${GENE_STARTS[$i]}"
        local gene_length="${GENE_LENGTHS[$i]}"
        
        if run_single_analysis "$sample_name" "$input_fasta" "$reference_fasta" "$gene_start" "$gene_length"; then
            ((success_count++))
        else
            failed_samples+=("$sample_name")
        fi
    done
    
    # Summary
    echo "=========================================="
    echo "BATCH ANALYSIS COMPLETE"
    echo "=========================================="
    echo "Total samples: $total_count"
    echo "Successful: $success_count"
    echo "Failed: $((total_count - success_count))"
    
    if [ ${#failed_samples[@]} -gt 0 ]; then
        echo "Failed samples: ${failed_samples[*]}"
    fi
    
    echo ""
    echo "Results saved in directories:"
    for sample_name in "${SAMPLE_NAMES[@]}"; do
        echo "  - aa_analysis_results_${sample_name}_batch/"
    done
}

# Run main function
main "$@"
