#!/bin/bash

# Batch AA analysis script using CSV configuration
# Usage: bash run_batch_aa_analysis_csv.sh [config_file.csv]
# Default config file: sample_config.csv

# Configuration
CONFIG_FILE="${1:-sample_config.csv}"
MIN_LENGTH=1500
FREQUENCY_CUTOFF=0.01
OTHER_SAMPLE_DIRS=""  # Space-separated list of directories for comparison

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE"
    echo "Usage: bash run_batch_aa_analysis_csv.sh [config_file.csv]"
    exit 1
fi

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
    echo "Starting batch AA analysis from CSV configuration..."
    echo "Configuration file: $CONFIG_FILE"
    echo "Other sample directories: ${OTHER_SAMPLE_DIRS:-'None'}"
    echo ""
    
    # Check if we have the required tools
    if ! command -v python &> /dev/null; then
        echo "Error: Python not found. Please ensure Python is available."
        exit 1
    fi
    
    # Read CSV and process each line (skip header)
    local line_count=0
    local success_count=0
    local failed_samples=()
    
    while IFS=',' read -r sample_name input_fasta reference_fasta gene_start gene_length; do
        # Skip header line
        if [ $line_count -eq 0 ]; then
            ((line_count++))
            continue
        fi
        
        # Skip empty lines
        if [ -z "$sample_name" ]; then
            continue
        fi
        
        ((line_count++))
        
        if run_single_analysis "$sample_name" "$input_fasta" "$reference_fasta" "$gene_start" "$gene_length"; then
            ((success_count++))
        else
            failed_samples+=("$sample_name")
        fi
    done < "$CONFIG_FILE"
    
    local total_count=$((line_count - 1))  # Subtract header
    
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
    # List output directories
    for sample_name in "${failed_samples[@]}"; do
        echo "  - aa_analysis_results_${sample_name}_batch/"
    done
    # Also show successful ones
    echo "  - (check for aa_analysis_results_*_batch/ directories)"
}

# Run main function
main "$@"
