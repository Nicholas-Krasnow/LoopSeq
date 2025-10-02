#!/bin/bash

# Local mutation set sharing analysis script
# This script analyzes mutation set sharing between samples with different reference fastas

# Set up conda environment
source activate loopseq

# Verify Python and packages
echo "Python version:"
python --version

echo "Testing package imports:"
python -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
print('All packages imported successfully!')
"

# Configuration
CONFIG_FILE="${1:-sample_config.csv}"
BASE_DIR="${2:-.}"
MIN_COUNT="${3:-100}"

echo "Starting mutation set sharing analysis..."
echo "Configuration file: $CONFIG_FILE"
echo "Base directory: $BASE_DIR"
echo "Minimum read count: $MIN_COUNT"
echo ""

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE"
    echo "Usage: bash run_mutation_sharing_analysis.sh [config_file.csv] [base_dir] [min_count]"
    exit 1
fi

# Run the analysis
python analyze_mutation_set_sharing.py \
    --config_file "$CONFIG_FILE" \
    --base_dir "$BASE_DIR" \
    --min_count "$MIN_COUNT"

echo ""
echo "Mutation set sharing analysis complete!"
echo "Check the _cluster_1based folders for:"
echo "  - mutation_set_sharing_histogram.png"
echo "  - mutation_set_sharing_histogram.pdf"
echo "  - mutation_set_sharing_analysis.csv"
