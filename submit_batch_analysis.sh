#!/bin/bash

# Helper script to submit batch AA analysis jobs
# This script submits individual jobs for each sample in sample_config.csv

CONFIG_FILE="sample_config.csv"
BASE_DIR="/broad/liulabdata/Nick_Krasnow/loopseq/pilot"
QSUB_OUT_DIR="$BASE_DIR/qsub_out"

# Create output directory for qsub logs
mkdir -p "$QSUB_OUT_DIR"

echo "Submitting batch AA analysis jobs..."
echo "Config file: $CONFIG_FILE"
echo "Base directory: $BASE_DIR"
echo ""

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file $CONFIG_FILE not found!" >&2
    exit 1
fi

# Process each sample in the config file
# Skip header line and process each sample
tail -n +2 "$CONFIG_FILE" | while IFS=',' read -r sample_name input_fasta reference_fasta gene_start gene_length; do
    # Skip empty lines
    if [ -z "$sample_name" ]; then
        continue
    fi
    
    echo "Submitting job for sample: $sample_name"
    
    # Create a temporary job script for this sample
    JOB_SCRIPT="temp_job_${sample_name}.sh"
    
    # Construct full paths
    INPUT_FASTA="$BASE_DIR/$input_fasta"
    REFERENCE_FASTA="$BASE_DIR/$reference_fasta"
    OUTPUT_DIR="$BASE_DIR/aa_analysis_results_${sample_name}_cluster_1based"
    
    # Create the job script
    cat > "$JOB_SCRIPT" << EOF
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
#$ -m ea
#$ -M nkrasnow@broadinstitute.org
#$ -j y

# Load modules
use Uger
reuse -q Anaconda3

# Set up conda environment
export PATH="/broad/software/conda/miniconda3/bin:$PATH"
source /broad/software/conda/miniconda3/etc/profile.d/conda.sh

# Create and activate conda environment if it doesn't exist
if ! conda env list | grep -q "loopseq"; then
    echo "Creating loopseq conda environment..."
    conda create -n loopseq python=3.9 -y
    conda activate loopseq
    conda install -c conda-forge pandas numpy matplotlib seaborn biopython -y
else
    echo "Activating existing loopseq conda environment..."
    conda activate loopseq
fi

# Verify Python and packages
python --version
python -c "import pandas, numpy, matplotlib, seaborn, Bio; print('All packages imported successfully')"

# Sample-specific parameters
SAMPLE_NAME="$sample_name"
INPUT_FASTA="$INPUT_FASTA"
REFERENCE_FASTA="$REFERENCE_FASTA"
GENE_START="$gene_start"
GENE_LENGTH="$gene_length"
MIN_LENGTH=1500
FREQUENCY_CUTOFF=0.01
OUTPUT_DIR="$OUTPUT_DIR"

# Create output directories
mkdir -p "\$OUTPUT_DIR"
mkdir -p "$QSUB_OUT_DIR"

# Set up log file names
LOG_FILE="$QSUB_OUT_DIR/aa_analysis_\${SAMPLE_NAME}_\$(date +%Y%m%d_%H%M%S).log"
ERROR_FILE="$QSUB_OUT_DIR/aa_analysis_\${SAMPLE_NAME}_\$(date +%Y%m%d_%H%M%S).err"

echo "Starting AA analysis for sample: \$SAMPLE_NAME"
echo "Log file: \$LOG_FILE"
echo "Error file: \$ERROR_FILE"
echo "Input FASTA: \$INPUT_FASTA"
echo "Reference FASTA: \$REFERENCE_FASTA"
echo "Gene start: \$GENE_START"
echo "Gene length: \$GENE_LENGTH"
echo "Output directory: \$OUTPUT_DIR"
echo ""

# Check if input files exist
if [ ! -f "\$INPUT_FASTA" ]; then
    echo "Error: Input FASTA not found: \$INPUT_FASTA" >&2
    exit 1
fi

if [ ! -f "\$REFERENCE_FASTA" ]; then
    echo "Error: Reference FASTA not found: \$REFERENCE_FASTA" >&2
    exit 1
fi

# Run analysis
echo "Running analysis for \$SAMPLE_NAME..."
python "$BASE_DIR/analyze_aa_mutations.py" \\
    --input_fasta "\$INPUT_FASTA" \\
    --reference_fasta "\$REFERENCE_FASTA" \\
    --min_length "\$MIN_LENGTH" \\
    --gene_start "\$GENE_START" \\
    --gene_length "\$GENE_LENGTH" \\
    --frequency_cutoff "\$FREQUENCY_CUTOFF" \\
    --output_dir "\$OUTPUT_DIR" \\
    > "\$LOG_FILE" 2> "\$ERROR_FILE"

# Check if analysis completed successfully
if [ \$? -eq 0 ]; then
    echo "Analysis completed successfully for \$SAMPLE_NAME"
else
    echo "Analysis failed for \$SAMPLE_NAME (check \$ERROR_FILE)"
    exit 1
fi

echo "Results saved to \$OUTPUT_DIR"
echo "Logs saved to \$LOG_FILE and \$ERROR_FILE"
EOF

    # Make the job script executable
    chmod +x "$JOB_SCRIPT"
    
    # Submit the job
    JOB_ID=$(qsub "$JOB_SCRIPT")
    echo "  Job submitted with ID: $JOB_ID"
    
    # Clean up the temporary job script
    rm "$JOB_SCRIPT"
    
    echo ""
done

echo "All batch analysis jobs submitted!"
echo "Monitor job status with: qstat"
echo "Check logs in: $QSUB_OUT_DIR"
