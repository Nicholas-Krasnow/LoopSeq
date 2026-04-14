#!/bin/bash

# Submit batch AA analysis jobs from sample_keys.xlsx (sheets "1-48" and "initial_set").
# Resolves paths relative to project root (parent of this script's directory).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
XLSX="${BASE_DIR}/sample_keys.xlsx"
ENUMERATOR="${SCRIPT_DIR}/enumerate_sample_keys.py"
ANALYZE_PY="${SCRIPT_DIR}/analyze_aa_mutations.py"
QSUB_OUT_DIR="${BASE_DIR}/qsub_out"

mkdir -p "$QSUB_OUT_DIR"

echo "Submitting batch AA analysis jobs..."
echo "Base directory: $BASE_DIR"
echo "Workbook: $XLSX"
echo ""

if [ ! -f "$XLSX" ]; then
    echo "Error: sample_keys.xlsx not found at $XLSX" >&2
    exit 1
fi

if [ ! -f "$ENUMERATOR" ]; then
    echo "Error: enumerator not found: $ENUMERATOR" >&2
    exit 1
fi

if [ ! -f "$ANALYZE_PY" ]; then
    echo "Error: analyze_aa_mutations.py not found: $ANALYZE_PY" >&2
    exit 1
fi

SAMPLES_TSV="$(mktemp)"
trap 'rm -f "$SAMPLES_TSV"' EXIT

if ! python3 "$ENUMERATOR" "$BASE_DIR" --xlsx "$XLSX" > "$SAMPLES_TSV"; then
    echo "Error: failed to enumerate samples from workbook" >&2
    exit 1
fi

if [ ! -s "$SAMPLES_TSV" ]; then
    echo "Error: no samples to submit (check workbook and directory layout)" >&2
    exit 1
fi

while IFS=$'\t' read -r sample_name input_fasta reference_fasta gene_length gene_start; do
    [ -z "$sample_name" ] && continue

    echo "Submitting job for sample: $sample_name"

    JOB_SCRIPT="$(mktemp "${TMPDIR:-/tmp}/aa_batch_job.XXXXXX.sh")"

    cat > "$JOB_SCRIPT" << EOF
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -l h_vmem=8G
#\$ -l h_rt=12:00:00
#\$ -m ea
#\$ -M nkrasnow@broadinstitute.org
#\$ -j y

use UGER
source /broad/software/scripts/useuse
reuse -q Anaconda3

export PATH="/broad/software/conda/miniconda3/bin:\$PATH"
source /broad/software/conda/miniconda3/etc/profile.d/conda.sh

if ! conda env list | grep -q "loopseq"; then
    echo "Creating loopseq conda environment..."
    conda create -n loopseq python=3.9 -y
    conda activate loopseq
    conda install -c conda-forge pandas numpy matplotlib seaborn biopython openpyxl -y
else
    echo "Activating existing loopseq conda environment..."
    conda activate loopseq
fi

echo "Active python (after conda activate): \$(which python)"
python --version
echo "Loopseq python (via conda run): \$(conda run -n loopseq python -c 'import sys; print(sys.executable)')"
conda run -n loopseq python --version

# If env exists but is missing dependencies, fix in-place.
conda run -n loopseq python - <<'PY'
import importlib, sys
mods = ["pandas","numpy","matplotlib","seaborn","Bio","openpyxl"]
missing = [m for m in mods if importlib.util.find_spec(m) is None]
if missing:
    print("Missing modules:", ", ".join(missing))
    sys.exit(2)
print("All packages importable")
PY
IMPORT_OK=\$?
if [ "\$IMPORT_OK" -eq 2 ]; then
    echo "Installing missing conda packages into loopseq env..."
    conda install -n loopseq -c conda-forge pandas numpy matplotlib seaborn biopython openpyxl -y
fi

SAMPLE_NAME="$sample_name"
INPUT_FASTA="$input_fasta"
REFERENCE_FASTA="$reference_fasta"
GENE_START="$gene_start"
GENE_LENGTH="$gene_length"
MIN_LENGTH=1500
FREQUENCY_CUTOFF=0.01
OUTPUT_DIR="$BASE_DIR/aa_analysis_results_\${SAMPLE_NAME}_cluster_1based"

mkdir -p "\$OUTPUT_DIR"
mkdir -p "$QSUB_OUT_DIR"

LOG_FILE="$QSUB_OUT_DIR/aa_analysis_\${SAMPLE_NAME}_\$(date +%Y%m%d_%H%M%S).log"
ERROR_FILE="$QSUB_OUT_DIR/aa_analysis_\${SAMPLE_NAME}_\$(date +%Y%m%d_%H%M%S).err"

echo "Starting AA analysis for sample: \$SAMPLE_NAME"
echo "Input FASTA: \$INPUT_FASTA"
echo "Reference FASTA: \$REFERENCE_FASTA"
echo "Gene length: \$GENE_LENGTH"
echo "Output directory: \$OUTPUT_DIR"

if [ ! -f "\$INPUT_FASTA" ]; then
    echo "Error: Input FASTA not found: \$INPUT_FASTA" >&2
    exit 1
fi

if [ ! -f "\$REFERENCE_FASTA" ]; then
    echo "Error: Reference FASTA not found: \$REFERENCE_FASTA" >&2
    exit 1
fi

echo "Running analysis for \$SAMPLE_NAME..."
conda run -n loopseq python "$ANALYZE_PY" \\
    --input_fasta "\$INPUT_FASTA" \\
    --reference_fasta "\$REFERENCE_FASTA" \\
    --min_length "\$MIN_LENGTH" \\
    --gene_start "\$GENE_START" \\
    --gene_length "\$GENE_LENGTH" \\
    --frequency_cutoff "\$FREQUENCY_CUTOFF" \\
    --output_dir "\$OUTPUT_DIR" \\
    > "\$LOG_FILE" 2> "\$ERROR_FILE"

if [ \$? -eq 0 ]; then
    echo "Analysis completed successfully for \$SAMPLE_NAME"
else
    echo "Analysis failed for \$SAMPLE_NAME (check \$ERROR_FILE)"
    exit 1
fi
EOF

    chmod +x "$JOB_SCRIPT"
    JOB_ID=$(qsub "$JOB_SCRIPT")
    echo "  Job submitted with ID: $JOB_ID"
    rm -f "$JOB_SCRIPT"
    echo ""
done < "$SAMPLES_TSV"

echo "All batch analysis jobs submitted!"
echo "Monitor job status with: qstat"
echo "Check logs in: $QSUB_OUT_DIR"
