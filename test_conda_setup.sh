#!/bin/bash

# Test script to verify conda environment setup
echo "Testing conda environment setup..."

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
echo "Python version:"
python --version

echo "Testing package imports:"
python -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
print('All packages imported successfully!')
print('Pandas version:', pd.__version__)
print('NumPy version:', np.__version__)
print('Matplotlib version:', plt.matplotlib.__version__)
print('Seaborn version:', sns.__version__)
print('Biopython version:', SeqIO.__version__)
"

echo "Conda environment test completed successfully!"
