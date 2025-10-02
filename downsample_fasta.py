#!/usr/bin/env python3
"""
Downsample FASTA file by randomly selecting 1/10th of the reads.
"""

import sys
import random
from Bio import SeqIO

def downsample_fasta(input_file, output_file, fraction=0.1, seed=42):
    """
    Downsample a FASTA file by randomly selecting a fraction of reads.
    
    Args:
        input_file: Input FASTA file path
        output_file: Output FASTA file path
        fraction: Fraction of reads to keep (default 0.1 for 10x downsampling)
        seed: Random seed for reproducibility
    """
    # Set random seed for reproducibility
    random.seed(seed)
    
    # Read all sequences
    print(f"Reading sequences from {input_file}...")
    sequences = list(SeqIO.parse(input_file, "fasta"))
    total_reads = len(sequences)
    print(f"Total reads: {total_reads}")
    
    # Calculate number of reads to keep
    reads_to_keep = int(total_reads * fraction)
    print(f"Keeping {reads_to_keep} reads ({fraction*100:.1f}% of total)")
    
    # Randomly sample reads
    sampled_sequences = random.sample(sequences, reads_to_keep)
    
    # Write downsampled sequences
    print(f"Writing downsampled sequences to {output_file}...")
    SeqIO.write(sampled_sequences, output_file, "fasta")
    
    print(f"Downsampling complete: {reads_to_keep}/{total_reads} reads kept")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python downsample_fasta.py <input_fasta> <output_fasta>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    downsample_fasta(input_file, output_file)
