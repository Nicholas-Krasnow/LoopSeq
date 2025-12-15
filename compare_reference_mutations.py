#!/usr/bin/env python3
"""
Compare reference sequences and identify amino acid mutations.
This script compares sp23_dna.fa and sp24_dna.fa and identifies amino acid mutations
in the +1 reading frame.
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def translate_dna(sequence, reading_frame=1):
    """Translate DNA sequence to amino acid sequence using specified reading frame."""
    # Adjust for reading frame (1-based)
    start_pos = reading_frame - 1
    frame_seq = sequence[start_pos:]
    
    # Remove any non-DNA characters and convert to uppercase
    clean_seq = ''.join(c for c in frame_seq.upper() if c in 'ATCG')
    
    # Ensure length is multiple of 3
    if len(clean_seq) % 3 != 0:
        clean_seq = clean_seq[:-(len(clean_seq) % 3)]
    
    if len(clean_seq) == 0:
        return ""
    
    # Translate using Biopython
    seq_obj = Seq(clean_seq)
    protein = seq_obj.translate()
    return str(protein)

def find_aa_mutations(ref_aa, query_aa, reading_frame=1):
    """
    Find amino acid mutations between reference and query sequences.
    Returns list of (position, reference_aa, mutated_aa) tuples.
    """
    mutations = []
    min_len = min(len(ref_aa), len(query_aa))
    
    for i in range(min_len):
        if ref_aa[i] != query_aa[i]:
            # Convert amino acid position to 1-based
            aa_pos = i + 1
            mutations.append((aa_pos, ref_aa[i], query_aa[i]))
    
    return mutations

def main():
    # File paths
    sp23_file = "sp23_dna.fa"
    sp24_file = "sp24_dna.fa"
    output_file = "reference_mutations_sp24_vs_sp23.txt"
    
    print("Loading reference sequences...")
    
    # Check if files exist
    if not os.path.exists(sp23_file):
        print(f"Error: {sp23_file} not found!")
        return 1
    
    if not os.path.exists(sp24_file):
        print(f"Error: {sp24_file} not found!")
        return 1
    
    # Load sequences
    try:
        sp23_record = next(SeqIO.parse(sp23_file, "fasta"))
        sp24_record = next(SeqIO.parse(sp24_file, "fasta"))
    except Exception as e:
        print(f"Error loading sequences: {e}")
        return 1
    
    sp23_seq = str(sp23_record.seq).upper()
    sp24_seq = str(sp24_record.seq).upper()
    
    print(f"SP23 sequence length: {len(sp23_seq)}")
    print(f"SP24 sequence length: {len(sp24_seq)}")
    
    # Translate to amino acids using +1 reading frame
    print("Translating sequences to amino acids (+1 reading frame)...")
    sp23_aa = translate_dna(sp23_seq, reading_frame=1)
    sp24_aa = translate_dna(sp24_seq, reading_frame=1)
    
    print(f"SP23 protein length: {len(sp23_aa)}")
    print(f"SP24 protein length: {len(sp24_aa)}")
    
    # Find mutations (SP24 as reference, SP23 as query)
    print("Finding amino acid mutations...")
    mutations = find_aa_mutations(sp24_aa, sp23_aa, reading_frame=1)
    
    print(f"Found {len(mutations)} amino acid mutations")
    
    # Save results
    print(f"Saving results to {output_file}...")
    with open(output_file, 'w') as f:
        f.write("Amino Acid Mutations: SP24 vs SP23 (Reference: SP24)\n")
        f.write("=" * 60 + "\n")
        f.write("Position\tSP24_AA\tSP23_AA\tMutation\n")
        f.write("-" * 60 + "\n")
        
        for pos, ref_aa, mut_aa in mutations:
            mutation_str = f"{ref_aa}{pos}{mut_aa}"
            f.write(f"{pos}\t{ref_aa}\t{mut_aa}\t{mutation_str}\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Total mutations: {len(mutations)}\n")
        f.write(f"SP23 protein length: {len(sp23_aa)}\n")
        f.write(f"SP24 protein length: {len(sp24_aa)}\n")
        f.write(f"Reading frame: +1\n")
    
    print(f"Results saved to {output_file}")
    print(f"Total mutations found: {len(mutations)}")
    
    # Print summary to console
    if mutations:
        print("\nFirst 10 mutations:")
        for i, (pos, ref_aa, mut_aa) in enumerate(mutations[:10]):
            print(f"  {pos}: {ref_aa} -> {mut_aa}")
        if len(mutations) > 10:
            print(f"  ... and {len(mutations) - 10} more")
    else:
        print("No mutations found between the sequences.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
