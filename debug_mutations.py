#!/usr/bin/env python3
"""
Debug script to understand why no mutations are being found.
"""

from Bio import SeqIO
from Bio.Seq import Seq

def translate_dna(sequence: str) -> str:
    """Translate DNA sequence to amino acid sequence using standard genetic code."""
    # Remove any non-DNA characters and convert to uppercase
    clean_seq = ''.join(c for c in sequence.upper() if c in 'ATCG')
    
    # Ensure length is multiple of 3
    if len(clean_seq) % 3 != 0:
        clean_seq = clean_seq[:-(len(clean_seq) % 3)]
    
    if len(clean_seq) == 0:
        return ""
    
    # Translate using Biopython
    seq_obj = Seq(clean_seq)
    protein = seq_obj.translate()
    return str(protein)

def find_gene_start_with_alignment(query_seq: str, gene_start: str, max_mismatches: int = 2) -> tuple:
    """
    Find gene_start in query sequence with alignment, allowing for mismatches.
    Returns (aligned_query, start_pos, end_pos, is_reverse_complement).
    """
    query_upper = query_seq.upper()
    gene_start_upper = gene_start.upper()
    
    # Try forward strand
    best_mismatches = max_mismatches + 1
    best_start = -1
    best_aligned = ""
    
    for i in range(len(query_upper) - len(gene_start_upper) + 1):
        window = query_upper[i:i + len(gene_start_upper)]
        mismatches = sum(1 for a, b in zip(window, gene_start_upper) if a != b)
        
        if mismatches <= max_mismatches and mismatches < best_mismatches:
            best_mismatches = mismatches
            best_start = i
            best_aligned = window
    
    if best_start != -1:
        return best_aligned, best_start, best_start + len(gene_start_upper), False
    
    # Try reverse complement
    query_rc = str(Seq(query_seq).reverse_complement())
    query_rc_upper = query_rc.upper()
    
    for i in range(len(query_rc_upper) - len(gene_start_upper) + 1):
        window = query_rc_upper[i:i + len(gene_start_upper)]
        mismatches = sum(1 for a, b in zip(window, gene_start_upper) if a != b)
        
        if mismatches <= max_mismatches and mismatches < best_mismatches:
            best_mismatches = mismatches
            best_start = i
            best_aligned = window
    
    if best_start != -1:
        return best_aligned, best_start, best_start + len(gene_start_upper), True
    
    return "", -1, -1, False

def extract_subread_with_alignment(query_seq: str, gene_start: str, gene_length: int, max_mismatches: int = 2) -> tuple:
    """
    Extract subread using gene_start alignment.
    Returns (subread, mismatches, is_reverse_complement).
    """
    aligned_query, start_pos, end_pos, is_reverse = find_gene_start_with_alignment(query_seq, gene_start, max_mismatches)
    
    if start_pos == -1:
        return "", 0, False
    
    # Calculate mismatches
    mismatches = sum(1 for a, b in zip(aligned_query, gene_start.upper()) if a != b)
    
    # Extract subread: start at the beginning of gene_start, extend for gene_length
    subread_start = start_pos
    subread_end = subread_start + gene_length
    
    if is_reverse:
        # Extract from reverse complement
        query_rc = str(Seq(query_seq).reverse_complement())
        if subread_start >= 0 and subread_end <= len(query_rc):
            subread = query_rc[subread_start:subread_end]
        else:
            return "", mismatches, is_reverse
    else:
        # Extract from forward strand
        if subread_start >= 0 and subread_end <= len(query_seq):
            subread = query_seq[subread_start:subread_end]
        else:
            return "", mismatches, is_reverse
    
    return subread, mismatches, is_reverse

def compare_subread_to_reference(subread: str, reference_seq: str, min_identity: float = 0.9) -> tuple:
    """
    Compare subread directly to reference sequence without alignment.
    Returns (ref_aligned, subread_aligned, start, end, identity).
    """
    if not subread or len(subread) == 0:
        return "", "", 0, 0, 0.0
    
    # Find best matching position in reference
    best_identity = 0.0
    best_start = 0
    best_ref = ""
    best_subread = ""
    
    for start in range(len(reference_seq) - len(subread) + 1):
        ref_window = reference_seq[start:start + len(subread)]
        
        # Calculate identity
        matches = sum(1 for a, b in zip(subread, ref_window) if a == b)
        identity = matches / len(subread) if len(subread) > 0 else 0.0
        
        if identity > best_identity:
            best_identity = identity
            best_start = start
            best_ref = ref_window
            best_subread = subread
    
    if best_identity >= min_identity:
        return best_ref, best_subread, best_start, best_start + len(subread), best_identity
    else:
        return "", "", 0, 0, 0.0

def find_aa_mutations(ref_aa: str, query_aa: str, start_pos: int) -> list:
    """
    Find amino acid mutations between reference and query sequences.
    Returns list of (position, reference_aa, mutated_aa) tuples.
    """
    mutations = []
    min_len = min(len(ref_aa), len(query_aa))
    
    for i in range(min_len):
        if ref_aa[i] != query_aa[i]:
            # Convert to 1-based indexing
            mutations.append((start_pos + i + 1, ref_aa[i], query_aa[i]))
    
    return mutations

def main():
    # Load reference
    reference_record = next(SeqIO.parse("sp24_dna.fa", "fasta"))
    reference_seq = str(reference_record.seq).upper()
    print(f"Reference length: {len(reference_seq)}")
    
    # Load a few reads
    reads = list(SeqIO.parse("L1e8e291d_LoopSeqSample_admin/L1e8e291d_output/L1e8e291d_sample_B3/L1e8e291d_sample_B3_contig_list_trimmed.fa", "fasta"))
    print(f"Total reads: {len(reads)}")
    
    # Parameters
    gene_start = "ATGCCCAAAATAAATACATTTAATT"
    gene_length = 1233
    min_length = 1500
    
    # Test first 5 reads
    for i, read in enumerate(reads[:5]):
        print(f"\n--- Read {i+1} ---")
        print(f"Read length: {len(read.seq)}")
        
        if len(read.seq) < min_length:
            print("Read too short, skipping")
            continue
            
        query_seq = str(read.seq).upper()
        
        # Extract subread
        subread, mismatches, is_reverse = extract_subread_with_alignment(query_seq, gene_start, gene_length, max_mismatches=2)
        
        if not subread:
            print("No gene start found")
            continue
            
        print(f"Subread length: {len(subread)}")
        print(f"Mismatches: {mismatches}")
        print(f"Is reverse: {is_reverse}")
        print(f"Subread start: {subread[:50]}...")
        print(f"Gene start: {gene_start}")
        
        if mismatches > 2:
            print("Too many mismatches, skipping")
            continue
            
        if len(subread) < gene_length:
            print("Subread too short, skipping")
            continue
            
        # Compare to reference
        ref_aligned, subread_aligned, start, end, identity = compare_subread_to_reference(subread, reference_seq, min_identity=0.0)
        
        print(f"Identity: {identity:.3f}")
        print(f"Ref aligned: {ref_aligned[:50]}...")
        print(f"Subread aligned: {subread_aligned[:50]}...")
        
        if not ref_aligned or not subread_aligned:
            print("No alignment found")
            continue
            
        if identity < 0.9:
            print("Identity too low, skipping")
            continue
            
        # Translate to amino acids
        ref_aa = translate_dna(ref_aligned)
        query_aa = translate_dna(subread_aligned)
        
        print(f"Reference AA length: {len(ref_aa)}")
        print(f"Query AA length: {len(query_aa)}")
        
        # Find mutations
        mutations = find_aa_mutations(ref_aa, query_aa, start)
        
        print(f"Mutations found: {len(mutations)}")
        for pos, ref_aa, mut_aa in mutations:
            print(f"  Position {pos}: {ref_aa} -> {mut_aa}")

if __name__ == "__main__":
    main()