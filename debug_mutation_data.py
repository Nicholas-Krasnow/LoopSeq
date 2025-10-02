#!/usr/bin/env python3
"""
Debug script to check mutation data format.
"""

from Bio import SeqIO
from Bio.Seq import Seq

def translate_dna(sequence: str) -> str:
    """Translate DNA sequence to amino acid sequence using standard genetic code."""
    clean_seq = ''.join(c for c in sequence.upper() if c in 'ATCG')
    if len(clean_seq) % 3 != 0:
        clean_seq = clean_seq[:-(len(clean_seq) % 3)]
    if len(clean_seq) == 0:
        return ""
    seq_obj = Seq(clean_seq)
    protein = seq_obj.translate()
    return str(protein)

def find_gene_start_with_alignment(query_seq: str, gene_start: str, max_mismatches: int = 2) -> tuple:
    query_upper = query_seq.upper()
    gene_start_upper = gene_start.upper()
    
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
    aligned_query, start_pos, end_pos, is_reverse = find_gene_start_with_alignment(query_seq, gene_start, max_mismatches)
    
    if start_pos == -1:
        return "", 0, False
    
    mismatches = sum(1 for a, b in zip(aligned_query, gene_start.upper()) if a != b)
    subread_start = start_pos
    subread_end = subread_start + gene_length
    
    if is_reverse:
        query_rc = str(Seq(query_seq).reverse_complement())
        if subread_start >= 0 and subread_end <= len(query_rc):
            subread = query_rc[subread_start:subread_end]
        else:
            return "", mismatches, is_reverse
    else:
        if subread_start >= 0 and subread_end <= len(query_seq):
            subread = query_seq[subread_start:subread_end]
        else:
            return "", mismatches, is_reverse
    
    return subread, mismatches, is_reverse

def compare_subread_to_reference(subread: str, reference_seq: str, min_identity: float = 0.9) -> tuple:
    if not subread or len(subread) == 0:
        return "", "", 0, 0, 0.0
    
    best_identity = 0.0
    best_start = 0
    best_ref = ""
    best_subread = ""
    
    for start in range(len(reference_seq) - len(subread) + 1):
        ref_window = reference_seq[start:start + len(subread)]
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
    mutations = []
    min_len = min(len(ref_aa), len(query_aa))
    
    for i in range(min_len):
        if ref_aa[i] != query_aa[i]:
            # Convert DNA position to amino acid position (1-based)
            aa_pos = (start_pos // 3) + i + 1
            mutations.append((aa_pos, ref_aa[i], query_aa[i]))
    
    return mutations

def main():
    # Load reference
    reference_record = next(SeqIO.parse("sp23_dna.fa", "fasta"))
    reference_seq = str(reference_record.seq).upper()
    print(f"Reference length: {len(reference_seq)}")
    
    # Load one read
    reads = list(SeqIO.parse("L1e8e291d_sample_A1_downsampled_10x.fa", "fasta"))
    read = reads[0]
    
    gene_start = "ATGCCCAAAATAAATACATTTAATT"
    gene_length = 1233
    
    query_seq = str(read.seq).upper()
    
    # Extract subread
    subread, mismatches, is_reverse = extract_subread_with_alignment(query_seq, gene_start, gene_length, max_mismatches=2)
    
    if not subread:
        print("No gene start found")
        return
    
    # Compare to reference
    ref_aligned, subread_aligned, start, end, identity = compare_subread_to_reference(subread, reference_seq, min_identity=0.9)
    
    if not ref_aligned or not subread_aligned:
        print("No alignment found")
        return
    
    print(f"DNA alignment start: {start}")
    print(f"DNA alignment end: {end}")
    
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
    
    # Test mutation formatting
    if mutations:
        mutation_set = tuple(sorted(mutations))
        print(f"\nMutation set: {mutation_set}")
        
        # Test the formatting used in the script
        mutation_str = "; ".join([f"{ref_aa}{pos}{mut_aa}" for pos, ref_aa, mut_aa in mutation_set])
        print(f"Formatted mutation string: {mutation_str}")

if __name__ == "__main__":
    main()
