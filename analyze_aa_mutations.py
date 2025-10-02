#!/usr/bin/env python3
"""
Amino acid mutation analysis script.

This script analyzes FASTA read files to identify amino acid mutations,
calculate mutation frequencies, and perform statistical analyses.
"""

import os
import sys
import argparse
from collections import Counter, defaultdict
from typing import List, Tuple, Dict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import warnings
warnings.filterwarnings('ignore')

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

def find_gene_start_with_alignment(query_seq: str, gene_start: str, max_mismatches: int = 2) -> Tuple[str, int, int, bool]:
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

def extract_subread_with_alignment(query_seq: str, gene_start: str, gene_length: int, max_mismatches: int = 2) -> Tuple[str, int, bool]:
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

def compare_subread_to_reference(subread: str, reference_seq: str, min_identity: float = 0.9) -> Tuple[str, str, int, int, float]:
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

def find_aa_mutations(ref_aa: str, query_aa: str, start_pos: int) -> List[Tuple[int, str, str]]:
    """
    Find amino acid mutations between reference and query sequences.
    Returns list of (position, reference_aa, mutated_aa) tuples.
    """
    mutations = []
    min_len = min(len(ref_aa), len(query_aa))
    
    for i in range(min_len):
        if ref_aa[i] != query_aa[i]:
            # Convert DNA position to amino acid position (1-based)
            # start_pos is DNA position, divide by 3 to get amino acid position
            aa_pos = (start_pos // 3) + i + 1
            mutations.append((aa_pos, ref_aa[i], query_aa[i]))
    
    return mutations

def run_frequency_filtering(all_mutations: List[List[Tuple[int, str, str]]], min_frequency: float = 0.05) -> List[int]:
    """Filter positions by mutation frequency and return list of positions that pass the filter."""
    print("Starting position frequency filtering...")
    
    # Count mutations per position
    position_counts = defaultdict(int)
    total_reads = len(all_mutations)
    
    for mutations in all_mutations:
        for pos, _, _ in mutations:
            position_counts[pos] += 1
    
    # Filter positions by frequency
    filtered_positions = []
    for pos in sorted(position_counts.keys()):
        frequency = position_counts[pos] / total_reads
        if frequency > min_frequency:
            filtered_positions.append(pos)
    
    n_positions = len(filtered_positions)
    print(f"Position filtering: {len(position_counts)} total positions, {n_positions} positions with mutations in >{min_frequency*100:.1f}% of reads")
    
    if n_positions < 2:
        print("Insufficient positions for pair analysis after filtering")
        print(f"Need at least 2 positions with >{min_frequency*100:.1f}% mutation frequency")
        return []
    
    print(f"Frequency filtering completed: {n_positions} positions passed filter")
    return filtered_positions

def format_mutation_label(mutation_label: str, reference_seq: str = None) -> str:
    """Format mutation label from Pos161_I-Pos224_S to S161I-F224S format."""
    if not reference_seq or 'Pos' not in mutation_label:
        return mutation_label
    
    parts = mutation_label.split('-')
    formatted_parts = []
    
    for part in parts:
        if 'Pos' in part and '_' in part:
            pos_part = part.replace('Pos', '')
            pos_num, mut_aa = pos_part.split('_')
            pos_idx = int(pos_num) - 1  # Convert to 0-based for indexing
            
            # Get reference amino acid at this position
            if pos_idx < len(reference_seq):
                ref_aa = reference_seq[pos_idx]
                formatted_parts.append(f"{ref_aa}{pos_num}{mut_aa}")
            else:
                formatted_parts.append(f"X{pos_num}{mut_aa}")
        else:
            formatted_parts.append(part)
    
    return '-'.join(formatted_parts)

def plot_observed_vs_expected(mutation_counts: Counter, output_dir: str, filtered_positions: List[int], min_frequency: float = 0.05, reference_seq: str = None):
    """Plot observed vs expected frequency for mutation pairs using only DCA-filtered positions."""
    print("Starting observed vs expected frequency calculation...")
    
    if len(mutation_counts) < 2:
        print("Insufficient mutation data for pair analysis")
        return
    
    observed_freqs = []
    expected_freqs = []
    pair_labels = []  # For individual mutation pairs (e.g., Pos161_I-Pos224_S)
    position_pairs = []  # For position pairs (e.g., Pos161-Pos224)
    
    total_reads_with_mutations = sum(mutation_counts.values())
    if total_reads_with_mutations == 0:
        print("No reads with mutations to analyze for pairs.")
        return
    
    # Calculate individual mutation frequencies
    individual_mutation_counts = defaultdict(int)
    for mutation_set, count in mutation_counts.items():
        for pos, ref_aa, mut_aa in mutation_set:
            individual_mutation_counts[(pos, ref_aa, mut_aa)] += count
    
    # Filter individual mutations by position frequency
    filtered_individual_mutations = {
        (pos, ref_aa, mut_aa): count
        for (pos, ref_aa, mut_aa), count in individual_mutation_counts.items()
        if pos in filtered_positions
    }
    
    if len(filtered_positions) < 2:
        print("Insufficient positions passed frequency filter for pair analysis (need at least 2).")
        return
    
    print(f"Analyzing mutation pairs for {len(filtered_positions)} DCA-filtered positions")
    
    # Iterate over all unique pairs of filtered positions
    for i in range(len(filtered_positions)):
        for j in range(i + 1, len(filtered_positions)):
            pos1 = filtered_positions[i]
            pos2 = filtered_positions[j]
            
            # Get all mutations at pos1 and pos2
            mutations_at_pos1 = [(ref, mut) for (p, ref, mut), count in filtered_individual_mutations.items() if p == pos1]
            mutations_at_pos2 = [(ref, mut) for (p, ref, mut), count in filtered_individual_mutations.items() if p == pos2]
            
            # Consider all combinations of mutations at these two positions
            for ref1, mut1 in mutations_at_pos1:
                for ref2, mut2 in mutations_at_pos2:
                    # Count reads with both mutations
                    observed_count = 0
                    for mutation_set, count in mutation_counts.items():
                        has_mut1 = any(p == pos1 and m_aa == mut1 for p, _, m_aa in mutation_set)
                        has_mut2 = any(p == pos2 and m_aa == mut2 for p, _, m_aa in mutation_set)
                        if has_mut1 and has_mut2:
                            observed_count += count
                    
                    observed_freq = observed_count / total_reads_with_mutations
                    
                    # Expected frequency assuming independence
                    freq_mut1 = individual_mutation_counts.get((pos1, ref1, mut1), 0) / total_reads_with_mutations
                    freq_mut2 = individual_mutation_counts.get((pos2, ref2, mut2), 0) / total_reads_with_mutations
                    expected_freq = freq_mut1 * freq_mut2
                    
                    observed_freqs.append(observed_freq)
                    expected_freqs.append(expected_freq)
                    pair_labels.append(f"Pos{pos1}_{mut1}-Pos{pos2}_{mut2}")
                    position_pairs.append(f"Pos{pos1}-Pos{pos2}")
    
    if not observed_freqs:
        print("No mutation pairs found for DCA-filtered positions")
        return
    
    # Save raw data to CSV file (only for positions that passed the frequency filter)
    raw_data_file = os.path.join(output_dir, "mutation_pairs_raw_data.csv")
    raw_data = []
    
    # Only include data for positions that passed the frequency filter
    for i, (obs, exp, pair) in enumerate(zip(observed_freqs, expected_freqs, position_pairs)):
        # Extract position numbers from pair string
        if '-' in pair and pair.count('-') == 1:
            pos1, pos2 = pair.replace('Pos', '').split('-')
            pos1_int, pos2_int = int(pos1), int(pos2)
            
            # Only include if both positions are in the filtered positions
            if pos1_int in filtered_positions and pos2_int in filtered_positions:
                raw_data.append({
                    'position_1': pos1_int,
                    'position_2': pos2_int,
                    'mutation_pair': pair_labels[i] if i < len(pair_labels) else f"Pos{pos1}-Pos{pos2}",
                    'observed_frequency': obs,
                    'expected_frequency': exp,
                    'position_pair': pair
                })
    
    raw_df = pd.DataFrame(raw_data)
    raw_df.to_csv(raw_data_file, index=False)
    print(f"Saved raw mutation pairs data to {raw_data_file} ({len(raw_data)} pairs from {len(filtered_positions)} filtered positions)")
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Create color map for position pairs
    unique_pairs = list(set(position_pairs))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_pairs)))
    pair_to_color = {pair: colors[i] for i, pair in enumerate(unique_pairs)}
    
    # Plot points with legend
    for pair in unique_pairs:
        pair_indices = [i for i, p in enumerate(position_pairs) if p == pair]
        if pair_indices:
            pair_obs = [observed_freqs[i] for i in pair_indices]
            pair_exp = [expected_freqs[i] for i in pair_indices]
            plt.scatter(pair_exp, pair_obs, c=[pair_to_color[pair]], alpha=0.7, s=50, label=pair)
            
            # Add labels for significant deviations with better spacing
            label_count = 0
            for i, idx in enumerate(pair_indices):
                obs = observed_freqs[idx]
                exp = expected_freqs[idx]
                if exp > 0.01:  # Only label points with expected frequency > 0.01
                    ratio = obs / exp
                    if ratio > 5 or ratio < 0.2:
                        mutation_label = pair_labels[idx] if idx < len(pair_labels) else f"Pos{pair.replace('Pos', '').replace('-', '-Pos')}"
                        formatted_label = format_mutation_label(mutation_label, reference_seq)
                        
                        # Better label positioning to avoid overlap
                        distance = 20 + (label_count % 3) * 15  # Vary distance
                        offset_x = distance * (1 if label_count % 2 == 0 else -1)  # Alternate left/right
                        offset_y = distance * (1 if label_count % 4 < 2 else -1)  # Alternate up/down
                        
                        plt.annotate(formatted_label, (exp, obs), 
                                   xytext=(offset_x, offset_y), 
                                   textcoords='offset points', 
                                   fontsize=6, alpha=0.8,
                                   bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
                                   rotation=0)  # All labels horizontal
                        label_count += 1
    
    # Add y=x line
    max_freq = max(max(observed_freqs), max(expected_freqs))
    plt.plot([0, max_freq], [0, max_freq], 'k--', alpha=0.5, label='y=x')
    
    plt.xlabel('Expected Frequency (Independence)')
    plt.ylabel('Observed Frequency')
    plt.title(f'Observed vs Expected Mutation Pair Frequencies\n(DCA-filtered positions, >{min_frequency*100:.1f}% mutation frequency)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plot_file = os.path.join(output_dir, "observed_vs_expected_frequencies.pdf")
    try:
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Saved observed vs expected plot to {plot_file}")
    except Exception as e:
        print(f"Warning: Could not save PDF plot: {e}")
        # Try saving as PNG instead
        png_file = os.path.join(output_dir, "observed_vs_expected_frequencies.png")
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        print(f"Saved observed vs expected plot as PNG to {png_file}")
    finally:
        plt.close()
    print("Observed vs expected frequency calculation completed")

def plot_observed_vs_expected_triplets(mutation_counts: Counter, output_dir: str, filtered_positions: List[int], min_frequency: float = 0.05, reference_seq: str = None):
    """Plot observed vs expected frequency for mutation triplets using only DCA-filtered positions."""
    print("Starting observed vs expected frequency calculation for triplets...")
    
    observed_freqs = []
    expected_freqs = []
    triplet_labels = []  # For individual mutation triplets (e.g., Pos1_A-Pos2_G-Pos3_L)
    position_triplets = []  # For position triplets (e.g., Pos1-Pos2-Pos3)
    
    total_reads_with_mutations = sum(mutation_counts.values())
    if total_reads_with_mutations == 0:
        print("No reads with mutations to analyze for triplets.")
        return
    
    # Calculate individual mutation frequencies
    individual_mutation_counts = defaultdict(int)
    for mutation_set, count in mutation_counts.items():
        for pos, ref_aa, mut_aa in mutation_set:
            individual_mutation_counts[(pos, ref_aa, mut_aa)] += count
    
    # Filter individual mutations by position frequency
    filtered_individual_mutations = {
        (pos, ref_aa, mut_aa): count
        for (pos, ref_aa, mut_aa), count in individual_mutation_counts.items()
        if pos in filtered_positions
    }
    
    if len(filtered_positions) < 3:
        print("Insufficient positions passed frequency filter for triplet analysis (need at least 3).")
        return
    
    print(f"Analyzing mutation triplets for {len(filtered_positions)} DCA-filtered positions")
    
    # Iterate over all unique triplets of filtered positions
    for i in range(len(filtered_positions)):
        for j in range(i + 1, len(filtered_positions)):
            for k in range(j + 1, len(filtered_positions)):
                pos1 = filtered_positions[i]
                pos2 = filtered_positions[j]
                pos3 = filtered_positions[k]
                
                # Get all mutations at pos1, pos2, and pos3
                mutations_at_pos1 = [(ref, mut) for (p, ref, mut), count in filtered_individual_mutations.items() if p == pos1]
                mutations_at_pos2 = [(ref, mut) for (p, ref, mut), count in filtered_individual_mutations.items() if p == pos2]
                mutations_at_pos3 = [(ref, mut) for (p, ref, mut), count in filtered_individual_mutations.items() if p == pos3]
                
                # Consider all combinations of mutations at these three positions
                for ref1, mut1 in mutations_at_pos1:
                    for ref2, mut2 in mutations_at_pos2:
                        for ref3, mut3 in mutations_at_pos3:
                            # Count reads with all three mutations
                            observed_count = 0
                            for mutation_set, count in mutation_counts.items():
                                has_mut1 = any(p == pos1 and m_aa == mut1 for p, _, m_aa in mutation_set)
                                has_mut2 = any(p == pos2 and m_aa == mut2 for p, _, m_aa in mutation_set)
                                has_mut3 = any(p == pos3 and m_aa == mut3 for p, _, m_aa in mutation_set)
                                if has_mut1 and has_mut2 and has_mut3:
                                    observed_count += count
                            
                            observed_freq = observed_count / total_reads_with_mutations
                            
                            # Expected frequency assuming independence
                            freq_mut1 = individual_mutation_counts.get((pos1, ref1, mut1), 0) / total_reads_with_mutations
                            freq_mut2 = individual_mutation_counts.get((pos2, ref2, mut2), 0) / total_reads_with_mutations
                            freq_mut3 = individual_mutation_counts.get((pos3, ref3, mut3), 0) / total_reads_with_mutations
                            expected_freq = freq_mut1 * freq_mut2 * freq_mut3
                            
                            observed_freqs.append(observed_freq)
                            expected_freqs.append(expected_freq)
                            triplet_labels.append(f"Pos{pos1}_{mut1}-Pos{pos2}_{mut2}-Pos{pos3}_{mut3}")
                            position_triplets.append(f"Pos{pos1}-Pos{pos2}-Pos{pos3}")
    
    if not observed_freqs:
        print("No mutation triplets found for DCA-filtered positions")
        return
    
    # Save raw data to CSV file (only for positions that passed the frequency filter)
    raw_data_file = os.path.join(output_dir, "mutation_triplets_raw_data.csv")
    raw_data = []
    
    # Only include data for positions that passed the frequency filter
    for i, (obs, exp, triplet) in enumerate(zip(observed_freqs, expected_freqs, position_triplets)):
        # Extract position numbers from triplet string
        if '-' in triplet and triplet.count('-') == 2:
            pos1, pos2, pos3 = triplet.replace('Pos', '').split('-')
            pos1_int, pos2_int, pos3_int = int(pos1), int(pos2), int(pos3)
            
            # Only include if all three positions are in the filtered positions
            if pos1_int in filtered_positions and pos2_int in filtered_positions and pos3_int in filtered_positions:
                raw_data.append({
                    'position_1': pos1_int,
                    'position_2': pos2_int,
                    'position_3': pos3_int,
                    'mutation_triplet': triplet_labels[i] if i < len(triplet_labels) else f"Pos{pos1}-Pos{pos2}-Pos{pos3}",
                    'observed_frequency': obs,
                    'expected_frequency': exp,
                    'position_triplet': triplet
                })
    
    raw_df = pd.DataFrame(raw_data)
    raw_df.to_csv(raw_data_file, index=False)
    print(f"Saved raw mutation triplets data to {raw_data_file} ({len(raw_data)} triplets from {len(filtered_positions)} filtered positions)")
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Create color map for position triplets
    unique_triplets = list(set(position_triplets))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_triplets)))
    triplet_to_color = {triplet: colors[i] for i, triplet in enumerate(unique_triplets)}
    
    # Plot points with legend
    for triplet in unique_triplets:
        triplet_indices = [i for i, t in enumerate(position_triplets) if t == triplet]
        if triplet_indices:
            triplet_obs = [observed_freqs[i] for i in triplet_indices]
            triplet_exp = [expected_freqs[i] for i in triplet_indices]
            plt.scatter(triplet_exp, triplet_obs, c=[triplet_to_color[triplet]], alpha=0.7, s=50, label=triplet)
            
            # Add labels for significant deviations with better spacing
            label_count = 0
            for i, idx in enumerate(triplet_indices):
                obs = observed_freqs[idx]
                exp = expected_freqs[idx]
                if exp > 0.01:  # Only label points with expected frequency > 0.01
                    ratio = obs / exp
                    if ratio > 5 or ratio < 0.2:
                        mutation_label = triplet_labels[idx] if idx < len(triplet_labels) else f"Pos{triplet.replace('Pos', '').replace('-', '-Pos')}"
                        formatted_label = format_mutation_label(mutation_label, reference_seq)
                        
                        # Better label positioning to avoid overlap
                        distance = 25 + (label_count % 4) * 10  # Vary distance
                        offset_x = distance * (1 if label_count % 2 == 0 else -1)  # Alternate left/right
                        offset_y = distance * (1 if label_count % 4 < 2 else -1)  # Alternate up/down
                        
                        plt.annotate(formatted_label, (exp, obs), 
                                   xytext=(offset_x, offset_y), 
                                   textcoords='offset points', 
                                   fontsize=5, alpha=0.8,
                                   bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
                                   rotation=0)  # All labels horizontal
                        label_count += 1
    
    # Add y=x line
    max_freq = max(max(observed_freqs), max(expected_freqs))
    plt.plot([0, max_freq], [0, max_freq], 'k--', alpha=0.5, label='y=x')
    
    plt.xlabel('Expected Frequency (Independence)')
    plt.ylabel('Observed Frequency')
    plt.title(f'Observed vs Expected Mutation Triplet Frequencies\n(DCA-filtered positions, >{min_frequency*100:.1f}% mutation frequency)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plot_file = os.path.join(output_dir, "observed_vs_expected_triplets.pdf")
    try:
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Saved observed vs expected triplets plot to {plot_file}")
    except Exception as e:
        print(f"Warning: Could not save PDF plot: {e}")
        # Try saving as PNG instead
        png_file = os.path.join(output_dir, "observed_vs_expected_triplets.png")
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        print(f"Saved observed vs expected triplets plot as PNG to {png_file}")
    finally:
        plt.close()
    print("Observed vs expected frequency calculation for triplets completed")

def calculate_entropy(mutation_counts: Counter) -> float:
    """Calculate normalized entropy of genotype distribution."""
    total_counts = sum(mutation_counts.values())
    if total_counts == 0:
        return 0.0
    
    # Calculate frequencies
    frequencies = [count / total_counts for count in mutation_counts.values()]
    
    # Calculate entropy
    entropy = -sum(p * np.log2(p) for p in frequencies if p > 0)
    
    # Normalize by maximum possible entropy (log2 of number of unique genotypes)
    max_entropy = np.log2(len(mutation_counts))
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0.0
    
    return normalized_entropy

def analyze_mutations(input_fasta: str, reference_fasta: str, min_length: int, gene_start: str, gene_length: int, output_dir: str, frequency_cutoff: float = 0.05):
    """Main analysis function."""
    print(f"Loading reference from {reference_fasta}")
    reference_record = next(SeqIO.parse(reference_fasta, "fasta"))
    reference_seq = str(reference_record.seq).upper()
    print(f"Reference length: {len(reference_seq)}")
    
    # Load the other reference for cross-reference checking
    other_reference_fasta = "sp24_dna.fa" if "sp23" in reference_fasta else "sp23_dna.fa"
    print(f"Loading other reference from {other_reference_fasta} for cross-reference checking")
    other_reference_record = next(SeqIO.parse(other_reference_fasta, "fasta"))
    other_reference_seq = str(other_reference_record.seq).upper()
    print(f"Other reference length: {len(other_reference_seq)}")
    
    print(f"Loading reads from {input_fasta}")
    reads = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Total reads: {len(reads)}")
    
    # Filter by length
    filtered_reads = [read for read in reads if len(read.seq) >= min_length]
    print(f"Reads after length filter (>{min_length}): {len(filtered_reads)}")
    
    if not filtered_reads:
        print("No reads passed length filter. Exiting.")
        return
    
    # Analyze mutations
    all_mutations = []
    mutation_counts = Counter()
    processed = 0
    skipped_no_gene_start = 0
    skipped_high_mismatch = 0
    skipped_short_subread = 0
    skipped_low_identity = 0
    skipped_cross_reference = 0
    no_mutations_count = 0
    
    print("Analyzing mutations...")
    for read in filtered_reads:
        processed += 1
        if processed % 100 == 0:
            print(f"Processed {processed} reads...")
        
        query_seq = str(read.seq).upper()
        
        # Extract subread using gene_start with alignment and mismatch threshold
        subread, mismatches, is_reverse = extract_subread_with_alignment(query_seq, gene_start, gene_length, max_mismatches=2)
        
        if not subread:
            skipped_no_gene_start += 1
            continue
        
        if mismatches > 2:
            skipped_high_mismatch += 1
            continue
        
        if len(subread) < gene_length:
            skipped_short_subread += 1
            continue
        
        # Compare subread directly to reference (no alignment)
        ref_aligned, subread_aligned, start, end, identity = compare_subread_to_reference(subread, reference_seq, min_identity=0.9)
        
        if not ref_aligned or not subread_aligned:
            continue
        
        if identity < 0.9:
            skipped_low_identity += 1
            continue
        
        # Cross-reference check: compare to other reference
        other_ref_aligned, other_subread_aligned, other_start, other_end, other_identity = compare_subread_to_reference(subread, other_reference_seq, min_identity=0.0)
        
        # If subread matches >90% to the other reference AND better than target reference, discard it
        # Temporarily disabled for debugging
        # if other_identity > 0.9 and other_identity > identity:
        #     skipped_cross_reference += 1
        #     continue
        
        # Translate DNA sequences to amino acid sequences
        ref_aa = translate_dna(ref_aligned)
        query_aa = translate_dna(subread_aligned)
        
        # Find amino acid mutations
        mutations = find_aa_mutations(ref_aa, query_aa, start)
        
        if mutations:
            all_mutations.append(mutations)
            # Create mutation set string - mutations are already in amino acid format
            mutation_set = tuple(sorted(mutations))
            mutation_counts[mutation_set] += 1
        else:
            no_mutations_count += 1
    
    print(f"Reads with mutations: {len(all_mutations)}")
    
    # Calculate entropy
    normalized_entropy = calculate_entropy(mutation_counts)
    print(f"Normalized entropy of genotype distribution: {normalized_entropy:.3f}")
    
    # Run frequency filtering for downstream analyses
    coupling_filtered_positions = run_frequency_filtering(all_mutations, min_frequency=frequency_cutoff)
    
    # Translate reference DNA to amino acids for proper labeling
    reference_aa_seq = translate_dna(reference_seq)
    
    # Create observed vs expected frequency plot using coupling-filtered positions
    plot_observed_vs_expected(mutation_counts, output_dir, coupling_filtered_positions, min_frequency=frequency_cutoff, reference_seq=reference_aa_seq)
    
    # Create observed vs expected frequency plot for mutation triplets
    plot_observed_vs_expected_triplets(mutation_counts, output_dir, coupling_filtered_positions, min_frequency=frequency_cutoff, reference_seq=reference_aa_seq)
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    
    # Save mutation table
    mutation_table_data = []
    mutation_table_data.append({"Mutation_Set": "No_mutations", "Count": no_mutations_count})
    
    for mutation_set, count in mutation_counts.most_common():
        mutation_str = "; ".join([f"{ref_aa}{pos}{mut_aa}" for pos, ref_aa, mut_aa in mutation_set])
        mutation_table_data.append({"Mutation_Set": mutation_str, "Count": count})
    
    mutation_df = pd.DataFrame(mutation_table_data)
    mutation_table_file = os.path.join(output_dir, "mutation_table.csv")
    mutation_df.to_csv(mutation_table_file, index=False)
    print(f"Saved mutation table to {mutation_table_file}")
    
    # Save summary statistics
    summary_file = os.path.join(output_dir, "summary_statistics.txt")
    with open(summary_file, 'w') as f:
        f.write(f"reads_total\t{len(reads)}\n")
        f.write(f"reads_length_pass\t{len(filtered_reads)}\n")
        f.write(f"reads_length_fail\t{len(reads) - len(filtered_reads)}\n")
        f.write(f"reads_gene_start_found\t{len(filtered_reads) - skipped_no_gene_start}\n")
        f.write(f"reads_discard_high_mismatch\t{skipped_high_mismatch}\n")
        f.write(f"reads_discard_short_subread\t{skipped_short_subread}\n")
        f.write(f"reads_discard_low_identity\t{skipped_low_identity}\n")
        f.write(f"reads_discard_cross_reference\t{skipped_cross_reference}\n")
        f.write(f"reads_no_mutations\t{no_mutations_count}\n")
        f.write(f"reads_with_mutations\t{len(all_mutations)}\n")
        f.write(f"unique_mutation_sets\t{len(mutation_counts)}\n")
        f.write(f"coupling_filtered_positions\t{len(coupling_filtered_positions)}\n")
        f.write(f"position_frequency_threshold\t{frequency_cutoff*100:.1f}\n")
        f.write(f"normalized_entropy\t{normalized_entropy:.3f}\n")
        f.write(f"coupling_analysis_skipped\tTrue\n")
    
    print(f"Saved summary statistics to {summary_file}")
    
    # Print summary
    print(f"\nSummary:")
    print(f"Total reads: {len(reads)}")
    print(f"Reads passing length filter: {len(filtered_reads)}")
    print(f"Reads failing length filter: {len(reads) - len(filtered_reads)}")
    print(f"Reads with gene_start found: {len(filtered_reads) - skipped_no_gene_start}")
    print(f"Reads discarded (high mismatch >2): {skipped_high_mismatch}")
    print(f"Reads discarded (short subread): {skipped_short_subread}")
    print(f"Reads discarded (low identity <90%): {skipped_low_identity}")
    print(f"Reads discarded (better match to other reference): {skipped_cross_reference}")
    print(f"Reads with no mutations: {no_mutations_count}")
    print(f"Reads with mutations: {len(all_mutations)}")
    print(f"Unique mutation sets: {len(mutation_counts)}")
    print(f"Coupling-filtered positions: {len(coupling_filtered_positions)} with >{frequency_cutoff*100:.1f}% mutation frequency")
    print(f"Normalized entropy: {normalized_entropy:.3f}")
    print(f"Coupling analysis: Skipped (EVcouplings requires pre-computed model files)")
    if mutation_counts:
        most_common = mutation_counts.most_common(1)[0]
        mutation_str = "; ".join([f"{ref_aa}{pos}{mut_aa}" for pos, ref_aa, mut_aa in most_common[0]])
        print(f"Most common mutation set: {mutation_str} (count: {most_common[1]})")
    
    print(f"\nAnalysis complete. Results saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Analyze amino acid mutations in FASTA reads')
    parser.add_argument('--input_fasta', required=True, help='Input FASTA file with reads')
    parser.add_argument('--reference_fasta', required=True, help='Reference FASTA file')
    parser.add_argument('--min_length', type=int, default=1500, help='Minimum read length (default: 1500)')
    parser.add_argument('--gene_start', required=True, help='Gene start sequence to align')
    parser.add_argument('--gene_length', type=int, required=True, help='Length of gene to analyze')
    parser.add_argument('--frequency_cutoff', type=float, default=0.05, help='Frequency cutoff for position filtering (default: 0.05)')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    print("Dependencies imported: pandas, numpy, biopython")
    print("Analysis initiated")
    print(f"Input FASTA: {args.input_fasta}")
    print(f"Reference FASTA: {args.reference_fasta}")
    print(f"Min length: {args.min_length}")
    print(f"Gene start: {args.gene_start}")
    print(f"Gene length: {args.gene_length}")
    print(f"Frequency cutoff: {args.frequency_cutoff*100:.1f}%")
    print(f"Output directory: {args.output_dir}")
    print()
    
    analyze_mutations(
        args.input_fasta,
        args.reference_fasta,
        args.min_length,
        args.gene_start,
        args.gene_length,
        args.output_dir,
        args.frequency_cutoff
    )

if __name__ == "__main__":
    main()
