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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import warnings
warnings.filterwarnings('ignore')

def translate_dna(sequence):
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

def find_gene_start_with_alignment(query_seq, gene_start, max_mismatches=2):
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

def extract_subread_with_alignment(query_seq, gene_start, gene_length, max_mismatches=2):
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

def compare_subread_to_reference(subread, reference_seq, min_identity=0.9):
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

def find_aa_mutations(ref_aa, query_aa, start_pos):
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

def run_frequency_filtering(all_mutations, min_frequency=0.05):
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
    print("Position filtering: {0} total positions, {1} positions with mutations in >{2:.1f}% of reads".format(len(position_counts), n_positions, min_frequency*100))
    
    if n_positions < 2:
        print("Insufficient positions for pair analysis after filtering")
        print("Need at least 2 positions with >{0:.1f}% mutation frequency".format(min_frequency*100))
        return []
    
    print("Frequency filtering completed: {0} positions passed filter".format(n_positions))
    return filtered_positions

def format_mutation_label(mutation_label, reference_seq=None):
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
                formatted_parts.append("{0}{1}{2}".format(ref_aa, pos_num, mut_aa))
            else:
                formatted_parts.append("X{0}{1}".format(pos_num, mut_aa))
        else:
            formatted_parts.append(part)
    
    return '-'.join(formatted_parts)

def plot_observed_vs_expected(mutation_counts, output_dir, filtered_positions, min_frequency=0.05, reference_seq=None):
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
    
    print("Analyzing mutation pairs for {0} DCA-filtered positions".format(len(filtered_positions)))
    
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
                    pair_labels.append("Pos{0}_{1}-Pos{2}_{3}".format(pos1, mut1, pos2, mut2))
                    position_pairs.append("Pos{0}-Pos{1}".format(pos1, pos2))
    
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
                    'mutation_pair': pair_labels[i] if i < len(pair_labels) else "Pos{0}-Pos{1}".format(pos1, pos2),
                    'observed_frequency': obs,
                    'expected_frequency': exp,
                    'position_pair': pair
                })
    
    raw_df = pd.DataFrame(raw_data)
    raw_df.to_csv(raw_data_file, index=False)
    print("Saved raw mutation pairs data to {0} ({1} pairs from {2} filtered positions)".format(raw_data_file, len(raw_data), len(filtered_positions)))
    
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
                # New labeling criteria: (ratio > 5 or < 0.2) AND (obs > 0.1 or exp > 0.1)
                ratio = obs / exp if exp > 0 else float('inf')
                if (ratio > 5 or ratio < 0.2) and (obs > 0.1 or exp > 0.1):
                        mutation_label = pair_labels[idx] if idx < len(pair_labels) else "Pos{0}".format(pair.replace('Pos', '').replace('-', '-Pos'))
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
                                   rotation=0,  # All labels horizontal
                                   arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5, lw=0.5))  # Add line to label
                        label_count += 1
    
    # Add y=x line
    max_freq = max(max(observed_freqs), max(expected_freqs))
    plt.plot([0, max_freq], [0, max_freq], 'k--', alpha=0.5, label='y=x')
    
    plt.xlabel('Expected Frequency (Independence)')
    plt.ylabel('Observed Frequency')
    plt.title("Observed vs Expected Mutation Pair Frequencies\n(DCA-filtered positions, >{0:.1f}% mutation frequency)".format(min_frequency*100))
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plot_file = os.path.join(output_dir, "observed_vs_expected_frequencies.pdf")
    try:
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print("Saved observed vs expected plot to {0}".format(plot_file))
    except Exception as e:
        print("Warning: Could not save PDF plot: {0}".format(e))
        # Try saving as PNG instead
        png_file = os.path.join(output_dir, "observed_vs_expected_frequencies.png")
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        print("Saved observed vs expected plot as PNG to {0}".format(png_file))
    finally:
        plt.close()
    print("Observed vs expected frequency calculation completed")


def calculate_entropy(mutation_counts):
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

def plot_mutation_set_distribution(mutation_counts, output_dir, other_sample_dirs=None, min_count=100):
    """Plot distribution of mutation set counts with >100 reads, colored by uniqueness."""
    print("Creating mutation set distribution plot...")
    
    # Filter mutation sets with >100 reads
    high_count_mutations = {mutation_set: count for mutation_set, count in mutation_counts.items() if count > min_count}
    
    if not high_count_mutations:
        print("No mutation sets found with >{0} reads".format(min_count))
        return
    
    # Get mutation sets from other samples for comparison
    other_sample_mutations = set()
    if other_sample_dirs:
        print("Checking for mutation sets in {0} other sample directories...".format(len(other_sample_dirs)))
        for sample_dir in other_sample_dirs:
            mutation_table_file = os.path.join(sample_dir, "mutation_table.csv")
            if os.path.exists(mutation_table_file):
                try:
                    other_df = pd.read_csv(mutation_table_file)
                    # Extract mutation sets from other samples (skip "No_mutations" row)
                    other_mutations = set()
                    for _, row in other_df.iterrows():
                        if row['Mutation_Set'] != 'No_mutations' and row['Count'] > min_count:
                            # Parse mutation set string back to tuple format
                            mutation_str = row['Mutation_Set']
                            mutations = []
                            for mut in mutation_str.split('; '):
                                if len(mut) >= 3:  # At least ref_aa + pos + mut_aa
                                    # Extract position and amino acids
                                    pos = int(''.join(c for c in mut if c.isdigit()))
                                    ref_aa = mut[0]
                                    mut_aa = mut[-1]
                                    mutations.append((pos, ref_aa, mut_aa))
                            if mutations:
                                other_mutations.add(tuple(sorted(mutations)))
                    other_sample_mutations.update(other_mutations)
                    print("Found {0} high-count mutation sets in {1}".format(len(other_mutations), sample_dir))
                except Exception as e:
                    print("Warning: Could not read mutation table from {0}: {1}".format(sample_dir, e))
    
    # Determine colors for each mutation set
    colors = []
    unique_count = 0
    shared_count = 0
    
    for mutation_set in high_count_mutations.keys():
        if mutation_set in other_sample_mutations:
            colors.append('gray')
            shared_count += 1
        else:
            colors.append('green')
            unique_count += 1
    
    # Create plot
    plt.figure(figsize=(15, 8))
    
    # Sort mutation sets by count for better visualization
    sorted_mutations = sorted(high_count_mutations.items(), key=lambda x: x[1], reverse=True)
    mutation_sets, counts = zip(*sorted_mutations)
    
    # Create mutation set labels for x-axis
    mutation_labels = []
    for mutation_set in mutation_sets:
        mutation_str = "; ".join(["{0}{1}{2}".format(ref_aa, pos, mut_aa) for pos, ref_aa, mut_aa in mutation_set])
        mutation_labels.append(mutation_str)
    
    # Create bars
    bars = plt.bar(range(len(mutation_sets)), counts, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Customize plot
    plt.xlabel('Mutation Set', fontsize=12)
    plt.ylabel('Read Count', fontsize=12)
    plt.title("Distribution of Mutation Sets with >{0} Reads\n(Green=Unique to this sample, Gray=Found in other samples)".format(min_count), fontsize=14)
    plt.xticks(range(len(mutation_sets)), mutation_labels, rotation=45, ha='right', fontsize=8)
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', alpha=0.7, label="Unique to this sample ({0})".format(unique_count)),
        Patch(facecolor='gray', alpha=0.7, label="Found in other samples ({0})".format(shared_count))
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    # Add count annotations on bars
    for i, (bar, count) in enumerate(zip(bars, counts)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01, 
                str(count), ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = os.path.join(output_dir, "mutation_set_distribution.pdf")
    try:
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print("Saved mutation set distribution plot to {0}".format(plot_file))
    except Exception as e:
        print("Warning: Could not save PDF plot: {0}".format(e))
        # Try saving as PNG instead
        png_file = os.path.join(output_dir, "mutation_set_distribution.png")
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        print("Saved mutation set distribution plot as PNG to {0}".format(png_file))
    finally:
        plt.close()
    
    print("Mutation set distribution plot complete: {0} unique, {1} shared with other samples".format(unique_count, shared_count))

def analyze_mutations(input_fasta, reference_fasta, min_length, gene_start, gene_length, output_dir, frequency_cutoff=0.05, other_sample_dirs=None):
    """Main analysis function."""
    print("Loading reference from {0}".format(reference_fasta))
    reference_record = next(SeqIO.parse(reference_fasta, "fasta"))
    reference_seq = str(reference_record.seq).upper()
    print("Reference length: {0}".format(len(reference_seq)))
    
    # Load the other reference for cross-reference checking
    other_reference_fasta = "sp24_dna.fa" if "sp23" in reference_fasta else "sp23_dna.fa"
    print("Loading other reference from {0} for cross-reference checking".format(other_reference_fasta))
    other_reference_record = next(SeqIO.parse(other_reference_fasta, "fasta"))
    other_reference_seq = str(other_reference_record.seq).upper()
    print("Other reference length: {0}".format(len(other_reference_seq)))
    
    print("Loading reads from {0}".format(input_fasta))
    reads = list(SeqIO.parse(input_fasta, "fasta"))
    print("Total reads: {0}".format(len(reads)))
    
    # Filter by length
    filtered_reads = [read for read in reads if len(read.seq) >= min_length]
    print("Reads after length filter (>{0}): {1}".format(min_length, len(filtered_reads)))
    
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
    skipped_chimeras = 0
    no_mutations_count = 0
    
    print("Analyzing mutations...")
    for read in filtered_reads:
        processed += 1
        if processed % 100 == 0:
            print("Processed {0} reads...".format(processed))
        
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
        if other_identity > 0.9 and other_identity > identity:
            skipped_cross_reference += 1
            continue
        
        # Chimera filter: discard sequences with >93% match to both references
        if identity > 0.93 and other_identity > 0.93:
            skipped_chimeras += 1
            continue
        
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
    
    print("Reads with mutations: {0}".format(len(all_mutations)))
    
    # Calculate entropy
    normalized_entropy = calculate_entropy(mutation_counts)
    print("Normalized entropy of genotype distribution: {0:.3f}".format(normalized_entropy))
    
    # Run frequency filtering for downstream analyses
    coupling_filtered_positions = run_frequency_filtering(all_mutations, min_frequency=frequency_cutoff)
    
    # Translate reference DNA to amino acids for proper labeling
    reference_aa_seq = translate_dna(reference_seq)
    
    # Create observed vs expected frequency plot using coupling-filtered positions
    plot_observed_vs_expected(mutation_counts, output_dir, coupling_filtered_positions, min_frequency=frequency_cutoff, reference_seq=reference_aa_seq)
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    
    # Save mutation table
    mutation_table_data = []
    mutation_table_data.append({"Mutation_Set": "No_mutations", "Count": no_mutations_count})
    
    for mutation_set, count in mutation_counts.most_common():
        mutation_str = "; ".join(["{0}{1}{2}".format(ref_aa, pos, mut_aa) for pos, ref_aa, mut_aa in mutation_set])
        mutation_table_data.append({"Mutation_Set": mutation_str, "Count": count})
    
    mutation_df = pd.DataFrame(mutation_table_data)
    mutation_table_file = os.path.join(output_dir, "mutation_table.csv")
    mutation_df.to_csv(mutation_table_file, index=False)
    print("Saved mutation table to {0}".format(mutation_table_file))
    
    # Create mutation set distribution plot
    plot_mutation_set_distribution(mutation_counts, output_dir, other_sample_dirs, min_count=100)
    
    # Save summary statistics
    summary_file = os.path.join(output_dir, "summary_statistics.txt")
    with open(summary_file, 'w') as f:
        f.write("reads_total\t{0}\n".format(len(reads)))
        f.write("reads_length_pass\t{0}\n".format(len(filtered_reads)))
        f.write("reads_length_fail\t{0}\n".format(len(reads) - len(filtered_reads)))
        f.write("reads_gene_start_found\t{0}\n".format(len(filtered_reads) - skipped_no_gene_start))
        f.write("reads_discard_high_mismatch\t{0}\n".format(skipped_high_mismatch))
        f.write("reads_discard_short_subread\t{0}\n".format(skipped_short_subread))
        f.write("reads_discard_low_identity\t{0}\n".format(skipped_low_identity))
        f.write("reads_discard_cross_reference\t{0}\n".format(skipped_cross_reference))
        f.write("reads_discard_chimeras\t{0}\n".format(skipped_chimeras))
        f.write("reads_no_mutations\t{0}\n".format(no_mutations_count))
        f.write("reads_with_mutations\t{0}\n".format(len(all_mutations)))
        f.write("unique_mutation_sets\t{0}\n".format(len(mutation_counts)))
        f.write("coupling_filtered_positions\t{0}\n".format(len(coupling_filtered_positions)))
        f.write("position_frequency_threshold\t{0:.1f}\n".format(frequency_cutoff*100))
        f.write("normalized_entropy\t{0:.3f}\n".format(normalized_entropy))
        f.write("coupling_analysis_skipped\tTrue\n")
    
    print("Saved summary statistics to {0}".format(summary_file))
    
    # Print summary
    print("\nSummary:")
    print("Total reads: {0}".format(len(reads)))
    print("Reads passing length filter: {0}".format(len(filtered_reads)))
    print("Reads failing length filter: {0}".format(len(reads) - len(filtered_reads)))
    print("Reads with gene_start found: {0}".format(len(filtered_reads) - skipped_no_gene_start))
    print("Reads discarded (high mismatch >2): {0}".format(skipped_high_mismatch))
    print("Reads discarded (short subread): {0}".format(skipped_short_subread))
    print("Reads discarded (low identity <90%): {0}".format(skipped_low_identity))
    print("Reads discarded (better match to other reference): {0}".format(skipped_cross_reference))
    print("Reads discarded (chimeras >93% match to both references): {0}".format(skipped_chimeras))
    print("Reads with no mutations: {0}".format(no_mutations_count))
    print("Reads with mutations: {0}".format(len(all_mutations)))
    print("Unique mutation sets: {0}".format(len(mutation_counts)))
    print("Coupling-filtered positions: {0} with >{1:.1f}% mutation frequency".format(len(coupling_filtered_positions), frequency_cutoff*100))
    print("Normalized entropy: {0:.3f}".format(normalized_entropy))
    print("Coupling analysis: Skipped (EVcouplings requires pre-computed model files)")
    if mutation_counts:
        most_common = mutation_counts.most_common(1)[0]
        mutation_str = "; ".join(["{0}{1}{2}".format(ref_aa, pos, mut_aa) for pos, ref_aa, mut_aa in most_common[0]])
        print("Most common mutation set: {0} (count: {1})".format(mutation_str, most_common[1]))
    
    print("\nAnalysis complete. Results saved to {0}".format(output_dir))

def main():
    parser = argparse.ArgumentParser(description='Analyze amino acid mutations in FASTA reads')
    parser.add_argument('--input_fasta', required=True, help='Input FASTA file with reads')
    parser.add_argument('--reference_fasta', required=True, help='Reference FASTA file')
    parser.add_argument('--min_length', type=int, default=1500, help='Minimum read length (default: 1500)')
    parser.add_argument('--gene_start', required=True, help='Gene start sequence to align')
    parser.add_argument('--gene_length', type=int, required=True, help='Length of gene to analyze')
    parser.add_argument('--frequency_cutoff', type=float, default=0.05, help='Frequency cutoff for position filtering (default: 0.05)')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--other_sample_dirs', nargs='*', default=None, help='List of other sample directories to compare mutation sets against')
    
    args = parser.parse_args()
    
    print("Dependencies imported: pandas, numpy, biopython")
    print("Analysis initiated")
    print("Input FASTA: {0}".format(args.input_fasta))
    print("Reference FASTA: {0}".format(args.reference_fasta))
    print("Min length: {0}".format(args.min_length))
    print("Gene start: {0}".format(args.gene_start))
    print("Gene length: {0}".format(args.gene_length))
    print("Frequency cutoff: {0:.1f}%".format(args.frequency_cutoff*100))
    print("Output directory: {0}".format(args.output_dir))
    print()
    
    analyze_mutations(
        args.input_fasta,
        args.reference_fasta,
        args.min_length,
        args.gene_start,
        args.gene_length,
        args.output_dir,
        args.frequency_cutoff,
        args.other_sample_dirs
    )

if __name__ == "__main__":
    main()
