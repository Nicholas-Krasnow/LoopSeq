#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Amino acid mutation analysis script.

## How sequences are processed and filtered

Per input read, the script applies these steps (in order) inside `analyze_mutations()`:

- **Load reads**: parses all records from `--input_fasta` (FASTA).
- **Length filter**: discards reads shorter than `--min_length` (default 1500 nt).
- **Gene-start anchoring (with mismatches)**:
  - Searches for the provided `--gene_start` sequence within the read, allowing up to 2 mismatches.
  - If not found on the forward strand, repeats the search on the read’s reverse-complement.
  - If found, uses the best (lowest-mismatch) hit to define the extraction start.
- **Subread extraction**: extracts exactly `--gene_length` nucleotides starting at the gene-start hit.
  - Discards the read if extraction would run off the end (no full-length subread).
- **Reference identity filter (no gapped alignment)**:
  - Slides the extracted subread across the full reference and selects the window with maximum identity.
  - Requires identity ≥ 0.90 to be considered a match.
- **Cross-reference / chimera filters**:
  - Compares the same subread to the “other” reference (sp23 vs sp24).
  - Discards reads that match the other reference at >90% identity *and* better than the target reference.
  - Discards likely chimeras that match both references at >93% identity.
- **Translation and mutation calling**:
  - Translates the matched reference window and the matched subread to amino acids (after stripping non-ATCG and truncating to a multiple of 3).
  - Calls AA substitutions position-by-position; AA positions are reported 1-based and offset by the chosen DNA start.
- **Downstream position filtering**:
  - After collecting mutation sets across reads, positions are retained for coupling-style summaries only if mutated in more than `--frequency_cutoff` fraction of mutated reads.

Outputs include a per-genotype count table, summary statistics, and plots in `--output_dir`.
"""

import os
import sys
import re
import argparse
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import warnings
warnings.filterwarnings('ignore')

try:
    import matplotlib.font_manager as _fm
    # Prefer Arial; fall back to common metrically-compatible sans fonts on Linux.
    _preferred = ["Arial", "Liberation Sans", "DejaVu Sans"]
    _available = {f.name for f in _fm.fontManager.ttflist}
    for _name in _preferred:
        if _name in _available:
            plt.rcParams["font.family"] = _name
            break
    else:
        # Leave Matplotlib default if none are available.
        pass
except Exception:
    # If font manager isn't available for any reason, leave defaults.
    pass

AXIS_LABEL_FONTSIZE = 24
NON_AXIS_FONTSIZE = 8

NO_MUTATIONS_LABEL = "No_mutations"

# 1-48 sheet: only cells with this exact text (after strip) define extra sp23 samples in the Shared universe.
SHARED_UNIVERSE_MARKER_1_48 = "sp23-2 last"


def _cell_str(val):
    if val is None:
        return ""
    if isinstance(val, float) and val == int(val):
        return str(int(val))
    return str(val).strip()


def _normalize_header_key(name):
    return re.sub(r"\s+", "_", str(name).strip().lower())


def _sharing_universe_sample_names_from_sample_keys(base_dir, xlsx_path):
    """
    Sample names (output dir tags) allowed to define \"Shared\" mutation sets for each reference:

    - All rows on sheet \"initial_set\" (by reference from reference_fasta filename).
    - Sheet \"1-48\": only table cells whose text equals SHARED_UNIVERSE_MARKER_1_48 (sp23).

    Returns None if the workbook cannot be read (caller treats as unrestricted).
    Returns (set sp23 names, set sp24 names) on success (sets may be empty).
    """
    try:
        from openpyxl import load_workbook
    except ImportError:
        print("Warning: openpyxl not available; cannot load sample_keys for Shared universe")
        return None

    if not xlsx_path or not os.path.isfile(xlsx_path):
        return None

    allowed_sp23 = set()
    allowed_sp24 = set()

    try:
        wb = load_workbook(xlsx_path, read_only=True, data_only=True)
    except Exception as e:
        print("Warning: could not open sample_keys for Shared universe: {0}".format(e))
        return None

    try:
        if "initial_set" in wb.sheetnames:
            ws = wb["initial_set"]
            rows = list(ws.iter_rows(values_only=True))
            if rows:
                headers = [_normalize_header_key(c) for c in rows[0]]
                try:
                    ix_folder = headers.index("folder")
                    ix_ref = headers.index("reference_fasta")
                except ValueError:
                    ix_folder = ix_ref = None
                if ix_folder is not None and ix_ref is not None:
                    for i, row in enumerate(rows[1:], start=2):
                        if not row:
                            continue

                        def get(j):
                            return _cell_str(row[j]) if j < len(row) else ""

                        folder = get(ix_folder)
                        ref_name = get(ix_ref)
                        if not folder or not ref_name:
                            continue
                        rlow = ref_name.lower()
                        if "sp24" in rlow:
                            tag = "sp24"
                        elif "sp23" in rlow:
                            tag = "sp23"
                        else:
                            continue
                        sample_name = "initial_set_r{0}_{1}".format(
                            i, re.sub(r"[^A-Za-z0-9._-]+", "_", folder)[:80]
                        )
                        if tag == "sp23":
                            allowed_sp23.add(sample_name)
                        else:
                            allowed_sp24.add(sample_name)

        if "1-48" in wb.sheetnames:
            ws = wb["1-48"]
            col_letters = "ABCDEFGH"
            for tr in range(12):
                for tc in range(8):
                    excel_row = 2 + tc
                    excel_col = 2 + tr
                    cell = ws.cell(row=excel_row, column=excel_col)
                    text = _cell_str(cell.value)
                    if text != SHARED_UNIVERSE_MARKER_1_48:
                        continue
                    sample_num = tr * 8 + tc + 1
                    coord = "{0}{1}".format(col_letters[tc], tr + 1)
                    sample_name = "1-48_{0}_n{1}".format(coord, sample_num)
                    allowed_sp23.add(sample_name)
    finally:
        try:
            wb.close()
        except Exception:
            pass

    print(
        "Shared-color universe from sample_keys: {0} sp23 sample names, {1} sp24 sample names".format(
            len(allowed_sp23), len(allowed_sp24)
        )
    )
    return allowed_sp23, allowed_sp24


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

def aa_contributing_dna_distance(ref_dna, query_dna):
    """
    Count nucleotide differences only in codons that change the translated amino acid.

    This ignores synonymous DNA changes and changes outside codons that alter AA identity.
    Returns an integer distance.
    """
    ref = "".join(c for c in str(ref_dna).upper() if c in "ATCG")
    qry = "".join(c for c in str(query_dna).upper() if c in "ATCG")
    n = min(len(ref), len(qry))
    n -= (n % 3)
    if n <= 0:
        return 0

    dist = 0
    for i in range(0, n, 3):
        rc = ref[i:i+3]
        qc = qry[i:i+3]
        if len(rc) < 3 or len(qc) < 3:
            continue
        # Translate codon-by-codon to decide if it contributes to AA mutation
        ra = str(Seq(rc).translate())
        qa = str(Seq(qc).translate())
        if ra != qa:
            dist += sum(1 for a, b in zip(rc, qc) if a != b)
    return dist

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
    
    plt.xlabel('Expected frequency (independence)', fontsize=AXIS_LABEL_FONTSIZE)
    plt.ylabel('Observed frequency', fontsize=AXIS_LABEL_FONTSIZE)
    # No plot title (requested)
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

def shannon_entropy_from_counts(counts):
    """Shannon entropy (bits) from a list of counts."""
    total = float(sum(counts))
    if total <= 0:
        return 0.0
    ps = [c / total for c in counts if c > 0]
    return float(-sum(p * np.log2(p) for p in ps))


def _reference_type_from_path(reference_fasta):
    base = os.path.basename(reference_fasta).lower()
    if "sp23" in base:
        return "sp23"
    if "sp24" in base:
        return "sp24"
    return "unknown"


def _load_other_sample_sets(
    output_parent_dir,
    exclude_dir,
    allowed_sp23_names=None,
    allowed_sp24_names=None,
):
    """
    Load mutation sets from other samples, split by reference type.

    Expects each sample directory to contain:
    - mutation_table.csv
    - reference_used.txt (single line containing reference fasta basename or path)

    If allowed_sp23_names / allowed_sp24_names are provided (not None), only directories whose
    result tag matches that universe are used (from sample_keys.xlsx rules). None = no restriction.
    """
    other_sp23 = set()
    other_sp24 = set()

    if not output_parent_dir or not os.path.isdir(output_parent_dir):
        return other_sp23, other_sp24

    for name in sorted(os.listdir(output_parent_dir)):
        if not name.startswith("aa_analysis_results_"):
            continue
        sample_dir = os.path.join(output_parent_dir, name)
        if not os.path.isdir(sample_dir):
            continue
        if os.path.abspath(sample_dir) == os.path.abspath(exclude_dir):
            continue

        m_tag = re.match(r"^aa_analysis_results_(.+)_cluster_1based$", name)
        sample_tag = m_tag.group(1) if m_tag else None

        ref_used_file = os.path.join(sample_dir, "reference_used.txt")
        mut_file = os.path.join(sample_dir, "mutation_table.csv")
        if not os.path.exists(ref_used_file) or not os.path.exists(mut_file):
            continue

        try:
            with open(ref_used_file, "r", encoding="utf-8", errors="replace") as f:
                ref_used = f.read().strip()
            rtype = _reference_type_from_path(ref_used)

            if rtype == "sp23" and allowed_sp23_names is not None:
                if sample_tag is None or sample_tag not in allowed_sp23_names:
                    continue
            if rtype == "sp24" and allowed_sp24_names is not None:
                if sample_tag is None or sample_tag not in allowed_sp24_names:
                    continue

            df = pd.read_csv(mut_file)
            # Presence/absence is based on *any* nonzero count in other-reference samples.
            # Exclude No_mutations from "mutation set" occurrence comparisons.
            sample_sets = set()
            for _, row in df.iterrows():
                if row.get("Mutation_Set") == NO_MUTATIONS_LABEL:
                    continue
                try:
                    cnt = int(row.get("Count", 0))
                except Exception:
                    continue
                if cnt <= 0:
                    continue

                mutation_str = str(row.get("Mutation_Set", "")).strip()
                if not mutation_str or mutation_str == "nan":
                    continue
                sample_sets.add(mutation_str)

            if rtype == "sp23":
                other_sp23.update(sample_sets)
            elif rtype == "sp24":
                other_sp24.update(sample_sets)
        except Exception as e:
            print("Warning: could not load comparison data from {0}: {1}".format(sample_dir, e))

    return other_sp23, other_sp24


def plot_mutation_set_distribution(
    mutation_counts,
    no_mutations_count,
    output_dir,
    reference_fasta,
    mutation_set_distance_counts=None,
    other_sample_dirs=None,
    min_count=100,
    sample_keys_xlsx=None,
):
    """
    Plot distribution of mutation set counts with >min_count reads.

    - Includes the No_mutations count as \"starting sequence\".
    - Bars are light gray (\"Shared\") only if the mutation set appears (any count) in the opposite
      reference type among samples listed in sample_keys: all of sheet \"initial_set\", plus 1-48
      cells whose text is exactly \"sp23-2 last\" (sp23). If sample_keys is missing/unreadable,
      all sibling result directories are used (legacy behavior).
    - If reference is sp23: bars for sets not seen in the restricted sp24 universe are green (\"D3-specific\").
    - If reference is sp24: bars for sets not seen in the restricted sp23 universe are dark gray (\"WT-specific\").
    - Returns summary dicts for downstream Excel reporting.
    """
    print("Creating mutation set distribution plot...")
    
    # Convert tuple-form mutation sets into human-readable strings (used for cross-sample comparisons)
    tuple_to_str = {}
    for mset in mutation_counts.keys():
        mutation_str = "; ".join(["{0}{1}{2}".format(ref_aa, pos, mut_aa) for pos, ref_aa, mut_aa in mset])
        tuple_to_str[mset] = mutation_str

    # Filter mutation sets with >min_count reads (excluding No_mutations, which is handled separately)
    high_count_mutations = {tuple_to_str[mset]: int(count) for mset, count in mutation_counts.items() if int(count) > min_count}
    
    if not high_count_mutations and (no_mutations_count is None or int(no_mutations_count) <= 0):
        print("No mutation sets found with >{0} reads and no starting-sequence reads".format(min_count))
        return
    
    # Determine reference type and comparison universe (sp23 vs sp24 groups)
    ref_type = _reference_type_from_path(reference_fasta)

    # Default: compare against sibling result directories.
    output_parent_dir = os.path.dirname(os.path.abspath(output_dir))
    keys_path = sample_keys_xlsx or os.path.join(output_parent_dir, "sample_keys.xlsx")
    universe = _sharing_universe_sample_names_from_sample_keys(output_parent_dir, keys_path)
    if universe is None:
        print(
            "Shared universe: unrestricted (sample_keys missing or unreadable at {0})".format(keys_path)
        )
        allowed_sp23 = None
        allowed_sp24 = None
    else:
        allowed_sp23, allowed_sp24 = universe

    sp23_sets, sp24_sets = _load_other_sample_sets(
        output_parent_dir, output_dir, allowed_sp23, allowed_sp24
    )

    # Optionally extend comparison with explicitly provided directories (legacy behavior)
    # Any directory in other_sample_dirs is treated as "other sample", but we cannot infer its reference type reliably
    # without reference_used.txt. If present, it will be used; otherwise it is ignored.
    if other_sample_dirs:
        for d in other_sample_dirs:
            try:
                extra_sp23, extra_sp24 = _load_other_sample_sets(
                    os.path.dirname(d), output_dir, allowed_sp23, allowed_sp24
                )
                sp23_sets.update(extra_sp23)
                sp24_sets.update(extra_sp24)
            except Exception:
                pass

    # Build the plotted list (include starting sequence first)
    plotted_labels = []
    plotted_counts = []
    plotted_categories = []  # starting_sequence / shared / D3-specific / WT-specific

    starting_count = int(no_mutations_count or 0)
    if starting_count > 0:
        plotted_labels.append("starting sequence")
        plotted_counts.append(starting_count)
        plotted_categories.append("starting_sequence")

    # Remaining high-count mutation sets sorted by count
    sorted_mutations = sorted(high_count_mutations.items(), key=lambda x: x[1], reverse=True)
    for label, count in sorted_mutations:
        plotted_labels.append(label)
        plotted_counts.append(int(count))

        if ref_type == "sp23":
            plotted_categories.append("D3-specific" if label not in sp24_sets else "Shared")
        elif ref_type == "sp24":
            plotted_categories.append("WT-specific" if label not in sp23_sets else "Shared")
        else:
            plotted_categories.append("Shared")
    
    # Create plot
    plt.figure(figsize=(15, 8))

    # Determine colors
    shared_color = "#BFC2BB"      # light gray
    d3_color = "#118040"          # green
    wt_color = "#464747"          # dark gray
    start_color = "#eeeeee"       # very light gray for starting sequence

    colors = []
    hatches = []
    for cat in plotted_categories:
        if cat == "starting_sequence":
            colors.append(start_color)
            hatches.append("//")
        elif cat == "D3-specific":
            colors.append(d3_color)
            hatches.append("")
        elif cat == "WT-specific":
            colors.append(wt_color)
            hatches.append("")
        else:
            colors.append(shared_color)
            hatches.append("")

    bars = plt.bar(range(len(plotted_labels)), plotted_counts, color=colors, alpha=0.85, edgecolor='black', linewidth=0.5)
    for bar, hatch in zip(bars, hatches):
        if hatch:
            bar.set_hatch(hatch)
    
    # Customize plot
    plt.xlabel('Mutation set', fontsize=AXIS_LABEL_FONTSIZE)
    plt.ylabel('Read count', fontsize=AXIS_LABEL_FONTSIZE)
    if ref_type == "sp23":
        subtitle = "Green=D3-specific (absent from restricted sp24 universe), Light gray=Shared"
    elif ref_type == "sp24":
        subtitle = "Dark gray=WT-specific (absent from restricted sp23 universe), Light gray=Shared"
    else:
        subtitle = "Light gray=Shared"
    # No plot title (requested)
    plt.xticks(range(len(plotted_labels)), plotted_labels, rotation=45, ha='right', fontsize=NON_AXIS_FONTSIZE)
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add legend
    from matplotlib.patches import Patch
    set_category_counts = Counter(plotted_categories)
    legend_elements = [Patch(facecolor=shared_color, alpha=0.85, label="Shared ({0})".format(set_category_counts.get("Shared", 0)))]
    if ref_type == "sp23":
        legend_elements.append(Patch(facecolor=d3_color, alpha=0.85, label="D3-specific ({0})".format(set_category_counts.get("D3-specific", 0))))
    elif ref_type == "sp24":
        legend_elements.append(Patch(facecolor=wt_color, alpha=0.85, label="WT-specific ({0})".format(set_category_counts.get("WT-specific", 0))))
    if starting_count > 0:
        legend_elements.append(Patch(facecolor=start_color, alpha=0.85, label="Starting sequence ({0})".format(set_category_counts.get("starting_sequence", 0))))
    plt.legend(handles=legend_elements, loc='upper right')
    
    # Add count annotations on bars
    for bar, count in zip(bars, plotted_counts):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(plotted_counts) * 0.01,
            str(count),
            ha='center',
            va='bottom',
            fontsize=NON_AXIS_FONTSIZE,
        )
    
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
    
    # Summaries for Excel reporting
    total_plotted_reads = int(sum(plotted_counts))
    category_counts = Counter()
    for cat, cnt in zip(plotted_categories, plotted_counts):
        category_counts[cat] += int(cnt)
    category_fractions = {k: (v / total_plotted_reads if total_plotted_reads > 0 else 0.0) for k, v in category_counts.items()}
    entropy_bits = shannon_entropy_from_counts(plotted_counts)

    print("Mutation set distribution plot complete: {0} plotted sets, {1} plotted reads".format(len(plotted_counts), total_plotted_reads))

    # DNA mutation distance distribution (for plotted sets)
    distance_summary = None
    if mutation_set_distance_counts is not None:
        dist_counts = Counter()
        # Map label -> distance counts
        # We stored mutation_set_distance_counts keyed by label strings including NO_MUTATIONS_LABEL
        for label, cnt in zip(plotted_labels, plotted_counts):
            if label == "starting sequence":
                key = NO_MUTATIONS_LABEL
            else:
                key = label
            if key in mutation_set_distance_counts:
                dist_counts.update(mutation_set_distance_counts[key])
        if dist_counts:
            plt.figure(figsize=(8, 4.5))
            xs = sorted(dist_counts.keys())
            ys = [dist_counts[x] for x in xs]
            plt.bar(xs, ys, color="#6baed6", edgecolor="black", linewidth=0.5)
            plt.xlabel("Mutation distance", fontsize=AXIS_LABEL_FONTSIZE)
            plt.ylabel("Read count", fontsize=AXIS_LABEL_FONTSIZE)
            ax_dist = plt.gca()
            ax_dist.xaxis.set_major_locator(mticker.MultipleLocator(1))
            ax_dist.xaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
            # No plot title (requested)
            plt.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            out_pdf = os.path.join(output_dir, "dna_mutation_distance_distribution.pdf")
            try:
                plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
            except Exception:
                out_png = os.path.join(output_dir, "dna_mutation_distance_distribution.png")
                plt.savefig(out_png, dpi=300, bbox_inches="tight")
            finally:
                plt.close()
            distance_summary = {"distance_counts": dict(dist_counts)}

    return {
        "ref_type": ref_type,
        "min_count": int(min_count),
        "plotted_labels": plotted_labels,
        "plotted_counts": plotted_counts,
        "plotted_categories": plotted_categories,
        "category_counts": dict(category_counts),
        "category_fractions": category_fractions,
        "entropy_bits": entropy_bits,
        "distance_summary": distance_summary,
    }

def analyze_mutations(
    input_fasta,
    reference_fasta,
    min_length,
    gene_start,
    gene_length,
    output_dir,
    frequency_cutoff=0.05,
    other_sample_dirs=None,
    sample_keys_xlsx=None,
):
    """Main analysis function."""
    print("Loading reference from {0}".format(reference_fasta))
    reference_record = next(SeqIO.parse(reference_fasta, "fasta"))
    reference_seq = str(reference_record.seq).upper()
    print("Reference length: {0}".format(len(reference_seq)))
    
    # Load the other reference for cross-reference checking
    reference_dir = os.path.dirname(os.path.abspath(reference_fasta)) or "."
    other_reference_name = "sp24_dna.fa" if "sp23" in os.path.basename(reference_fasta).lower() else "sp23_dna.fa"
    other_reference_fasta = os.path.join(reference_dir, other_reference_name)
    print("Loading other reference from {0} for cross-reference checking".format(other_reference_fasta))
    if not os.path.exists(other_reference_fasta):
        raise FileNotFoundError(
            "Other reference FASTA not found: {0} (expected alongside --reference_fasta in {1})".format(
                other_reference_fasta, reference_dir
            )
        )
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
    mutation_set_distance_counts = defaultdict(Counter)  # label string -> Counter(distance -> count)
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

        # DNA-level mutation distance that *contributes* to AA changes.
        # If there are no AA mutations, force distance=0 so the 0-distance bin equals starting sequence count.
        if mutations:
            dna_distance = aa_contributing_dna_distance(ref_aligned, subread_aligned)
        else:
            dna_distance = 0
        
        if mutations:
            all_mutations.append(mutations)
            # Create mutation set string - mutations are already in amino acid format
            mutation_set = tuple(sorted(mutations))
            mutation_counts[mutation_set] += 1
            mut_label = "; ".join(["{0}{1}{2}".format(ref_aa, pos, mut_aa) for pos, ref_aa, mut_aa in mutation_set])
            mutation_set_distance_counts[mut_label][dna_distance] += 1
        else:
            no_mutations_count += 1
            mutation_set_distance_counts[NO_MUTATIONS_LABEL][dna_distance] += 1
    
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

    # Record reference used for cross-sample comparison
    try:
        with open(os.path.join(output_dir, "reference_used.txt"), "w", encoding="utf-8") as f:
            f.write(str(reference_fasta))
            f.write("\n")
    except Exception as e:
        print("Warning: could not write reference_used.txt: {0}".format(e))
    
    # Create mutation set distribution plot
    dist_summary = plot_mutation_set_distribution(
        mutation_counts,
        no_mutations_count,
        output_dir,
        reference_fasta,
        mutation_set_distance_counts=mutation_set_distance_counts,
        other_sample_dirs=other_sample_dirs,
        min_count=100,
        sample_keys_xlsx=sample_keys_xlsx,
    )

    # Write per-sample Excel summary
    try:
        xlsx_file = os.path.join(output_dir, "sample_results.xlsx")
        with pd.ExcelWriter(xlsx_file, engine="openpyxl") as writer:
            if dist_summary:
                # Category table
                cat_rows = []
                total_plotted = int(sum(dist_summary["plotted_counts"]))
                for cat, cnt in dist_summary["category_counts"].items():
                    cat_rows.append(
                        {
                            "sample_output_dir": os.path.basename(output_dir),
                            "reference_type": dist_summary["ref_type"],
                            "min_count": dist_summary["min_count"],
                            "category": cat,
                            "read_count": int(cnt),
                            "fraction_of_plotted_reads": float(dist_summary["category_fractions"].get(cat, 0.0)),
                            "total_plotted_reads": total_plotted,
                            "shannon_entropy_bits_plotted_sets": float(dist_summary["entropy_bits"]),
                        }
                    )
                pd.DataFrame(cat_rows).to_excel(writer, sheet_name="mutation_set_categories", index=False)

                # Plotted set table
                pd.DataFrame(
                    {
                        "mutation_set_label": dist_summary["plotted_labels"],
                        "read_count": dist_summary["plotted_counts"],
                        "category": dist_summary["plotted_categories"],
                    }
                ).to_excel(writer, sheet_name="plotted_mutation_sets", index=False)

                # DNA distance distribution table
                if dist_summary.get("distance_summary") and dist_summary["distance_summary"].get("distance_counts"):
                    dcounts = dist_summary["distance_summary"]["distance_counts"]
                    xs = sorted(dcounts.keys())
                    total = sum(dcounts.values())
                    df_dist = pd.DataFrame(
                        {
                            "dna_distance": xs,
                            "read_count": [dcounts[x] for x in xs],
                            "fraction": [(dcounts[x] / total) if total > 0 else 0.0 for x in xs],
                        }
                    )
                    df_dist.to_excel(writer, sheet_name="dna_distance_distribution", index=False)
        print("Saved per-sample Excel summary to {0}".format(xlsx_file))
    except Exception as e:
        print("Warning: could not write sample_results.xlsx: {0}".format(e))
    
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
    parser.add_argument(
        '--sample_keys',
        default=None,
        help='sample_keys.xlsx path for Shared mutation-set coloring (default: <output_dir_parent>/sample_keys.xlsx)',
    )
    
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
        args.other_sample_dirs,
        sample_keys_xlsx=args.sample_keys,
    )

if __name__ == "__main__":
    main()
