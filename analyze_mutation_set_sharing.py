#!/usr/bin/env python3
"""
Analyze mutation set sharing between samples with different reference fastas.
This script compares mutation sets with >100 reads across samples and categorizes them
as "Starting point-specific" or "Shared" based on whether they appear in samples
with alternative reference fastas.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import argparse

def load_sample_config(config_file):
    """Load sample configuration from CSV file."""
    df = pd.read_csv(config_file)
    return df

def load_mutation_table(sample_dir, sample_name):
    """Load mutation table for a sample from its _cluster_1based folder."""
    mutation_table_path = os.path.join(sample_dir, f"aa_analysis_results_{sample_name}_cluster_1based", "mutation_table.csv")
    
    if not os.path.exists(mutation_table_path):
        print(f"Warning: Mutation table not found for {sample_name}: {mutation_table_path}")
        return None
    
    try:
        df = pd.read_csv(mutation_table_path)
        return df
    except Exception as e:
        print(f"Error loading mutation table for {sample_name}: {e}")
        return None

def get_high_count_mutation_sets(mutation_df, min_count=100):
    """Extract mutation sets with count > min_count."""
    if mutation_df is None:
        return set()
    
    # Filter for high-count mutation sets (exclude "No_mutations")
    high_count_sets = set()
    for _, row in mutation_df.iterrows():
        if row['Mutation_Set'] != 'No_mutations' and row['Count'] > min_count:
            high_count_sets.add(row['Mutation_Set'])
    
    return high_count_sets

def get_all_mutation_sets(mutation_df):
    """Extract all mutation sets (any count > 0, exclude "No_mutations")."""
    if mutation_df is None:
        return set()
    
    # Get all mutation sets with any count > 0 (exclude "No_mutations")
    all_sets = set()
    for _, row in mutation_df.iterrows():
        if row['Mutation_Set'] != 'No_mutations' and row['Count'] > 0:
            all_sets.add(row['Mutation_Set'])
    
    return all_sets

def categorize_mutation_sets(sample_config, base_dir):
    """Categorize mutation sets as Starting point-specific or Shared."""
    results = {}
    
    # Group samples by reference fasta
    sp23_samples = sample_config[sample_config['reference_fasta'] == 'sp23_dna.fa']['sample_name'].tolist()
    sp24_samples = sample_config[sample_config['reference_fasta'] == 'sp24_dna.fa']['sample_name'].tolist()
    
    print(f"SP23 samples: {sp23_samples}")
    print(f"SP24 samples: {sp24_samples}")
    
    # Process each sample
    for _, row in sample_config.iterrows():
        sample_name = row['sample_name']
        reference_fasta = row['reference_fasta']
        
        print(f"\nProcessing sample: {sample_name} (reference: {reference_fasta})")
        
        # Load mutation table for this sample
        mutation_df = load_mutation_table(base_dir, sample_name)
        if mutation_df is None:
            continue
        
        # Get high-count mutation sets for this sample
        high_count_sets = get_high_count_mutation_sets(mutation_df)
        print(f"Found {len(high_count_sets)} mutation sets with >100 reads")
        
        # Determine alternative reference samples
        if reference_fasta == 'sp23_dna.fa':
            alternative_samples = sp24_samples
        else:  # sp24_dna.fa
            alternative_samples = sp23_samples
        
        # Check sharing with alternative reference samples
        sample_results = []
        for mutation_set in high_count_sets:
            is_shared = False
            
            # Check if this mutation set appears in any alternative reference sample (any count)
            for alt_sample in alternative_samples:
                alt_mutation_df = load_mutation_table(base_dir, alt_sample)
                if alt_mutation_df is not None:
                    alt_all_sets = get_all_mutation_sets(alt_mutation_df)
                    if mutation_set in alt_all_sets:
                        is_shared = True
                        break
            
            # Categorize the mutation set
            if is_shared:
                category = "Shared"
            else:
                category = "Starting point-specific"
            
            # Get count for this mutation set
            count = mutation_df[mutation_df['Mutation_Set'] == mutation_set]['Count'].iloc[0]
            
            sample_results.append({
                'mutation_set': mutation_set,
                'count': count,
                'category': category,
                'reference_fasta': reference_fasta
            })
        
        results[sample_name] = sample_results
        print(f"Categorized {len(sample_results)} mutation sets")
    
    return results

def create_histogram(sample_name, sample_results, output_dir):
    """Create histogram for a sample's mutation sets."""
    if not sample_results:
        print(f"No mutation sets to plot for {sample_name}")
        return
    
    # Prepare data for plotting
    mutation_sets = [r['mutation_set'] for r in sample_results]
    counts = [r['count'] for r in sample_results]
    categories = [r['category'] for r in sample_results]
    reference_fasta = sample_results[0]['reference_fasta']
    
    # Determine colors based on category and reference
    colors = []
    for cat in categories:
        if cat == "Shared":
            colors.append('lightgray')
        elif cat == "Starting point-specific":
            if reference_fasta == 'sp23_dna.fa':
                colors.append('green')
            else:  # sp24_dna.fa
                colors.append('#2F2F2F')  # Very dark gray
        else:
            colors.append('blue')  # fallback
    
    # Create the plot
    plt.figure(figsize=(15, 8))
    
    # Sort by count for better visualization
    sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)
    sorted_counts = [counts[i] for i in sorted_indices]
    sorted_colors = [colors[i] for i in sorted_indices]
    sorted_sets = [mutation_sets[i] for i in sorted_indices]
    sorted_categories = [categories[i] for i in sorted_indices]
    
    # Create bars
    bars = plt.bar(range(len(sorted_counts)), sorted_counts, color=sorted_colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Customize plot
    plt.xlabel('Mutation Set', fontsize=12)
    plt.ylabel('Read Count', fontsize=12)
    plt.title(f"Mutation Sets with >100 Reads - {sample_name}\n(Reference: {reference_fasta})", fontsize=14)
    
    # Label bars with full mutation sets
    plt.xticks(range(len(sorted_sets)), sorted_sets, rotation=45, ha='right', fontsize=8)
    plt.grid(True, alpha=0.3, axis='y')
    
    # Create legend
    from matplotlib.patches import Patch
    legend_elements = []
    
    # Count categories for legend
    shared_count = sum(1 for cat in sorted_categories if cat == "Shared")
    sp23_specific_count = sum(1 for cat, ref in zip(sorted_categories, [sample_results[i]['reference_fasta'] for i in sorted_indices]) 
                             if cat == "Starting point-specific" and ref == 'sp23_dna.fa')
    sp24_specific_count = sum(1 for cat, ref in zip(sorted_categories, [sample_results[i]['reference_fasta'] for i in sorted_indices]) 
                             if cat == "Starting point-specific" and ref == 'sp24_dna.fa')
    
    if shared_count > 0:
        legend_elements.append(Patch(facecolor='lightgray', alpha=0.7, label=f"Shared ({shared_count})"))
    if sp23_specific_count > 0:
        legend_elements.append(Patch(facecolor='green', alpha=0.7, label=f"D3-specific ({sp23_specific_count})"))
    if sp24_specific_count > 0:
        legend_elements.append(Patch(facecolor='#2F2F2F', alpha=0.7, label=f"WT-specific ({sp24_specific_count})"))
    
    if legend_elements:
        plt.legend(handles=legend_elements, loc='upper right')
    
    # Add count annotations on bars
    for i, (bar, count) in enumerate(zip(bars, sorted_counts)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(sorted_counts)*0.01, 
                str(count), ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    # Save plots
    cluster_dir = os.path.join(output_dir, f"aa_analysis_results_{sample_name}_cluster_1based")
    os.makedirs(cluster_dir, exist_ok=True)
    
    # Save as PNG
    png_file = os.path.join(cluster_dir, "mutation_set_sharing_histogram.png")
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    print(f"Saved PNG histogram to {png_file}")
    
    # Save as PDF
    pdf_file = os.path.join(cluster_dir, "mutation_set_sharing_histogram.pdf")
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    print(f"Saved PDF histogram to {pdf_file}")
    
    plt.close()
    
    # Save detailed results to CSV
    results_df = pd.DataFrame([
        {
            'mutation_set': sorted_sets[i],
            'count': sorted_counts[i],
            'category': sorted_categories[i],
            'reference_fasta': reference_fasta
        }
        for i in range(len(sorted_sets))
    ])
    
    csv_file = os.path.join(cluster_dir, "mutation_set_sharing_analysis.csv")
    results_df.to_csv(csv_file, index=False)
    print(f"Saved detailed results to {csv_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze mutation set sharing between samples with different reference fastas')
    parser.add_argument('--config_file', default='sample_config.csv', help='Sample configuration CSV file')
    parser.add_argument('--base_dir', default='.', help='Base directory containing analysis results')
    parser.add_argument('--min_count', type=int, default=100, help='Minimum read count threshold')
    
    args = parser.parse_args()
    
    print("Loading sample configuration...")
    sample_config = load_sample_config(args.config_file)
    print(f"Loaded {len(sample_config)} samples")
    
    print("\nCategorizing mutation sets...")
    results = categorize_mutation_sets(sample_config, args.base_dir)
    
    print("\nCreating histograms...")
    for sample_name, sample_results in results.items():
        if sample_results:
            print(f"\nCreating histogram for {sample_name}...")
            create_histogram(sample_name, sample_results, args.base_dir)
        else:
            print(f"No mutation sets found for {sample_name}")
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
