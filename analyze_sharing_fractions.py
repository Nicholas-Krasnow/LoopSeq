#!/usr/bin/env python3
"""
Analyze mutation set sharing fractions across all _cluster_1based folders.
Calculates the fraction of reads that are Shared vs Starting point-specific.
"""

import os
import pandas as pd
import glob

def analyze_sharing_fractions(base_dir="."):
    """Analyze sharing fractions for all _cluster_1based folders."""
    
    # Find all mutation_set_sharing_analysis.csv files in _cluster_1based folders
    pattern = os.path.join(base_dir, "aa_analysis_results_*_cluster_1based", "mutation_set_sharing_analysis.csv")
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        print(f"No mutation_set_sharing_analysis.csv files found matching pattern: {pattern}")
        return None
    
    print(f"Found {len(csv_files)} mutation_set_sharing_analysis.csv files")
    
    results = []
    
    for csv_file in sorted(csv_files):
        # Extract folder name from path
        folder_path = os.path.dirname(csv_file)
        folder_name = os.path.basename(folder_path)
        
        print(f"\nAnalyzing: {folder_name}")
        
        try:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Calculate total reads
            total_reads = df['count'].sum()
            
            # Calculate reads by category
            shared_reads = df[df['category'] == 'Shared']['count'].sum()
            specific_reads = df[df['category'] == 'Starting point-specific']['count'].sum()
            
            # Calculate fractions
            fraction_shared = shared_reads / total_reads if total_reads > 0 else 0.0
            fraction_specific = specific_reads / total_reads if total_reads > 0 else 0.0
            
            # Get reference fasta (should be the same for all rows)
            reference_fasta = df['reference_fasta'].iloc[0] if len(df) > 0 else "unknown"
            
            results.append({
                'folder': folder_name,
                'reference_fasta': reference_fasta,
                'total_reads': total_reads,
                'shared_reads': shared_reads,
                'specific_reads': specific_reads,
                'fraction_shared': fraction_shared,
                'fraction_specific': fraction_specific
            })
            
            print(f"  Total reads: {total_reads:,}")
            print(f"  Shared reads: {shared_reads:,} ({fraction_shared:.1%})")
            print(f"  Starting point-specific reads: {specific_reads:,} ({fraction_specific:.1%})")
            
        except Exception as e:
            print(f"  Error reading {csv_file}: {e}")
            continue
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(results)
    
    # Save to CSV
    output_file = "mutation_sharing_fractions_summary.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\n{'='*60}")
    print(f"Summary saved to: {output_file}")
    print(f"{'='*60}")
    
    # Print summary table
    print("\nSummary Table:")
    print(summary_df.to_string(index=False))
    
    return summary_df

if __name__ == "__main__":
    analyze_sharing_fractions()
