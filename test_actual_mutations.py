#!/usr/bin/env python3
"""
Test script to check actual mutation data from the analysis.
"""

import pandas as pd

# Read the mutation table
df = pd.read_csv("aa_analysis_results_A1_test_10x_v2/mutation_table.csv")

print("First 10 rows of mutation table:")
print(df.head(10))

# Check if the mutations are DNA or amino acid format
print("\nChecking mutation format:")
for i, row in df.head(5).iterrows():
    mutation_set = row['Mutation_Set']
    if mutation_set != 'No_mutations':
        print(f"Row {i}: {mutation_set}")
        # Check if it looks like DNA (G485T) or amino acid (S162I)
        if 'G' in mutation_set and 'T' in mutation_set and any(char.isdigit() for char in mutation_set):
            print("  -> Looks like DNA format")
        elif any(char in mutation_set for char in 'ACDEFGHIKLMNPQRSTVWY') and any(char.isdigit() for char in mutation_set):
            print("  -> Looks like amino acid format")
        else:
            print("  -> Unknown format")
