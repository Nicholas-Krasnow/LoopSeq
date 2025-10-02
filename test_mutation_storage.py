#!/usr/bin/env python3
"""
Test script to check how mutations are stored in mutation_counts.
"""

# Simulate the mutation data storage
mutations = [(162, 'S', 'I'), (163, 'S', 'F')]
mutation_set = tuple(sorted(mutations))

print(f"Original mutations: {mutations}")
print(f"Mutation set: {mutation_set}")

# Test the formatting used in the script
mutation_str = "; ".join([f"{ref_aa}{pos}{mut_aa}" for pos, ref_aa, mut_aa in mutation_set])
print(f"Formatted mutation string: {mutation_str}")

# Test with Counter
from collections import Counter
mutation_counts = Counter()
mutation_counts[mutation_set] += 1

print(f"Mutation counts: {mutation_counts}")

# Test the most common formatting
most_common = mutation_counts.most_common(1)[0]
print(f"Most common: {most_common}")

mutation_str_common = "; ".join([f"{ref_aa}{pos}{mut_aa}" for pos, ref_aa, mut_aa in most_common[0]])
print(f"Most common formatted: {mutation_str_common}")
