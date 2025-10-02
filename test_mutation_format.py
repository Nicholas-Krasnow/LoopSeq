#!/usr/bin/env python3
"""
Test script to verify mutation labeling format.
"""

# Test the mutation formatting
mutation_set = [(6, 'S', 'T'), (32, 'S', 'A'), (56, 'H', 'L')]

# Test the old format
old_format = "; ".join([f"{ref_aa}>{mut_aa}" for _, ref_aa, mut_aa in mutation_set])
print(f"Old format: {old_format}")

# Test the new format
new_format = "; ".join([f"{ref_aa}{pos}{mut_aa}" for pos, ref_aa, mut_aa in mutation_set])
print(f"New format: {new_format}")

# Expected format should be: S6T; S32A; H56L
print(f"Expected: S6T; S32A; H56L")
