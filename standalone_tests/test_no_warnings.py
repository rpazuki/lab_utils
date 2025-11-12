"""Test that warnings are removed by default but show in verbose mode"""

import pandas as pd

from labUtils.amn_mappings import (create_supplement_exchange_matrix,
                                   load_default_iml1515_mapping)

# Create sample data with typos/variations
df = pd.DataFrame({
    'well': ['A1', 'A2', 'A3'],
    'supplements': ['Glucose', 'Glu; Adenin', 'Riboze'],
    'mu_max': [0.2, 0.15, 0.18],
    'success': [True, True, True]
})

# Get mapping
mapping = load_default_iml1515_mapping()

print("Test 1: Without verbose (should see no fuzzy match warnings)")
print("=" * 70)
matrix1 = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max',
    fuzzy_threshold=0.6,
    verbose=False  # Explicitly False
)
print("Matrix created successfully (no warnings shown)")
print(matrix1)

print("\n\n")

print("Test 2: With verbose=True (should show fuzzy match details)")
print("=" * 70)
matrix2 = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max',
    fuzzy_threshold=0.6,
    verbose=True
)
print("\nMatrix created with verbose output")
