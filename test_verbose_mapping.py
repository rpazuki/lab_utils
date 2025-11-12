"""Test verbose mapping output"""

import pandas as pd

from labUtils.amn_mappings import (create_supplement_exchange_matrix,
                                   load_default_iml1515_mapping)

# Create sample data
df = pd.DataFrame({
    'well': ['A1', 'A2', 'A3'],
    'supplements': ['Glucose', 'Fructose; Adenine', 'Unknown; Ribose'],
    'mu_max': [0.2, 0.15, 0.18],
    'success': [True, True, True]
})

# Get mapping
mapping = load_default_iml1515_mapping()

print("Testing verbose=True:")
print("-" * 70)

# Create matrix with verbose output
matrix = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max',
    verbose=True
)

print("\nResulting matrix:")
print(matrix)
