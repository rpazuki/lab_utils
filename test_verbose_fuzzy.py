"""Test verbose mapping with fuzzy matches"""

import pandas as pd

from labUtils.amn_mappings import (create_supplement_exchange_matrix,
                                   load_default_iml1515_mapping)

# Create sample data with some typos/variations
df = pd.DataFrame({
    'well': ['A1', 'A2', 'A3'],
    'supplements': ['Glucose', 'Glu; Adenin', 'Unknown_Chemical; Riboze'],
    'mu_max': [0.2, 0.15, 0.18],
    'success': [True, True, True]
})

# Get mapping
mapping = load_default_iml1515_mapping()

print("Testing verbose=True with fuzzy matching:")
print("-" * 70)

# Create matrix with verbose output
matrix = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max',
    fuzzy_threshold=0.6,
    verbose=True
)

print("\nResulting matrix:")
print(matrix)
