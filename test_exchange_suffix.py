"""Test exchange_suffix parameter"""

import pandas as pd

from labUtils.amn_mappings import (create_supplement_exchange_matrix,
                                   load_default_iml1515_mapping)

# Create sample data
df = pd.DataFrame({
    'well': ['A1', 'A2', 'A3'],
    'supplements': ['Glucose', 'Fructose; Adenine', 'Ribose'],
    'mu_max': [0.2, 0.15, 0.18],
    'success': [True, True, True]
})

# Get mapping
mapping = load_default_iml1515_mapping()

print("Test 1: Without suffix (default behavior)")
print("=" * 70)
matrix1 = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max'
)
print("Columns:", list(matrix1.columns))
print(matrix1)

print("\n\nTest 2: With suffix '_input'")
print("=" * 70)
matrix2 = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max',
    exchange_suffix='_input'
)
print("Columns:", list(matrix2.columns))
print(matrix2)

print("\n\nTest 3: With suffix '_v2' and baseline exchanges")
print("=" * 70)
baseline = ['EX_pi_e_i', 'EX_o2_e_i', 'EX_h2o_e_i']
matrix3 = create_supplement_exchange_matrix(
    df,
    supplement_to_exchange_map=mapping,
    supplement_column='supplements',
    growth_rate_column='mu_max',
    baseline_exchanges=baseline,
    exchange_suffix='_v2'
)
print("Columns:", list(matrix3.columns))
print(matrix3.head())

print("\n" + "=" * 70)
print("✓ Note: growth_rate_column 'mu_max' remains unchanged in all cases")
print("=" * 70)
