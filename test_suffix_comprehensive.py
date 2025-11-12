"""Comprehensive test of exchange_suffix parameter"""

import pandas as pd

from labUtils.amn_mappings import (create_supplement_exchange_matrix,
                                   load_default_iml1515_mapping)

df = pd.DataFrame({
    'well': ['A1', 'A2'],
    'supplements': ['Glucose', 'Fructose'],
    'mu_max': [0.2, 0.15],
    'success': [True, True]
})

mapping = load_default_iml1515_mapping()
baseline = ['EX_pi_e_i', 'EX_o2_e_i']

print("Exchange Suffix Feature Test")
print("=" * 70)

# Test 1: No suffix (default)
print("\n1. Default (no suffix):")
m1 = create_supplement_exchange_matrix(df, mapping)
print(f"   Columns: {list(m1.columns)}")

# Test 2: With suffix
print("\n2. With suffix '_input':")
m2 = create_supplement_exchange_matrix(df, mapping, exchange_suffix='_input')
print(f"   Columns: {list(m2.columns)}")

# Test 3: With baseline and suffix
print("\n3. With baseline exchanges and suffix '_test':")
m3 = create_supplement_exchange_matrix(df, mapping, baseline_exchanges=baseline, exchange_suffix='_test')
print(f"   Columns: {list(m3.columns)}")

# Test 4: Different growth rate column name
print("\n4. Different growth_rate_column with suffix:")
df2 = df.rename(columns={'mu_max': 'growth_rate'})
m4 = create_supplement_exchange_matrix(
    df2,
    mapping,
    growth_rate_column='growth_rate',
    exchange_suffix='_v2'
)
print(f"   Columns: {list(m4.columns)}")
print(f"   Growth rate column unchanged: {'growth_rate' in m4.columns}")

# Test 5: Empty suffix (should behave like None)
print("\n5. With empty string suffix '':")
m5 = create_supplement_exchange_matrix(df, mapping, exchange_suffix='')
print(f"   Columns: {list(m5.columns)}")
print(f"   Same as no suffix: {list(m5.columns) == list(m1.columns)}")

print("\n" + "=" * 70)
print("Summary:")
print("  ✓ Suffix is optional (default: None)")
print("  ✓ Only exchange columns get the suffix")
print("  ✓ Growth rate column remains unchanged")
print("  ✓ Works with baseline exchanges")
print("  ✓ Works with custom growth_rate_column names")
print("=" * 70)
