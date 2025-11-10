"""
Example demonstrating compatibility between parse() and calculate_replicate_statistics()
with value_column_name parameter
"""
from pathlib import Path

from labUtils.media_bot import calculate_replicate_statistics, parse

# Paths to test data
tests_dir = Path(__file__).parent
raw_data_path = tests_dir / "rw_data.csv"
meta_data_path = tests_dir / "protocol_metadata.csv"

print("="*80)
print("Compatibility Example: parse() + calculate_replicate_statistics()")
print("="*80)

# Example 1: Default value_column_name="od"
print("\n" + "="*80)
print("Example 1: Default behavior (value_column_name='od')")
print("="*80)

df1 = parse(raw_data_path, meta_data_path)  # Default: value_column_name="od"
print("✓ Parsed data with column: 'od'")

stats1 = calculate_replicate_statistics(
    df1,
    direction="numerical",
    sample_size=3
)  # Default: value_column_name="od"
print("✓ Calculated statistics with columns: 'od_mean', 'od_std'")
print(f"  Result shape: {stats1.shape}")
print(f"  Output columns: {list(stats1.columns[:7])}")

# Example 2: Custom value_column_name="od600"
print("\n" + "="*80)
print("Example 2: Using value_column_name='od600'")
print("="*80)

df2 = parse(raw_data_path, meta_data_path, value_column_name="od600")
print("✓ Parsed data with column: 'od600'")

stats2 = calculate_replicate_statistics(
    df2,
    direction="numerical",
    sample_size=3,
    value_column_name="od600"  # MUST match parse()
)
print("✓ Calculated statistics with columns: 'od600_mean', 'od600_std'")
print(f"  Result shape: {stats2.shape}")
print(f"  Output columns: {list(stats2.columns[:7])}")

# Example 3: Custom value_column_name="absorbance"
print("\n" + "="*80)
print("Example 3: Using value_column_name='absorbance'")
print("="*80)

df3 = parse(raw_data_path, meta_data_path, value_column_name="absorbance")
print("✓ Parsed data with column: 'absorbance'")

stats3 = calculate_replicate_statistics(
    df3,
    direction="numerical",
    sample_size=3,
    value_column_name="absorbance"  # MUST match parse()
)
print("✓ Calculated statistics with columns: 'absorbance_mean', 'absorbance_std'")
print(f"  Result shape: {stats3.shape}")
print(f"  Output columns: {list(stats3.columns[:7])}")

# Show actual data comparison
print("\n" + "="*80)
print("Data Comparison (first time point, first group)")
print("="*80)
print(f"\nWith 'od':         {stats1.iloc[0]['od_mean']:.4f} ± {stats1.iloc[0]['od_std']:.4f}")
print(f"With 'od600':      {stats2.iloc[0]['od600_mean']:.4f} ± {stats2.iloc[0]['od600_std']:.4f}")
print(f"With 'absorbance': {stats3.iloc[0]['absorbance_mean']:.4f} ± {stats3.iloc[0]['absorbance_std']:.4f}")
print("\n(All values are identical - only column names differ)")

# Important note
print("\n" + "="*80)
print("IMPORTANT: The value_column_name parameter MUST match between")
print("parse() and calculate_replicate_statistics() for compatibility!")
print("="*80)

# Example of incompatibility error
print("\n" + "="*80)
print("Example 4: Demonstrating incompatibility error")
print("="*80)

try:
    df_mismatch = parse(raw_data_path, meta_data_path, value_column_name="od600")
    print("✓ Parsed with value_column_name='od600'")

    # Trying to use different column name will fail
    stats_mismatch = calculate_replicate_statistics(
        df_mismatch,
        direction="numerical",
        sample_size=3,
        value_column_name="od"  # MISMATCH! DataFrame has 'od600', not 'od'
    )
except ValueError as e:
    print(f"✗ Error as expected: {e}")
    print("\nSolution: Use the same value_column_name in both functions!")

print("\n" + "="*80)
print("Best Practice:")
print("="*80)
print("""
# Method 1: Use a variable
value_col = "od600"
df = parse(raw_path, meta_path, value_column_name=value_col)
stats = calculate_replicate_statistics(df, ..., value_column_name=value_col)

# Method 2: Use default everywhere (recommended for simplicity)
df = parse(raw_path, meta_path)  # Uses "od"
stats = calculate_replicate_statistics(df, ...)  # Uses "od"
""")
