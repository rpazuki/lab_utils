"""
Test the calculate_replicate_statistics function with sample data
"""
from pathlib import Path

from labUtils.media_bot import calculate_replicate_statistics, parse

# Paths to test data
tests_dir = Path(__file__).parent
raw_data_path = tests_dir / "rw_data.csv"
meta_data_path = tests_dir / "protocol_metadata.csv"

# Parse the data
df = parse(raw_data_path, meta_data_path)

print("Original data shape:", df.shape)
print("\nUnique wells:", sorted(df["well"].unique()))
print("\nSample of original data:")
print(df[["well", "well_row", "well_col", "time_h", "od"]].head(20))

# Test 1: Alphabetical direction with sample_size=3
print("\n" + "="*80)
print("Test 1: Alphabetical direction, sample_size=3")
print("Expected groups: A1-A2, B1-B2, C1-C2 (wells along rows)")
print("="*80)
try:
    stats_alpha = calculate_replicate_statistics(
        df,
        direction="alphabetical",
        sample_size=2,  # Changed to 2 since we have 2 columns per row
        ddof=1  # Unbiased std
    )
    print("\nResult shape:", stats_alpha.shape)
    print("\nGroup IDs:", sorted(stats_alpha["group_id"].unique()))
    print("\nSample of averaged data:")
    print(stats_alpha[["group_id",
                       "wells",
                       "well_rows",
                       "well_cols",
                       "time_h",
                       "od_mean",
                       "od_std",
                       "n_replicates"]].head(10))
except (ValueError, KeyError, FileNotFoundError, RuntimeError) as e:
    print(f"Error: {e}")

# Test 2: Numerical direction with sample_size=3
print("\n" + "="*80)
print("Test 2: Numerical direction, sample_size=3")
print("Expected groups: A1-B1-C1, A2-B2-C2 (wells along columns)")
print("="*80)
try:
    stats_num = calculate_replicate_statistics(
        df,
        direction="numerical",
        sample_size=3,
        ddof=1  # Unbiased std
    )
    print("\nResult shape:", stats_num.shape)
    print("\nGroup IDs:", sorted(stats_num["group_id"].unique()))
    print("\nSample of averaged data:")
    print(stats_num[["group_id",
                     "wells",
                     "well_rows",
                     "well_cols",
                     "time_h",
                     "od_mean",
                     "od_std",
                     "n_replicates"]].head(10))
except (ValueError, KeyError, FileNotFoundError, RuntimeError) as e:
    print(f"Error: {e}")

# Test 3: Population std (ddof=0)
print("\n" + "="*80)
print("Test 3: Same as Test 2 but with population std (ddof=0)")
print("="*80)
try:
    stats_pop = calculate_replicate_statistics(
        df,
        direction="numerical",
        sample_size=3,
        ddof=0  # Biased/population std
    )
    print("\nComparing std values (ddof=1 vs ddof=0) for first time point:")
    comparison = stats_num[stats_num["time_h"] == 0.0][["group_id", "od_std"]].merge(
        stats_pop[stats_pop["time_h"] == 0.0][["group_id", "od_std"]],
        on="group_id",
        suffixes=("_unbiased", "_biased")
    )
    print(comparison)
except (ValueError, KeyError, FileNotFoundError, RuntimeError) as e:
    print(f"Error: {e}")

print("\n" + "="*80)
print("All tests completed!")
print("="*80)

# Test 4: Using custom value_column_name
print("\n" + "="*80)
print("Test 4: Custom value_column_name compatibility test")
print("="*80)
try:
    # Parse with custom column name
    df_custom = parse(raw_data_path, meta_data_path, value_column_name="absorbance")
    print("\nParsed with custom column name: 'absorbance'")
    print(f"Columns: {list(df_custom.columns)}")

    # Calculate statistics with matching column name
    stats_custom = calculate_replicate_statistics(
        df_custom,
        direction="numerical",
        sample_size=3,
        value_column_name="absorbance"  # Must match parse()
    )
    print("\nStatistics calculated successfully!")
    print(f"Output columns: {list(stats_custom.columns)}")
    print("\nSample of averaged data with custom column name:")
    print(stats_custom[["group_id", "wells", "time_h", "absorbance_mean", "absorbance_std"]].head(5))
    print("\n✓ Success: parse() and calculate_replicate_statistics() are compatible!")
except (ValueError, KeyError, FileNotFoundError, RuntimeError) as e:
    print(f"✗ Error: {e}")

print("\n" + "="*80)
print("Compatibility test completed!")
print("="*80)
