"""
Example demonstrating the use of value_column_name parameter in parse function
"""
from pathlib import Path

from labUtils.media_bot import parse

# Paths to test data
tests_dir = Path(__file__).parent
raw_data_path = tests_dir / "rw_data.csv"
meta_data_path = tests_dir / "protocol_metadata.csv"

print("="*80)
print("Example 1: Default value_column_name (od)")
print("="*80)
# Parse with default column name "od"
df_default = parse(raw_data_path, meta_data_path)
print(f"\nColumns in DataFrame: {list(df_default.columns)}")
print("\nValue column name: 'od'")
print(df_default[["well", "time_h", "od"]].head(10))

print("\n" + "="*80)
print("Example 2: Custom value_column_name (od600)")
print("="*80)
# Parse with custom column name "od600"
df_custom = parse(raw_data_path, meta_data_path, value_column_name="od600")
print(f"\nColumns in DataFrame: {list(df_custom.columns)}")
print("\nValue column name: 'od600'")
print(df_custom[["well", "time_h", "od600"]].head(10))

print("\n" + "="*80)
print("Example 3: Another custom value_column_name (absorbance)")
print("="*80)
# Parse with another custom column name "absorbance"
df_absorbance = parse(raw_data_path, meta_data_path, value_column_name="absorbance")
print(f"\nColumns in DataFrame: {list(df_absorbance.columns)}")
print("\nValue column name: 'absorbance'")
print(df_absorbance[["well", "time_h", "absorbance"]].head(10))

print("\n" + "="*80)
print("Benefits of value_column_name parameter:")
print("="*80)
print("1. Flexibility: Use 'od', 'od600', 'absorbance', or any other name")
print("2. Clarity: Column name reflects the actual measurement type")
print("3. Compatibility: Default 'od' maintains backward compatibility")
print("4. Consistency: Same parameter across parse_raw_CLARIOstar_export and parse")
