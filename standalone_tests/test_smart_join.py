"""Test script for smart_join function."""

import numpy as np
import pandas as pd

from src.labUtils.utils import smart_join

print("Testing smart_join function\n")
print("=" * 80)

# Test 1: Identical columns (should keep only one)
print("\nTest 1: Identical columns in both DataFrames")
print("-" * 80)
df1 = pd.DataFrame({"id": [1, 2, 3], "name": ["Alice", "Bob", "Charlie"], "status": ["active", "active", "inactive"]})
df2 = pd.DataFrame(
    {
        "id": [1, 2, 3],
        "name": ["Alice", "Bob", "Charlie"],  # Same values as df1
        "score": [95, 87, 92],
    }
)

print("Left DataFrame:")
print(df1)
print("\nRight DataFrame:")
print(df2)

result1 = smart_join(df1, df2, on_cols=["id"])
print("\nResult (should have only one 'name' column):")
print(result1)
print(f"Columns: {list(result1.columns)}")

# Test 2: Different columns (should keep both with _right suffix)
print("\n\n" + "=" * 80)
print("\nTest 2: Different values in duplicate columns")
print("-" * 80)
df3 = pd.DataFrame({"id": [1, 2, 3], "status": ["active", "inactive", "pending"]})
df4 = pd.DataFrame(
    {
        "id": [1, 2, 3],
        "status": ["complete", "review", "approved"],  # Different values
    }
)

print("Left DataFrame:")
print(df3)
print("\nRight DataFrame:")
print(df4)

result2 = smart_join(df3, df4, on_cols=["id"])
print("\nResult (should have 'status' and 'status_right'):")
print(result2)
print(f"Columns: {list(result2.columns)}")

# Test 3: Mixed case - some identical, some different
print("\n\n" + "=" * 80)
print("\nTest 3: Mixed - some columns identical, some different")
print("-" * 80)
df5 = pd.DataFrame({"id": [1, 2], "city": ["NYC", "LA"], "country": ["USA", "USA"], "value_a": [100, 200]})
df6 = pd.DataFrame(
    {
        "id": [1, 2],
        "city": ["Boston", "Seattle"],  # Different
        "country": ["USA", "USA"],  # Same
        "value_b": [10, 20],
    }
)

print("Left DataFrame:")
print(df5)
print("\nRight DataFrame:")
print(df6)

result3 = smart_join(df5, df6, on_cols=["id"])
print("\nResult ('country' once, 'city' and 'city_right'):")
print(result3)
print(f"Columns: {list(result3.columns)}")

# Test 4: NaN handling
print("\n\n" + "=" * 80)
print("\nTest 4: NaN handling (NaN == NaN should be treated as identical)")
print("-" * 80)
df7 = pd.DataFrame({"id": [1, 2, 3], "optional": ["A", np.nan, "C"]})
df8 = pd.DataFrame(
    {
        "id": [1, 2, 3],
        "optional": ["A", np.nan, "C"],  # Same, including NaN
    }
)

print("Left DataFrame:")
print(df7)
print("\nRight DataFrame:")
print(df8)

result4 = smart_join(df7, df8, on_cols=["id"])
print("\nResult (should have only one 'optional' column):")
print(result4)
print(f"Columns: {list(result4.columns)}")

# Test 5: Multiple join keys
print("\n\n" + "=" * 80)
print("\nTest 5: Multiple join keys")
print("-" * 80)
df9 = pd.DataFrame({"year": [2023, 2023, 2024], "month": [1, 2, 1], "temp": [20, 22, 18]})
df10 = pd.DataFrame(
    {
        "year": [2023, 2023, 2024],
        "month": [1, 2, 1],
        "temp": [20, 22, 18],  # Same
        "humidity": [65, 70, 60],
    }
)

print("Left DataFrame:")
print(df9)
print("\nRight DataFrame:")
print(df10)

result5 = smart_join(df9, df10, on_cols=["year", "month"])
print("\nResult (should have only one 'temp' column):")
print(result5)
print(f"Columns: {list(result5.columns)}")

print("\n\n" + "=" * 80)
print("All tests completed!")
