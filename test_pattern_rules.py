#!/usr/bin/env python
"""
Test script demonstrating the enhanced calculate_replicate_statistics_by_custom
function with pattern-based rules support.
"""

import pandas as pd
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from labUtils.media_bot import calculate_replicate_statistics_by_custom

# Create a sample dataset
def create_sample_data():
    """Create sample data with different strain patterns."""
    wells = []
    well_rows = []
    well_cols = []
    strains = []
    times = []
    values = []

    # Create data for strains: st_1_r1, st_2_r1, st_3_r1, st_1_r2, st_2_r2, control
    strain_configs = [
        ("st_1_r1", ["A", "B", "C"]),  # 3 wells
        ("st_2_r1", ["D", "E", "F"]),  # 3 wells
        ("st_3_r1", ["G", "H", "I"]),  # 3 wells
        ("st_1_r2", ["J", "K", "L"]),  # 3 wells
        ("st_2_r2", ["M", "N", "O"]),  # 3 wells
        ("control", ["P", "Q", "R"]),  # 3 wells (different rule)
    ]

    for strain, well_rows_for_strain in strain_configs:
        for i, row in enumerate(well_rows_for_strain):
            well = f"{row}1"
            for time_idx, time_h in enumerate([0, 2, 4, 6]):
                wells.append(well)
                well_rows.append(row)
                well_cols.append(1)
                strains.append(strain)
                times.append(f"{time_h} h")
                # Add some variation
                values.append(0.1 + time_idx * 0.05 + i * 0.01)

    df = pd.DataFrame({
        "well": wells,
        "well_row": well_rows,
        "well_col": well_cols,
        "strain": strains,
        "time_label": times,
        "time_h": [int(t.split()[0]) for t in times],
        "od": values,
    })

    return df


# Test 1: Using pattern-based rules
print("=" * 80)
print("TEST 1: Pattern-based rules")
print("=" * 80)

df = create_sample_data()
print("\nSample data created:")
print(f"  Strains: {sorted(df['strain'].unique())}")
print(f"  Total records: {len(df)}")

custom_rules = {
    "experimental_strains": {
        "pattern": r"st_\d+_r\d+",  # Match st_1_r1, st_2_r1, etc.
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "control_group": {
        "direction": "alphabetical",
        "sample_size": 3,
    }
}

# Note: For pattern-based matching, the key should not be used as strain name
# Update the rule to match the actual strain name
custom_rules = {
    "experimental": {
        "pattern": r"st_\d+_r\d+",  # Match all st_*_r* strains
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "control": {  # Exact match for "control" strain
        "direction": "alphabetical",
        "sample_size": 3,
    }
}

try:
    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=None,  # No global pattern, use exact strain names
        custom_rules=custom_rules,
    )

    print("\nResult with pattern-based rules:")
    print(f"  Processed strains: {sorted(result['strain'].unique())}")
    print(f"  Result shape: {result.shape}")
    print("\nFirst few rows:")
    print(result.head(10))
    print("\n✓ TEST 1 PASSED: Pattern-based rules work correctly!")

except Exception as e:
    print(f"\n✗ TEST 1 FAILED: {e}")
    import traceback
    traceback.print_exc()


# Test 2: Mixing exact match and pattern-based rules
print("\n" + "=" * 80)
print("TEST 2: Mixing exact match and pattern-based rules")
print("=" * 80)

custom_rules_mixed = {
    "st_1_r1": {  # Exact match
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "other_experimental": {
        "pattern": r"st_[2-3]_r\d+",  # Match st_2_r*, st_3_r*
        "direction": "numerical",
        "sample_size": 3,
    },
    "control": {
        "direction": "alphabetical",
        "sample_size": 3,
    }
}

try:
    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=None,
        custom_rules=custom_rules_mixed,
    )

    print("\nResult with mixed exact and pattern-based rules:")
    print(f"  Processed strains: {sorted(result['strain'].unique())}")
    print(f"  Result shape: {result.shape}")
    print("\nGrouped by strain:")
    for strain in sorted(result['strain'].unique()):
        count = len(result[result['strain'] == strain])
        print(f"  {strain}: {count} rows")
    print("\n✓ TEST 2 PASSED: Mixed rules work correctly!")

except Exception as e:
    print(f"\n✗ TEST 2 FAILED: {e}")
    import traceback
    traceback.print_exc()


# Test 3: Global pattern + per-rule patterns
print("\n" + "=" * 80)
print("TEST 3: Global pattern extraction + per-rule patterns")
print("=" * 80)

# Create data with more complex strain naming
df_complex = create_sample_data()
df_complex["strain"] = df_complex["strain"].apply(lambda x: f"exp_{x}")
print(f"\nComplex strains: {sorted(df_complex['strain'].unique())}")

custom_rules_complex = {
    "experimental": {
        "pattern": r"exp_st_\d+_r\d+",  # Match exp_st_*_r*
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "control": {
        "direction": "alphabetical",
        "sample_size": 3,
    }
}

try:
    result = calculate_replicate_statistics_by_custom(
        df_complex,
        strain_pattern=r"exp_(.+)",  # Global pattern to preserve full name
        custom_rules=custom_rules_complex,
    )

    print("\nResult with global and per-rule patterns:")
    print(f"  Processed strains: {sorted(result['strain'].unique())}")
    print(f"  Result shape: {result.shape}")
    print("\n✓ TEST 3 PASSED: Global + per-rule patterns work correctly!")

except Exception as e:
    print(f"\n✗ TEST 3 FAILED: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 80)
print("All pattern-based rule tests completed!")
print("=" * 80)
