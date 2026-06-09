import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))
from labUtils.media_bot import calculate_replicate_statistics_by_custom

# Load parsed data
df = pd.read_csv(r'h:\ROBOT_SCIENTIST\E_coli\Growth_rates\Anthony-Five-Replicates\processed\replicates\260206_EGMB_61-72_96hr\parsed_data.csv')

print("Testing corrected custom rules...")
print(f"Input shape: {df.shape}")
print(f"Unique strains: {df['strain'].nunique()}")
print()

custom_rules = {
    "EGMB": {
        "pattern": r"EGMB_\d+",
        "direction": "alphabetical",
        "sample_size": 5,
    },
    "MG": {
        "pattern": r"MG1655",
        "direction": "alphabetical",
        "sample_size": 4,
    },
    "Blank (empty_well)": {
        "direction": "alphabetical",
        "sample_size": 1,
    },
}

try:
    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=r'([A-Za-z0-9]+(?:_\d+)?)',  # Extract base strain
        custom_rules=custom_rules,
        value_column_name="od600",
    )

    print(f"✓ SUCCESS!")
    print(f"Output shape: {result.shape}")
    print(f"Unique strains in output: {sorted(result['strain'].unique())}")
    print()
    print("Output columns:", list(result.columns))
    print()
    print("First few rows:")
    print(result.head(15))

except Exception as e:
    print(f"✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
