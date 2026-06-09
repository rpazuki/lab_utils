import pandas as pd
import re

# Load parsed data to see actual strain values
df = pd.read_csv(r'h:\ROBOT_SCIENTIST\E_coli\Growth_rates\Anthony-Five-Replicates\processed\replicates\260206_EGMB_61-72_96hr\parsed_data.csv')

print("Unique strains in parsed_data.csv:")
print(sorted(df['strain'].unique()))
print()

# Test the patterns you provided
patterns = {
    "EGMB (user)": r"EGMB_*_R*",
    "EGMB (corrected)": r"EGMB_\d+_R\d+",
    "EGMB (alt)": r"EGMB_.*_R.*",
    "MG (user)": r"MG1655_R.*",
    "MG (alt)": r"MG\d+_R.*",
}

for strain in sorted(df['strain'].unique()):
    print(f"Strain: {strain!r}")
    for rule_name, pattern_str in patterns.items():
        try:
            pat = re.compile(pattern_str)
            match = pat.search(str(strain))
            print(f"  {rule_name:30s} {pattern_str:20s} -> {'✓ MATCH' if match else '✗ NO MATCH'}")
        except re.error as e:
            print(f"  {rule_name:30s} {pattern_str:20s} -> ERROR: {e}")
    print()
