import pandas as pd

# Load your parsed data
df = pd.read_csv(r'h:\ROBOT_SCIENTIST\E_coli\Growth_rates\Anthony-Five-Replicates\processed\replicates\260206_EGMB_61-72_96hr\parsed_data.csv')

print("Sample of data:")
print(df[['well', 'strain']].head(20))
print()

# Check strain + well combinations
print("Unique strain-well pairs per strain:")
for strain in sorted(df['strain'].unique())[:5]:  # Just first 5 strains
    wells = sorted(df[df['strain'] == strain]['well'].unique())
    print(f"  {strain:20s} -> wells: {wells}")
print()

# Check how many unique wells per strain
print("Count of unique wells per strain:")
strain_well_counts = df.groupby('strain')['well'].nunique()
print(strain_well_counts)
print()

# Check rule matching
import re
custom_rules = {
    "EGMB": {"pattern": r"EGMB_*_R*", "direction": "alphabetical", "sample_size": 5},
    "MG": {"pattern": r"MG1655_R.*", "direction": "alphabetical", "sample_size": 5},
    "nan": {"direction": "alphabetical", "sample_size": 1},
    "BLANK": {"direction": "alphabetical", "sample_size": 1},
}

print("Which rules match which strains:")
df_strain_copy = df[['strain']].drop_duplicates().copy()
matched = set()
for rule_key, rules in custom_rules.items():
    if "pattern" in rules:
        pat = re.compile(rules["pattern"])
        matching = df['strain'].unique()
        matching = [s for s in matching if s and pat.search(str(s))]
        print(f"  Rule '{rule_key}' (pattern): {matching}")
        matched.update(matching)
    else:
        # Exact match
        if rule_key in df['strain'].values:
            print(f"  Rule '{rule_key}' (exact match): {[rule_key]}")
            matched.add(rule_key)
        else:
            print(f"  Rule '{rule_key}' (exact match): NO MATCH (available: {list(df['strain'].unique())})")

print()
print(f"Total unique strains in data: {df['strain'].nunique()}")
print(f"Total matched strains: {len(matched)}")
print(f"Unmatched strains: {set(df['strain'].unique()) - matched}")
