# Metabolic Mapping Module

Convert growth rate data with supplement information into binary exchange reaction matrices for metabolic modeling.

## Overview

The `metabolic_mapping` module provides tools to:
1. **Parse SBML metabolic models** (via COBRApy) to extract exchange reactions
2. **Map supplement names to exchange reactions** (exact + fuzzy matching)
3. **Create binary matrices** suitable for machine learning with metabolic models
4. **Audit mappings** to identify unmapped supplements

## Installation

```bash
pip install cobra  # Required for SBML parsing
```

COBRApy is already installed in this workspace.

## Quick Start

```python
import pandas as pd
from labUtils.metabolic_mapping import (
    create_supplement_exchange_matrix,
    parse_sbml_exchanges,
    get_supplement_mapping,
    load_default_iml1515_mapping,
    load_minimal_media_exchanges,
)

# Your growth rate data (from fit_modified_gompertz_per_series)
growth_df = pd.DataFrame({
    'well': ['A1', 'A2', 'A3'],
    'supplements': ['Glucose', 'Maltose', 'Glucose; Ribose'],
    'mu_max': [0.2, 0.15, 0.22],
})

# Option 1: Parse SBML file first, then create matrix
sbml_mapping = parse_sbml_exchanges("path/to/iML1515.xml")
matrix = create_supplement_exchange_matrix(
    growth_df,
    supplement_to_exchange_map=sbml_mapping,
    supplement_column="supplements",
    growth_rate_column="mu_max",
    baseline_exchanges=load_minimal_media_exchanges(),
)

# Option 2: Use curated mapping (no SBML needed)
mapping = load_default_iml1515_mapping()
matrix = create_supplement_exchange_matrix(
    growth_df,
    supplement_to_exchange_map=mapping,
    baseline_exchanges=load_minimal_media_exchanges(),
)
```

## Functions

### `parse_sbml_exchanges()`

Parse SBML metabolic model file to extract exchange reactions.

**Parameters:**
- `sbml_path`: Path to SBML model file (str or Path)

**Returns:** Dict mapping metabolite names to exchange reaction IDs

**Example:**
```python
# Parse SBML first
sbml_mapping = parse_sbml_exchanges("iML1515.xml")
print(sbml_mapping["d-glucose"])  # Output: EX_glc__D_e

# Use the mapping
matrix = create_supplement_exchange_matrix(
    growth_df,
    supplement_to_exchange_map=sbml_mapping
)
```

### `create_supplement_exchange_matrix()`

Create a binary matrix mapping supplements to exchange reactions.

**Parameters:**
- `growth_rates_df`: DataFrame with growth rates and supplement info
- `supplement_to_exchange_map`: Dict mapping supplement names to exchange IDs (REQUIRED)
- `supplement_column`: Column name with supplement data (default: "supplements")
- `growth_rate_column`: Column name with growth rates (default: "mu_max")
- `baseline_exchanges`: List of always-present exchange IDs (optional)
- `separator`: Supplement separator character (default: ";")
- `fuzzy_threshold`: Fuzzy matching threshold 0-1 (default: 0.6)

**Returns:** DataFrame with exchange reactions as columns, growth rate as last column

**Example:**
```python
# Parse SBML first
mapping = parse_sbml_exchanges("iML1515.xml")

# Create matrix
matrix = create_supplement_exchange_matrix(
    growth_df,
    supplement_to_exchange_map=mapping,
    baseline_exchanges=["EX_o2_e_i", "EX_pi_e_i"],
    fuzzy_threshold=0.7
)
```

### `get_supplement_mapping()`

Audit supplement-to-reaction mappings before creating the matrix.

**Parameters:**
- `growth_rates_df`: DataFrame with supplement info
- `supplement_to_exchange_map`: Dict mapping supplement names to exchange IDs (REQUIRED)
- `supplement_column`: Column name with supplement data (default: "supplements")
- `separator`: Supplement separator character (default: ";")
- `fuzzy_threshold`: Fuzzy matching threshold 0-1 (default: 0.6)

**Returns:** Dict with three keys:
- `"exact_matches"`: Dict of exact supplement->reaction mappings
- `"fuzzy_matches"`: Dict of fuzzy supplement->reaction mappings (with match info)
- `"unmapped"`: List of supplements that couldn't be mapped

**Example:**
```python
# Parse SBML first
mapping = parse_sbml_exchanges("iML1515.xml")

# Audit mappings
mapping_info = get_supplement_mapping(
    growth_df,
    supplement_to_exchange_map=mapping,
    fuzzy_threshold=0.6
)

print("Exact matches:", mapping_info["exact_matches"])
print("Fuzzy matches:", mapping_info["fuzzy_matches"])
print("Unmapped:", mapping_info["unmapped"])

# Fix unmapped supplements by merging with custom mappings
custom_mapping = {
    **mapping,
    "custom_supplement": "EX_custom_e"
}
```

### `load_default_iml1515_mapping()`

Load curated mapping for common E. coli supplements.

**Returns:** Dict mapping supplement names to exchange reaction IDs

**Includes:**
- Sugars: glucose, fructose, galactose, ribose, maltose, trehalose, glycerol
- Amino acids: alanine, proline, threonine, glycine
- Organic acids: acetate, pyruvate, succinate, lactate
- Nucleobases: adenine, uracil, guanine, cytosine, thymine

**Example:**
```python
mapping = load_default_iml1515_mapping()
print(mapping["glucose"])  # Output: EX_glc__D_e
```

### `load_minimal_media_exchanges()`

Load standard minimal media exchange reactions (baseline components).

**Returns:** List of exchange reaction IDs for minerals, ions, gases

**Includes:** phosphate, oxygen, ammonium, sulfate, trace metals, etc. (23 total)

**Example:**
```python
baseline = load_minimal_media_exchanges()
# These will be set to 1 for all rows in the matrix
```

## Workflow Example

### Step 1: Parse SBML or Load Curated Mapping

```python
# Option A: Parse SBML file
from labUtils.metabolic_mapping import parse_sbml_exchanges
mapping = parse_sbml_exchanges("iML1515.xml")

# Option B: Use curated mapping
from labUtils.metabolic_mapping import load_default_iml1515_mapping
mapping = load_default_iml1515_mapping()
```

### Step 2: Audit Your Data

```python
# Check what will be mapped
mapping_info = get_supplement_mapping(
    growth_df,
    supplement_to_exchange_map=mapping,
    fuzzy_threshold=0.6
)

# Review results
for supp in mapping_info["unmapped"]:
    print(f"⚠️  '{supp}' needs manual mapping")
```

### Step 3: Fix Unmapped Supplements

```python
# Create custom mapping for unmapped supplements
custom_mapping = {
    **mapping,  # Start with SBML or curated mapping
    "weird_supplement": "EX_weird_e",
    "custom_carbon": "EX_custom_e"
}
```

### Step 4: Create Matrix

```python
matrix = create_supplement_exchange_matrix(
    growth_df,
    supplement_to_exchange_map=custom_mapping,
    baseline_exchanges=load_minimal_media_exchanges(),
    growth_rate_column="mu_max"
)

# Save if needed
matrix.to_csv("exchange_matrix.csv", index=False)
```

## Fuzzy Matching

The module uses fuzzy string matching when exact matches fail:

```python
# "glucose" might fuzzy-match to "d-glucose" in SBML
# Warnings show what was matched:
# UserWarning: Supplement 'glucose' fuzzy-matched to 'd-glucose' -> EX_glc__D_e
```

**Adjust threshold:**
- Lower (0.5): More permissive, more matches (higher false positive risk)
- Higher (0.8): More strict, fewer matches (may miss valid matches)
- Default (0.6): Balanced

## Output Format

The resulting matrix has:
- **Rows**: One per row in input dataframe
- **Columns**:
  - Exchange reactions (sorted alphabetically)
  - Growth rate column (last, preserves original name)
- **Values**:
  - `1` = metabolite present
  - `0` = metabolite absent

**Example output:**
```
   EX_glc__D_e  EX_malt_e_i  EX_o2_e_i  EX_pi_e_i  mu_max
0            1            0          1          1    0.20
1            0            1          1          1    0.15
2            1            0          1          1    0.18
```

## SBML Parsing

The module provides explicit SBML parsing functions:

```python
from labUtils.metabolic_mapping import parse_sbml_exchanges

# Parse SBML to get exchange mappings
mapping = parse_sbml_exchanges("path/to/iML1515.xml")

# Inspect what was found
print(f"Found {len(mapping)} metabolites")
for name, exchange_id in list(mapping.items())[:5]:
    print(f"  {name:20} -> {exchange_id}")

# Use the mapping
matrix = create_supplement_exchange_matrix(
    growth_df,
    supplement_to_exchange_map=mapping
)
```

**Two-step workflow benefits:**
1. **Transparency**: See exactly what's in the SBML before creating matrix
2. **Reusability**: Parse once, use multiple times
3. **Customization**: Merge SBML with manual mappings easily
4. **Debugging**: Inspect mappings before matrix creation

**Requirements:**
- COBRApy package installed (for primary parser)
- Fallback XML parser available if COBRApy fails
- Valid SBML file (XML format)
- Exchange reactions must have "EX_" prefix

**Downloading models:**
```python
# iML1515 from BiGG Models: http://bigg.ucsd.edu/models/iML1515
# Or use the test example which auto-downloads
```

## Examples

See the `examples/` directory:
- `example_metabolic_mapping.py` - Basic usage patterns
- `example_audit_mapping.py` - Auditing and fixing unmapped supplements
- `test_sbml_parsing.py` - SBML parsing validation with iML1515

Run any example:
```bash
python examples/example_audit_mapping.py
```

## Tips

1. **Always audit first**: Use `get_supplement_mapping()` before creating the matrix
2. **Review fuzzy matches**: Check warnings to ensure correct mappings
3. **Use baseline exchanges**: Include minimal media components for biological accuracy
4. **Custom mappings override SBML**: Manually specify problematic mappings
5. **Case insensitive**: "Glucose", "glucose", "GLUCOSE" all match the same

## Troubleshooting

**Issue:** Supplements not mapping
- Check spelling in your data
- Lower fuzzy_threshold
- Add manual mapping

**Issue:** Wrong exchange reaction matched
- Increase fuzzy_threshold
- Use manual mapping to override

**Issue:** COBRApy import error
- Install: `pip install cobra`

**Issue:** SBML file not found
- Check file path
- Download from BiGG Models: http://bigg.ucsd.edu/models/iML1515

## API Reference

Full API documentation is in the docstrings:

```python
help(create_supplement_exchange_matrix)
help(get_supplement_mapping)
```
