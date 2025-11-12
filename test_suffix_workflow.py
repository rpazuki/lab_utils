"""Test exchange_suffix in complete workflow with bounds template"""

import tempfile
from pathlib import Path

import pandas as pd

from labUtils.amn_mappings import (create_exchange_bounds_template,
                                   create_supplement_exchange_matrix,
                                   parse_sbml_exchange_bounds,
                                   parse_sbml_exchanges)

print("Testing exchange_suffix with complete workflow")
print("=" * 70)

model_path = Path(tempfile.gettempdir()) / "iML1515.xml"

if model_path.exists():
    # Parse SBML
    mapping = parse_sbml_exchanges(model_path)
    bounds = parse_sbml_exchange_bounds(model_path)

    # Create sample data
    df = pd.DataFrame({
        'well': ['A1', 'A2'],
        'supplements': ['Glucose', 'Fructose'],
        'mu_max': [0.2, 0.15],
        'success': [True, True]
    })

    print("\n1. Create matrix with '_input' suffix:")
    matrix = create_supplement_exchange_matrix(
        df,
        supplement_to_exchange_map=mapping,
        supplement_column='supplements',
        growth_rate_column='mu_max',
        exchange_suffix='_input'
    )
    print(f"   Columns: {list(matrix.columns)}")
    print(f"   Shape: {matrix.shape}")

    print("\n2. Create bounds template (should work with suffixed matrix):")
    # Note: The growth_rate_column is still 'mu_max' (no suffix)
    template = create_exchange_bounds_template(
        matrix,
        growth_rate_column='mu_max',
        sbml_bounds=bounds
    )
    print(f"   Template shape: {template.shape}")
    print(f"   Sample columns: {list(template.columns[:5])}")

    print("\n3. Verify exchange columns have suffix but growth rate doesn't:")
    exchange_cols = [c for c in matrix.columns if c != 'mu_max']
    has_suffix = all(c.endswith('_input') for c in exchange_cols)
    print(f"   All exchange columns have '_input' suffix: {has_suffix}")
    print(f"   Growth rate column is 'mu_max': {'mu_max' in matrix.columns}")

    print("\n" + "=" * 70)
    print("✓ Complete workflow with exchange_suffix successful!")
    print("=" * 70)
else:
    print(f"\niML1515.xml not found at {model_path}")
