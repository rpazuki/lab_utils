"""Complete workflow demonstration with SBML bounds parsing"""

import os
import tempfile

import pandas as pd

from labUtils.metabolic_mapping import (create_exchange_bounds_template,
                                        create_supplement_exchange_matrix,
                                        parse_sbml_exchange_bounds,
                                        parse_sbml_exchanges)

print("="*70)
print("Complete Workflow with SBML Bounds")
print("="*70)

model_path = os.path.join(tempfile.gettempdir(), 'iML1515.xml')

if os.path.exists(model_path):
    print("\nStep 1: Parse SBML for exchange names")
    name_mapping = parse_sbml_exchanges(model_path)
    print(f"  ✓ Found {len(name_mapping)} metabolite names")

    print("\nStep 2: Parse SBML for exchange bounds")
    bounds_mapping = parse_sbml_exchange_bounds(model_path)
    print(f"  ✓ Found {len(bounds_mapping)} exchange bounds")

    print("\nStep 3: Create growth data")
    df = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['Glucose', 'Fructose', 'Glucose; Ribose'],
        'mu_max': [0.2, 0.15, 0.22],
        'success': [True, True, True]
    })
    print(df)

    print("\nStep 4: Create exchange matrix")
    matrix = create_supplement_exchange_matrix(df, name_mapping)
    print(matrix)

    print("\nStep 5: Create bounds template from SBML")
    template_sbml = create_exchange_bounds_template(matrix, sbml_bounds=bounds_mapping)
    print(template_sbml)

    print("\nStep 6: Override specific bounds with custom values")
    print("  (Glucose gets special treatment: level=5, max_value=250)")
    custom = {'EX_glc__D_e': (5, 250)}
    template_custom = create_exchange_bounds_template(
        matrix,
        sbml_bounds=bounds_mapping,
        custom_bounds=custom
    )
    print(template_custom)

    print("\nStep 7: Full precedence demonstration")
    print("  Defaults: level=2, max_value=999")
    print("  SBML: Fructose gets (1, 1000) from file")
    print("  Custom: Glucose overridden to (10, 100)")
    template_full = create_exchange_bounds_template(
        matrix,
        default_level=2,
        default_max_value=999,
        sbml_bounds=bounds_mapping,
        custom_bounds={'EX_glc__D_e': (10, 100)}
    )
    print(template_full)

    print("\n" + "="*70)
    print("✓ Complete workflow successful!")
    print("="*70)

else:
    print(f"\niML1515.xml not found at {model_path}")
    print("Please run examples/test_sbml_parsing.py first to download it.")
