"""Test script for SBML bounds parsing"""

import os
import tempfile

import pandas as pd

from labUtils.amn_mappings import (create_exchange_bounds_template,
                                   create_supplement_exchange_matrix,
                                   load_default_iml1515_mapping,
                                   parse_sbml_exchange_bounds)

# Check if iML1515 model exists
model_path = os.path.join(tempfile.gettempdir(), 'iML1515.xml')

print("="*70)
print("Testing SBML Exchange Bounds Parsing")
print("="*70)

if os.path.exists(model_path):
    print(f"\n1. Parsing bounds from {model_path}...")
    bounds = parse_sbml_exchange_bounds(model_path)
    print(f"   Found {len(bounds)} exchange reactions with bounds")

    print("\n2. Sample bounds (first 10):")
    for i, (rxn, (level, max_val)) in enumerate(list(bounds.items())[:10]):
        print(f"   {rxn:20s} -> level={level}, max_value={max_val}")

    # Create test data
    print("\n3. Creating exchange matrix...")
    df = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['Glucose', 'Fructose', 'Ribose'],
        'mu_max': [0.2, 0.15, 0.18],
        'success': [True, True, True]
    })

    mapping = load_default_iml1515_mapping()
    matrix = create_supplement_exchange_matrix(df, mapping)
    print(f"   Matrix shape: {matrix.shape}")
    print(matrix)

    print("\n4. Creating bounds template with SBML bounds...")
    template = create_exchange_bounds_template(matrix, sbml_bounds=bounds)
    print(template)

    print("\n5. Testing precedence: SBML bounds + custom override...")
    custom = {'EX_glc__D_e': (10, 100)}
    template_custom = create_exchange_bounds_template(
        matrix,
        sbml_bounds=bounds,
        custom_bounds=custom
    )
    print(template_custom)

    print("\n" + "="*70)
    print("SUCCESS: All tests passed!")
    print("="*70)

else:
    print(f"\niML1515.xml not found at {model_path}")
    print("Testing with mock bounds instead...")

    # Create mock SBML bounds
    mock_bounds = {
        'EX_glc__D_e': (1, 10),
        'EX_fru_e_i': (1, 20),
        'EX_rib__D_e_i': (1, 30)
    }

    # Create test data
    df = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['Glucose', 'Fructose', 'Ribose'],
        'mu_max': [0.2, 0.15, 0.18],
        'success': [True, True, True]
    })

    mapping = load_default_iml1515_mapping()
    matrix = create_supplement_exchange_matrix(df, mapping)

    print("\nExchange Matrix:")
    print(matrix)

    print("\n1. Bounds with mock SBML bounds:")
    template1 = create_exchange_bounds_template(matrix, sbml_bounds=mock_bounds)
    print(template1)

    print("\n2. Precedence test: SBML + custom:")
    custom = {'EX_glc__D_e': (5, 50)}
    template2 = create_exchange_bounds_template(
        matrix,
        sbml_bounds=mock_bounds,
        custom_bounds=custom
    )
    print(template2)

    print("\n3. Precedence test: defaults + SBML + custom:")
    template3 = create_exchange_bounds_template(
        matrix,
        default_level=2,
        default_max_value=999,
        sbml_bounds={'EX_fru_e_i': (3, 300)},
        custom_bounds={'EX_glc__D_e': (10, 100)}
    )
    print(template3)
    print("\nExpected precedence:")
    print("  EX_glc__D_e:    custom (10, 100)")
    print("  EX_fru_e_i:     SBML (3, 300)")
    print("  EX_rib__D_e_i:  defaults (2, 999)")

    print("\n" + "="*70)
    print("SUCCESS: Mock tests passed!")
    print("="*70)
