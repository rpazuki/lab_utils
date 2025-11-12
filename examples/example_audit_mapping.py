"""
Example usage of get_supplement_mapping utility.

This demonstrates how to audit supplement mappings before creating
the exchange matrix.
"""

import pandas as pd

from labUtils.metabolic_mapping import (create_supplement_exchange_matrix,
                                        get_supplement_mapping,
                                        load_default_iml1515_mapping,
                                        load_minimal_media_exchanges)


def example_audit_mapping():
    """Example: Audit supplement mappings before creating matrix"""

    # Sample growth data
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3', 'A4'],
        'supplements': [
            'Glucose',
            'Glucose; Adenine; Uracil',
            'Maltose',
            'Glucose; Ribose; Weird_Supplement'  # This one won't map
        ],
        'mu_max': [0.169, 0.134, 0.189, 0.156],
    })

    print("Sample growth data:")
    print(growth_data)
    print("\n" + "="*60)

    # Get mapping information for auditing
    mapping = load_default_iml1515_mapping()

    mapping_info = get_supplement_mapping(
        growth_data,
        supplement_to_exchange_map=mapping,
        supplement_column='supplements',
        fuzzy_threshold=0.6
    )

    print("\n1. EXACT MATCHES:")
    print("-" * 60)
    for supp, rxn in mapping_info["exact_matches"].items():
        print(f"  ✓ {supp:20s} -> {rxn}")

    print("\n2. FUZZY MATCHES:")
    print("-" * 60)
    if mapping_info["fuzzy_matches"]:
        for supp, info in mapping_info["fuzzy_matches"].items():
            print(f"  ~ {supp:20s} -> {info}")
    else:
        print("  (none)")

    print("\n3. UNMAPPED SUPPLEMENTS:")
    print("-" * 60)
    if mapping_info["unmapped"]:
        for supp in mapping_info["unmapped"]:
            print(f"  ✗ {supp}")
        print("\n  ⚠️  These supplements need manual mapping!")
    else:
        print("  (none)")

    print("\n" + "="*60)
    print("\n4. CREATING MATRIX:")
    print("-" * 60)

    # Now create the matrix with the same settings
    baseline = load_minimal_media_exchanges()

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=baseline,
        supplement_column='supplements',
        growth_rate_column='mu_max',
        fuzzy_threshold=0.6
    )

    print(f"\nMatrix created with shape: {matrix.shape}")
    print("\nFirst few rows:")
    print(matrix.head())

    # Show which columns correspond to supplements
    baseline_set = set(baseline)
    supplement_cols = [col for col in matrix.columns
                      if col not in baseline_set and col != 'mu_max']

    print(f"\nSupplement columns ({len(supplement_cols)}):")
    for col in supplement_cols:
        print(f"  • {col}")

    return matrix, mapping_info


def example_fix_unmapped():
    """Example: Fix unmapped supplements by adding manual mapping"""

    growth_data = pd.DataFrame({
        'well': ['A1', 'A2'],
        'supplements': ['Glucose', 'Custom_Carbon_Source'],
        'mu_max': [0.2, 0.15],
    })

    print("\n" + "="*60)
    print("FIXING UNMAPPED SUPPLEMENTS")
    print("="*60)

    # First attempt - will show unmapped
    mapping = load_default_iml1515_mapping()

    print("\nFirst attempt (before fixing):")
    info1 = get_supplement_mapping(
        growth_data,
        supplement_to_exchange_map=mapping,
    )
    print(f"  Unmapped: {info1['unmapped']}")

    # Add manual mapping for the unmapped supplement
    custom_mapping = {
        **mapping,
        "custom_carbon_source": "EX_custom_e"
    }

    print("\nSecond attempt (after adding custom mapping):")
    info2 = get_supplement_mapping(
        growth_data,
        supplement_to_exchange_map=custom_mapping,
    )
    print(f"  Exact matches: {list(info2['exact_matches'].keys())}")
    print(f"  Unmapped: {info2['unmapped']}")

    # Now create matrix with fixed mapping
    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=custom_mapping,
        growth_rate_column='mu_max',
    )

    print(f"\nMatrix columns: {list(matrix.columns)}")
    print(matrix)


if __name__ == "__main__":
    print("\n" + "="*60)
    print("EXAMPLE 1: Audit Mapping Before Creating Matrix")
    print("="*60)
    matrix, info = example_audit_mapping()

    print("\n\n" + "="*60)
    print("EXAMPLE 2: Fix Unmapped Supplements")
    print("="*60)
    example_fix_unmapped()

    print("\n\n" + "="*60)
    print("DONE! Use get_supplement_mapping() to audit your data")
    print("before creating the exchange matrix.")
    print("="*60)
