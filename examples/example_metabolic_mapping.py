"""
Example usage of metabolic mapping functions.

This script demonstrates how to create a supplement-to-exchange reaction matrix
from growth rate data.
"""

from pathlib import Path

import pandas as pd

from labUtils.metabolic_mapping import (create_supplement_exchange_matrix,
                                        load_default_iml1515_mapping,
                                        load_minimal_media_exchanges)


def example_with_manual_mapping():
    """Example using manual mapping (no SBML file required)"""

    # Example growth rates dataframe (output from fit_modified_gompertz_per_series)
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3', 'A4'],
        'strain': ['MG1655', 'MG1655', 'MG1655', 'MG1655'],
        'supplements': [
            'Glucose',
            'Glucose; Adenine; Uracil',
            'Maltose',
            'Glucose; Ribose'
        ],
        'mu_max': [0.169, 0.134, 0.189, 0.199],
        'r2': [0.99, 0.98, 0.99, 0.98],
    })

    # Use the default iML1515 mapping
    mapping = load_default_iml1515_mapping()
    baseline = load_minimal_media_exchanges()

    # Create the matrix
    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=baseline,
        growth_rate_column='mu_max'
    )

    print("Exchange Reaction Matrix:")
    print(matrix)
    print(f"\nMatrix shape: {matrix.shape}")
    print(f"Columns: {list(matrix.columns)}")

    # Optionally save to CSV if needed
    # matrix.to_csv('exchange_matrix.csv', index=False)

    return matrix


def example_with_custom_mapping():
    """Example with custom supplement mapping"""

    # Custom mapping for specific supplements
    custom_mapping = {
        "glucose": "EX_glc__D_e",
        "ribose": "EX_rib__D_e_i",
        "maltose": "EX_malt_e_i",
        "adenine": "EX_ade_e",
        "uracil": "EX_ura_e",
        "fructose": "EX_fru_e_i",
        "galactose": "EX_gal_e_i",
    }

    # Example data
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['Glucose', 'Maltose; Adenine', 'Fructose; Uracil'],
        'mu_max': [0.2, 0.15, 0.18],
    })

    # Only include specific baseline exchanges
    baseline = ["EX_pi_e_i", "EX_co2_e_i", "EX_nh4_e_i", "EX_o2_e_i"]

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=custom_mapping,
        baseline_exchanges=baseline,
    )

    print("\nCustom Mapping Matrix:")
    print(matrix)

    return matrix


def example_from_csv_data():
    """Example loading from CSV file (like the test data)"""

    # This would load your actual growth rates CSV
    # growth_data = pd.read_csv('path/to/growth_rates.csv')

    # For this example, create sample data
    growth_data = pd.DataFrame({
        'well': ['A1', 'B1', 'C1'],
        'strain': ['MG1655', 'MG1655', 'MG1655'],
        'media_type': ['M9', 'M9', 'M9'],
        'supplements': [
            'Glucose',
            'Fructose; Adenine',
            'Maltose; Trehalose'
        ],
        'mu_max': [0.17, 0.13, 0.19],
        'lambda': [2.5, 3.0, 2.8],
        'r2': [0.99, 0.98, 0.99],
        'success': [True, True, True],
    })

    # Use default mappings
    mapping = load_default_iml1515_mapping()
    baseline = load_minimal_media_exchanges()

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=baseline,
    )

    print("\nFrom CSV Data Matrix:")
    print(matrix)

    # Show which supplements are present
    supplement_cols = [col for col in matrix.columns
                      if col != 'GR_AVG' and col not in baseline]
    print(f"\nVariable supplement columns: {supplement_cols}")

    return matrix


def example_minimal_no_baseline():
    """Example without baseline (only variable supplements)"""

    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['Glucose', 'Maltose', 'Glucose; Ribose'],
        'mu_max': [0.2, 0.15, 0.22],
    })

    mapping = {
        "glucose": "EX_glc__D_e",
        "ribose": "EX_rib__D_e_i",
        "maltose": "EX_malt_e_i",
    }

    # No baseline exchanges
    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=None,  # Only variable supplements
    )

    print("\nMinimal Matrix (no baseline):")
    print(matrix)

    return matrix


if __name__ == "__main__":
    print("=" * 60)
    print("Metabolic Mapping Examples")
    print("=" * 60)

    print("\n1. Using default iML1515 mapping with baseline:")
    example_with_manual_mapping()

    print("\n" + "=" * 60)
    print("\n2. Using custom mapping:")
    example_with_custom_mapping()

    print("\n" + "=" * 60)
    print("\n3. Processing CSV-like data:")
    example_from_csv_data()

    print("\n" + "=" * 60)
    print("\n4. Minimal matrix without baseline:")
    example_minimal_no_baseline()
