"""
Test SBML parsing with iML1515 E. coli model.

This script downloads the iML1515 model from BiGG Models database
and validates the SBML parsing functionality.
"""
# pylint: disable=broad-except  # Test script: continue on errors to see all results

import tempfile
import traceback
import urllib.request
from pathlib import Path

import pandas as pd

from labUtils.amn_mappings import (create_supplement_exchange_matrix,
                                   get_supplement_mapping,
                                   load_minimal_media_exchanges,
                                   parse_sbml_exchanges)


def download_iml1515_model():
    """Download iML1515 model from BiGG Models database."""
    url = "http://bigg.ucsd.edu/static/models/iML1515.xml"
    temp_dir = Path(tempfile.gettempdir())
    output_path = temp_dir / "iML1515.xml"

    if output_path.exists():
        print(f"Model already exists at {output_path}")
        return output_path

    print(f"Downloading iML1515 model from {url}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"Model saved to {output_path}")
        return output_path
    except Exception as e:
        print(f"Failed to download model: {e}")
        return None


def test_sbml_parsing():
    """Test SBML parsing with iML1515 model."""
    print("=" * 70)
    print("Testing SBML Parsing with iML1515 Model")
    print("=" * 70)

    # Download model
    model_path = download_iml1515_model()
    if model_path is None:
        print("Cannot proceed without model file.")
        return

    # Create sample growth data
    print("\n1. Creating sample growth data...")
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3', 'A4', 'A5'],
        'strain': ['MG1655'] * 5,
        'supplements': [
            'Glucose',
            'Glucose; Adenine',
            'Fructose',
            'Maltose; Trehalose',
            'Ribose; Glycerol'
        ],
        'mu_max': [0.169, 0.145, 0.134, 0.189, 0.156],
        'lambda': [2.5, 2.8, 3.0, 2.7, 2.9],
        'r2': [0.99, 0.98, 0.99, 0.99, 0.98],
    })
    print(growth_data)

    # Parse SBML model to get mapping
    print("\n2. Parsing SBML model...")
    try:
        sbml_mapping = parse_sbml_exchanges(model_path)
        print(f"   Found {len(sbml_mapping)} metabolite mappings from SBML")
    except Exception as e:
        print(f"   Error parsing SBML: {e}")
        return

    # Test get_supplement_mapping
    print("\n3. Testing get_supplement_mapping...")
    try:
        mapping_info = get_supplement_mapping(
            growth_data,
            supplement_to_exchange_map=sbml_mapping,
            supplement_column='supplements',
            fuzzy_threshold=0.6
        )

        print("\n   Exact matches:")
        for supp, rxn in mapping_info["exact_matches"].items():
            print(f"     {supp:20s} -> {rxn}")

        print("\n   Fuzzy matches:")
        for supp, info in mapping_info["fuzzy_matches"].items():
            print(f"     {supp:20s} -> {info}")

        print("\n   Unmapped supplements:")
        for supp in mapping_info["unmapped"]:
            print(f"     {supp}")

    except Exception as e:
        print(f"   Error in get_supplement_mapping: {e}")
        traceback.print_exc()

    # Test create_supplement_exchange_matrix
    print("\n4. Testing create_supplement_exchange_matrix...")
    try:
        baseline = load_minimal_media_exchanges()

        matrix = create_supplement_exchange_matrix(
            growth_data,
            supplement_to_exchange_map=sbml_mapping,
            supplement_column='supplements',
            growth_rate_column='mu_max',
            baseline_exchanges=baseline[:5],  # Use just first 5 for cleaner output
            separator=';',
            fuzzy_threshold=0.6
        )

        print(f"\n   Matrix shape: {matrix.shape}")
        print("\n   First few rows:")
        print(matrix.head())

        print(f"\n   Column names ({len(matrix.columns)} total):")
        for i, col in enumerate(matrix.columns):
            print(f"     {i+1:2d}. {col}")
            if i >= 15:  # Show first 15 columns
                print(f"     ... and {len(matrix.columns) - 16} more")
                break

        # Check which supplements got mapped
        print("\n   Supplement columns (non-baseline):")
        baseline_set = set(baseline[:5])
        supp_cols = [col for col in matrix.columns if col not in baseline_set and col != 'mu_max']
        for col in supp_cols:
            print(f"     {col}")

    except Exception as e:
        print(f"   Error in create_supplement_exchange_matrix: {e}")
        traceback.print_exc()

    # Test with manual mapping override
    print("\n5. Testing with manual mapping override...")
    try:
        # Combine SBML mapping with manual overrides
        manual_map = {
            **sbml_mapping,  # Start with SBML
            "glucose": "EX_glc__D_e",  # Override
            "adenine": "EX_ade_e",  # Override
            "custom_supplement": "EX_custom_e"  # Add custom
        }

        matrix2 = create_supplement_exchange_matrix(
            growth_data,
            supplement_to_exchange_map=manual_map,
            supplement_column='supplements',
            growth_rate_column='mu_max',
        )

        print(f"\n   With manual override, matrix shape: {matrix2.shape}")
        print(f"   Columns: {list(matrix2.columns)}")

    except Exception as e:
        print(f"   Error with manual override: {e}")
        traceback.print_exc()

    print("\n" + "=" * 70)
    print("SBML Parsing Test Complete!")
    print("=" * 70)


if __name__ == "__main__":
    test_sbml_parsing()
