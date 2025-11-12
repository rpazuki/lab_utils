"""
Tests for metabolic mapping functionality.
"""

import pandas as pd
import pytest

from labUtils.metabolic_mapping import (create_supplement_exchange_matrix,
                                        load_default_iml1515_mapping,
                                        load_minimal_media_exchanges)


def test_basic_mapping():
    """Test basic supplement to exchange mapping"""
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2'],
        'supplements': ['Glucose', 'Maltose'],
        'mu_max': [0.2, 0.15],
    })

    mapping = {
        'glucose': 'EX_glc__D_e',
        'maltose': 'EX_malt_e_i',
    }

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=None,
    )

    # Check shape
    assert matrix.shape == (2, 3)  # 2 exchanges + mu_max

    # Check column names
    assert 'EX_glc__D_e' in matrix.columns
    assert 'EX_malt_e_i' in matrix.columns
    assert 'mu_max' in matrix.columns

    # Check values
    assert matrix.loc[0, 'EX_glc__D_e'] == 1
    assert matrix.loc[0, 'EX_malt_e_i'] == 0
    assert matrix.loc[1, 'EX_glc__D_e'] == 0
    assert matrix.loc[1, 'EX_malt_e_i'] == 1

    # Check growth rates
    assert matrix.loc[0, 'mu_max'] == 0.2
    assert matrix.loc[1, 'mu_max'] == 0.15


def test_multiple_supplements():
    """Test handling of multiple semicolon-separated supplements"""
    growth_data = pd.DataFrame({
        'well': ['A1'],
        'supplements': ['Glucose; Adenine; Uracil'],
        'mu_max': [0.134],
    })

    mapping = {
        'glucose': 'EX_glc__D_e',
        'adenine': 'EX_ade_e',
        'uracil': 'EX_ura_e',
    }

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
    )

    # All three should be 1
    assert matrix.loc[0, 'EX_glc__D_e'] == 1
    assert matrix.loc[0, 'EX_ade_e'] == 1
    assert matrix.loc[0, 'EX_ura_e'] == 1


def test_baseline_exchanges():
    """Test that baseline exchanges are always set to 1"""
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2'],
        'supplements': ['Glucose', 'Maltose'],
        'mu_max': [0.2, 0.15],
    })

    mapping = {
        'glucose': 'EX_glc__D_e',
        'maltose': 'EX_malt_e_i',
    }

    baseline = ['EX_o2_e_i', 'EX_nh4_e_i']

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=baseline,
    )

    # Baseline should be 1 for all rows
    assert matrix.loc[0, 'EX_o2_e_i'] == 1
    assert matrix.loc[0, 'EX_nh4_e_i'] == 1
    assert matrix.loc[1, 'EX_o2_e_i'] == 1
    assert matrix.loc[1, 'EX_nh4_e_i'] == 1


def test_case_insensitive_matching():
    """Test that supplement matching is case-insensitive"""
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['GLUCOSE', 'glucose', 'GlUcOsE'],
        'mu_max': [0.2, 0.2, 0.2],
    })

    mapping = {'glucose': 'EX_glc__D_e'}

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
    )

    # All should match
    assert matrix['EX_glc__D_e'].sum() == 3


def test_whitespace_handling():
    """Test that whitespace is properly handled"""
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2'],
        'supplements': ['  Glucose  ', 'Glucose ; Maltose  '],
        'mu_max': [0.2, 0.15],
    })

    mapping = {
        'glucose': 'EX_glc__D_e',
        'maltose': 'EX_malt_e_i',
    }

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
    )

    assert matrix.loc[0, 'EX_glc__D_e'] == 1
    assert matrix.loc[1, 'EX_glc__D_e'] == 1
    assert matrix.loc[1, 'EX_malt_e_i'] == 1


def test_missing_supplement_column():
    """Test handling when supplement column is missing"""
    growth_data = pd.DataFrame({
        'well': ['A1'],
        'mu_max': [0.2],
    })

    mapping = {'glucose': 'EX_glc__D_e'}

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
    )

    # Should only have mu_max column when supplement column is missing
    assert 'mu_max' in matrix.columns
    assert matrix.shape[0] == 1
    assert matrix.loc[0, 'mu_max'] == 0.2


def test_empty_supplements():
    """Test handling of empty or NaN supplement values"""
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3'],
        'supplements': ['Glucose', '', None],
        'mu_max': [0.2, 0.15, 0.18],
    })

    mapping = {'glucose': 'EX_glc__D_e'}

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
    )

    assert matrix.loc[0, 'EX_glc__D_e'] == 1
    assert matrix.loc[1, 'EX_glc__D_e'] == 0
    assert matrix.loc[2, 'EX_glc__D_e'] == 0


def test_column_order():
    """Test that mu_max is always the last column"""
    growth_data = pd.DataFrame({
        'well': ['A1'],
        'supplements': ['Glucose; Maltose; Ribose'],
        'mu_max': [0.2],
    })

    mapping = {
        'glucose': 'EX_glc__D_e',
        'maltose': 'EX_malt_e_i',
        'ribose': 'EX_rib__D_e_i',
    }

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
    )

    # mu_max should be last
    assert matrix.columns[-1] == 'mu_max'

    # Other columns should be alphabetically sorted
    other_cols = list(matrix.columns[:-1])
    assert other_cols == sorted(other_cols)


def test_load_default_mapping():
    """Test that default iML1515 mapping loads correctly"""
    mapping = load_default_iml1515_mapping()

    assert isinstance(mapping, dict)
    assert len(mapping) > 0

    # Check some expected entries
    assert 'glucose' in mapping
    assert 'ribose' in mapping
    assert 'maltose' in mapping

    # Check values are exchange reaction IDs
    assert mapping['glucose'].startswith('EX_')


def test_load_minimal_media():
    """Test that minimal media exchanges load correctly"""
    baseline = load_minimal_media_exchanges()

    assert isinstance(baseline, list)
    assert len(baseline) > 0

    # Check some expected entries
    assert 'EX_o2_e_i' in baseline
    assert 'EX_nh4_e_i' in baseline
    assert 'EX_pi_e_i' in baseline


def test_no_mapping_provided():
    """Test that error is raised when no mapping is provided"""
    growth_data = pd.DataFrame({
        'well': ['A1'],
        'supplements': ['Glucose'],
        'mu_max': [0.2],
    })

    with pytest.raises(TypeError):
        # Missing required parameter supplement_to_exchange_map
        create_supplement_exchange_matrix(growth_data, None) # type: ignore


def test_custom_separator():
    """Test using a custom separator for supplements"""
    growth_data = pd.DataFrame({
        'well': ['A1'],
        'supplements': ['Glucose, Maltose, Ribose'],
        'mu_max': [0.2],
    })

    mapping = {
        'glucose': 'EX_glc__D_e',
        'maltose': 'EX_malt_e_i',
        'ribose': 'EX_rib__D_e_i',
    }

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        separator=',',
    )

    assert matrix.loc[0, 'EX_glc__D_e'] == 1
    assert matrix.loc[0, 'EX_malt_e_i'] == 1
    assert matrix.loc[0, 'EX_rib__D_e_i'] == 1


def test_real_world_scenario():
    """Test a realistic scenario matching the example CSV"""
    growth_data = pd.DataFrame({
        'well': ['A1', 'A2', 'A3', 'A4'],
        'strain': ['MG1655'] * 4,
        'media_type': ['M9'] * 4,
        'supplements': [
            'Glucose',
            'Fructose',
            'Maltose',
            'Melibiose',
        ],
        'mu_max': [0.169, 0.134, 0.189, 0.199],
        'lambda': [2.5, 3.0, 2.8, 3.2],
        'r2': [0.99, 0.98, 0.99, 0.98],
    })

    mapping = load_default_iml1515_mapping()
    baseline = load_minimal_media_exchanges()

    matrix = create_supplement_exchange_matrix(
        growth_data,
        supplement_to_exchange_map=mapping,
        baseline_exchanges=baseline,
    )

    # Check dimensions
    assert matrix.shape[0] == 4  # 4 rows
    assert 'mu_max' in matrix.columns

    # Check baseline columns exist
    for ex in baseline:
        assert ex in matrix.columns
        assert matrix[ex].sum() == 4  # All 1s

    # Check supplement-specific exchanges
    assert matrix.loc[0, 'EX_glc__D_e'] == 1
    assert matrix.loc[1, 'EX_fru_e_i'] == 1
    assert matrix.loc[2, 'EX_malt_e_i'] == 1
    assert matrix.loc[3, 'EX_melib_e_i'] == 1

    # Check growth rates are preserved
    assert matrix['mu_max'].tolist() == [0.169, 0.134, 0.189, 0.199]
