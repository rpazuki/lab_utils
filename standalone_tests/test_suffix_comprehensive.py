"""Comprehensive exchange suffix tests for current AMN API."""

import pandas as pd
import pytest

from labUtils.amn_mappings import (
    MediumSource,
    build_AMN_inputs_dataframe,
    build_AMN_levels_dataframe,
    build_mappings,
)


def _build_dataframes(suffix: str | None):
    growth_df = pd.DataFrame(
        {
            "well": ["A1", "A2"],
            "supplements": ["Glucose", "Fructose"],
            "mu_max": [0.2, 0.15],
            "success": [True, True],
        }
    )
    mappings_df = build_mappings(
        growth_df,
        supplement_to_exchange_map={
            "glucose": "EX_glc__D_e",
            "fructose": "EX_fru_e_i",
            "oxygen": "EX_o2_e_i",
        },
        custom_mapping={
            "oxygen": {
                "exchange_name": "EX_o2_e_i",
                "source": MediumSource.FIXED.value,
                "flux_upper_bound": 12,
                "mass_per_litre": 0.0,
            }
        },
    )
    inputs_df = build_AMN_inputs_dataframe(growth_df, mappings_df, exchange_suffix=suffix)

    flux_columns = [c for c in inputs_df.columns if c != "mu_max"]
    flux_df = pd.DataFrame({c: [1.0, 2.0] for c in flux_columns} | {"mu_max": [0.2, 0.15]})
    levels_df = build_AMN_levels_dataframe(
        inputs_df,
        mappings_df,
        flux_df,
        exchange_suffix=suffix,
    )
    return inputs_df, levels_df


@pytest.mark.parametrize("suffix", [None, "", "_input", "_v2"])
def test_suffix_variants_keep_growth_column_unchanged(suffix):
    inputs_df, _ = _build_dataframes(suffix)
    assert "mu_max" in inputs_df.columns


def test_levels_columns_match_input_exchange_columns_with_suffix():
    inputs_df, levels_df = _build_dataframes("_input")
    exchange_cols = {c for c in inputs_df.columns if c != "mu_max"}
    assert set(levels_df.columns) - {"name"} == exchange_cols
