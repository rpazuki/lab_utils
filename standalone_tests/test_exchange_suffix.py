"""Tests for exchange_suffix support in current AMN API."""

import pandas as pd

from labUtils.amn_mappings import MediumSource, build_AMN_inputs_dataframe, build_mappings


def _inputs_with_suffix(exchange_suffix: str | None) -> pd.DataFrame:
    growth_df = pd.DataFrame(
        {
            "well": ["A1", "A2"],
            "supplements": ["Glucose", ""],
            "mu_max": [0.2, 0.1],
            "success": [True, True],
        }
    )
    mappings_df = build_mappings(
        growth_df,
        supplement_to_exchange_map={
            "glucose": "EX_glc__D_e",
            "oxygen": "EX_o2_e_i",
        },
        custom_mapping={
            "oxygen": {
                "exchange_name": "EX_o2_e_i",
                "source": MediumSource.FIXED.value,
                "flux_upper_bound": 10,
                "mass_per_litre": 0.0,
            }
        },
    )
    return build_AMN_inputs_dataframe(growth_df, mappings_df, exchange_suffix=exchange_suffix)


def test_exchange_suffix_applies_only_to_exchange_columns():
    inputs = _inputs_with_suffix("_input")

    assert "mu_max" in inputs.columns
    assert "mu_max_input" not in inputs.columns
    assert all(col.endswith("_input") for col in inputs.columns if col != "mu_max")


def test_exchange_suffix_none_keeps_original_exchange_names():
    inputs = _inputs_with_suffix(None)
    assert "EX_glc__D_e" in inputs.columns
    assert "EX_o2_e_i" in inputs.columns
