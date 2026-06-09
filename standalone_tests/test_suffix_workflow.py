"""End-to-end suffix workflow tests against the current API."""

import pandas as pd

from labUtils.amn_mappings import (
    MediumSource,
    build_AMN_inputs_dataframe,
    build_AMN_levels_dataframe,
    build_mappings,
)


def test_suffix_workflow_from_mappings_to_levels():
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
                "flux_upper_bound": 9,
                "mass_per_litre": 0.0,
            }
        },
    )

    inputs_df = build_AMN_inputs_dataframe(
        growth_df,
        mappings_df,
        exchange_suffix="_input",
    )
    flux_cols = [col for col in inputs_df.columns if col != "mu_max"]
    flux_df = pd.DataFrame({col: [1.0, 2.0] for col in flux_cols} | {"mu_max": [0.2, 0.15]})
    levels_df = build_AMN_levels_dataframe(
        inputs_df,
        mappings_df,
        flux_df,
        exchange_suffix="_input",
    )

    assert "mu_max" in inputs_df.columns
    assert all(col.endswith("_input") for col in flux_cols)
    assert set(levels_df.columns) - {"name"} == set(flux_cols)
