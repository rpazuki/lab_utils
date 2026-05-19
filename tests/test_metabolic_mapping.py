"""Tests for the current AMN mapping API."""

import pandas as pd

from labUtils.amn_mappings import (
    MediumSource,
    build_AMN_inputs_dataframe,
    build_AMN_levels_dataframe,
    build_mappings,
    load_default_iml1515_mapping,
    load_minimal_media_exchanges,
)


def _mapping_row(mappings_df: pd.DataFrame, name: str) -> pd.Series:
    return mappings_df.loc[mappings_df["name"] == name].iloc[0]


def test_load_default_mapping():
    """Default curated mapping loads and includes common supplements."""
    mapping = load_default_iml1515_mapping()

    assert isinstance(mapping, dict)
    assert len(mapping) > 0

    assert "glucose" in mapping
    assert "ribose" in mapping
    assert "maltose" in mapping
    assert mapping["glucose"].startswith("EX_")


def test_load_minimal_media():
    """Minimal media helper returns expected baseline exchanges."""
    baseline = load_minimal_media_exchanges()

    assert isinstance(baseline, list)
    assert len(baseline) > 0

    assert "EX_o2_e_i" in baseline
    assert "EX_nh4_e_i" in baseline
    assert "EX_pi_e_i" in baseline


def test_build_mappings_sets_source_for_observed_supplements():
    growth_data = pd.DataFrame(
        {
            "well": ["A1", "A2", "A3"],
            "supplements": ["Glucose; Adenine", "Fructose", ""],
            "mu_max": [0.2, 0.1, 0.05],
        }
    )

    mappings_df = build_mappings(
        growth_data,
        supplement_to_exchange_map=load_default_iml1515_mapping(),
    )

    assert _mapping_row(mappings_df, "glucose")["source"] == MediumSource.SUPPLEMENT.value
    assert _mapping_row(mappings_df, "adenine")["source"] == MediumSource.SUPPLEMENT.value
    assert _mapping_row(mappings_df, "fructose")["source"] == MediumSource.SUPPLEMENT.value
    assert _mapping_row(mappings_df, "trehalose")["source"] == MediumSource.UNSTATED.value


def test_build_amn_inputs_dataframe_respects_sources_and_suffix():
    growth_data = pd.DataFrame(
        {
            "well": ["A1", "A2", "A3"],
            "supplements": ["Glucose", "Fructose", ""],
            "mu_max": [0.22, 0.11, 0.05],
            "success": [True, False, True],
        }
    )

    base_mapping = {
        "glucose": "EX_glc__D_e",
        "fructose": "EX_fru_e_i",
        "oxygen": "EX_o2_e_i",
    }
    custom_mapping = {
        "oxygen": {
            "exchange_name": "EX_o2_e_i",
            "source": MediumSource.FIXED.value,
            "flux_upper_bound": 10,
            "mass_per_litre": 0.0,
        },
        "fructose": {
            "exchange_name": "EX_fru_e_i",
            "source": MediumSource.MEDIUM.value,
            "mass_per_litre": 0.0,
        },
    }

    mappings_df = build_mappings(
        growth_data,
        supplement_to_exchange_map=base_mapping,
        custom_mapping=custom_mapping,
    )
    inputs_df = build_AMN_inputs_dataframe(
        growth_data,
        mappings_df,
        exchange_suffix="_input",
    )

    assert len(inputs_df) == 2
    assert "mu_max" in inputs_df.columns
    assert "mu_max_input" not in inputs_df.columns
    assert all(col.endswith("_input") for col in inputs_df.columns if col != "mu_max")

    assert inputs_df["EX_o2_e_i_input"].tolist() == [1, 1]
    assert inputs_df["EX_fru_e_i_input"].tolist() == [0, 0]
    assert inputs_df["EX_glc__D_e_input"].tolist() == [1, 0]


def test_build_amn_levels_dataframe_precedence_with_suffix():
    growth_data = pd.DataFrame(
        {
            "well": ["A1", "A2"],
            "supplements": ["Glucose", ""],
            "mu_max": [0.2, 0.1],
            "success": [True, True],
        }
    )
    base_mapping = {
        "glucose": "EX_glc__D_e",
        "oxygen": "EX_o2_e_i",
    }
    custom_mapping = {
        "oxygen": {
            "exchange_name": "EX_o2_e_i",
            "source": MediumSource.FIXED.value,
            "flux_upper_bound": 30,
            "mass_per_litre": 0.0,
        }
    }

    mappings_df = build_mappings(
        growth_data,
        supplement_to_exchange_map=base_mapping,
        custom_mapping=custom_mapping,
    )
    inputs_df = build_AMN_inputs_dataframe(
        growth_data,
        mappings_df,
        exchange_suffix="_input",
    )
    flux_df = pd.DataFrame(
        {
            "EX_glc__D_e_input": [0.0, 5.0],
            "EX_o2_e_i_input": [7.0, 8.0],
            "mu_max": [0.2, 0.1],
        }
    )

    levels_df = build_AMN_levels_dataframe(
        inputs_df,
        mappings_df,
        flux_df,
        exchange_suffix="_input",
        custom_bounds={"EX_glc__D_e_input": (7, 70)},
    ).set_index("name")

    assert levels_df.at["level", "EX_glc__D_e_input"] == 7
    assert levels_df.at["max_value", "EX_glc__D_e_input"] == 70
    assert levels_df.at["level", "EX_o2_e_i_input"] == 1
    assert levels_df.at["max_value", "EX_o2_e_i_input"] == 8
