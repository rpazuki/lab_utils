"""Tests for verbose behavior using the current mapping API."""

import pandas as pd

from labUtils.amn_mappings import build_mappings, load_default_iml1515_mapping


def test_verbose_build_mappings_keeps_unknown_supplements():
    growth_df = pd.DataFrame(
        {
            "well": ["A1", "A2"],
            "supplements": ["Glucose", "UnknownCompound"],
            "mu_max": [0.2, 0.1],
        }
    )

    mappings_df = build_mappings(
        growth_df,
        supplement_to_exchange_map=load_default_iml1515_mapping(),
        verbose=True,
    )

    assert mappings_df.loc[mappings_df["name"] == "unknowncompound"].empty


def test_verbose_build_mappings_sets_known_supplement_source():
    growth_df = pd.DataFrame(
        {
            "well": ["A1"],
            "supplements": ["Glucose"],
            "mu_max": [0.2],
        }
    )

    mappings_df = build_mappings(
        growth_df,
        supplement_to_exchange_map=load_default_iml1515_mapping(),
        verbose=True,
    )

    glucose_row = mappings_df.loc[mappings_df["name"] == "glucose"].iloc[0]
    assert glucose_row["source"] == "supplement"
    assert glucose_row["exchange_reaction"] == "EX_glc__D_e"
