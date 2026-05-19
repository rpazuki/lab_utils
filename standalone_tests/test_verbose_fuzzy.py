"""Fuzzy-match tests for the current mapping API."""

import pandas as pd

from labUtils.amn_mappings import build_mappings, load_default_iml1515_mapping


def test_build_mappings_fuzzy_match_updates_aliases():
    growth_df = pd.DataFrame(
        {
            "well": ["A1"],
            "supplements": ["Glucoze; Adenin"],
            "mu_max": [0.2],
        }
    )

    mappings_df = build_mappings(
        growth_df,
        supplement_to_exchange_map=load_default_iml1515_mapping(),
        fuzzy_threshold=0.6,
        verbose=True,
    )

    glucose_row = mappings_df.loc[mappings_df["exchange_reaction"] == "EX_glc__D_e"].iloc[0]
    adenine_row = mappings_df.loc[mappings_df["exchange_reaction"] == "EX_ade_e"].iloc[0]

    assert glucose_row["source"] == "supplement"
    assert "glucoze" in glucose_row["other_names"]
    assert adenine_row["source"] == "supplement"
    assert "adenin" in adenine_row["other_names"]
