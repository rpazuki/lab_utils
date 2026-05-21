"""Characterization tests for `labUtils.amn_mappings`.

These pin the *current* behavior so the upcoming refactor cannot drift.
They do NOT assert correctness — they assert reproducibility.

No network: PubChem lookups go through the existing `pubchempy not available`
fallback (MW = None), so any hand-crafted mappings_df sets `mmol_concentration`
directly to avoid network dependence.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import pytest

from labUtils.amn_mappings import (
    MediumSource,
    RecordOrigin,
    SearchResult,
    _parse_sbml_exchange_bounds_fallback,
    _parse_sbml_exchanges_fallback,
    build_AMN_inputs_dataframe,
    build_AMN_levels_dataframe,
    build_AMN_levels_for_different_levels,
    build_mappings,
    build_supplement_flux_dataframe,
    load_default_iml1515_mapping,
    load_minimal_media_exchanges,
)

GLUCOSE_MW = 180.156


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


def test_enum_values_match_string_literals():
    assert MediumSource.UNSTATED.value == "unstated"
    assert MediumSource.SUPPLEMENT.value == "supplement"
    assert MediumSource.MEDIUM.value == "medium"
    assert MediumSource.FIXED.value == "fixed"
    assert RecordOrigin.SMBL.value == "smbl"
    assert RecordOrigin.UPDATED.value == "updated"
    assert {m.value for m in SearchResult} == {
        "name", "iupac", "fuzzy_name", "fuzzy_iupac", "not_found"
    }


# ---------------------------------------------------------------------------
# build_mappings — columns, defaults, precedence
# ---------------------------------------------------------------------------


def _minimal_growth_df(supps: str = "") -> pd.DataFrame:
    return pd.DataFrame({
        "well": ["A1"],
        "supplements": [supps],
        "mu_max": [0.0],
    })


def test_build_mappings_returns_expected_columns_and_order():
    df = build_mappings(_minimal_growth_df(), supplement_to_exchange_map={"glucose": "EX_glc__D_e"})
    assert list(df.columns) == [
        "name",
        "iupac_name",
        "other_names",
        "exchange_reaction",
        "mass_per_litre",
        "mmol_concentration",
        "flux_upper_bound",
        "source",
        "record_origin",
    ]


def test_build_mappings_default_yaml_path_is_loaded():
    """No `custom_mapping_file` arg -> default custom_exchange_mapping.yaml is loaded.

    CURRENT BEHAVIOR: the default yaml is an *overlay* on top of the
    `supplement_to_exchange_map` — entries are only applied when their
    exchange_reaction already exists in the base map. So to verify the yaml
    is loaded we register `glucose` in the base map and confirm yaml-provided
    `mass_per_litre=4.0` lands on the glucose row.
    """
    df = build_mappings(
        _minimal_growth_df(), supplement_to_exchange_map={"glucose": "EX_glc__D_e"}
    )
    glc = df.loc[df["name"] == "glucose"].iloc[0]
    assert glc["mass_per_litre"] == 4.0


def test_build_mappings_custom_mapping_overrides_base():
    base = {"glucose": "EX_glc__D_e"}
    custom = {
        "glucose": {
            "exchange_name": "EX_glc__D_e",
            "source": MediumSource.FIXED.value,
            "flux_upper_bound": 99.0,
            "mass_per_litre": 0.0,
        }
    }
    df = build_mappings(_minimal_growth_df(), supplement_to_exchange_map=base, custom_mapping=custom)
    row = df.loc[df["name"] == "glucose"].iloc[0]
    assert row["source"] == MediumSource.FIXED.value
    assert row["flux_upper_bound"] == 99.0


def test_build_mappings_unique_supplements_marked_SUPPLEMENT_with_UPDATED_origin():
    """Supplements observed in `growth_rates_df` get source=SUPPLEMENT, origin=UPDATED."""
    df = build_mappings(_minimal_growth_df("Glucose"), supplement_to_exchange_map={"glucose": "EX_glc__D_e"})
    glc = df.loc[df["name"] == "glucose"].iloc[0]
    assert glc["source"] == MediumSource.SUPPLEMENT.value
    assert glc["record_origin"] == RecordOrigin.UPDATED.value


def test_build_mappings_unstated_supplement_has_source_UNSTATED():
    df = build_mappings(_minimal_growth_df(""), supplement_to_exchange_map={"glucose": "EX_glc__D_e"})
    glc = df.loc[df["name"] == "glucose"].iloc[0]
    # No `growth_rates_df` supplement column reference -> stays UNSTATED.
    assert glc["source"] == MediumSource.UNSTATED.value


def test_build_mappings_other_names_is_sorted_semicolon_string():
    """Synonyms collected via custom_mapping appear in `other_names` as `;`-joined sorted string."""
    custom = {
        "glucose": {
            "exchange_name": "EX_glc__D_e",
            "iupac_name": "(2R,3S,4R,5R)-2,3,4,5,6-Pentahydroxyhexanal",
            "source": MediumSource.SUPPLEMENT.value,
            "mass_per_litre": 0.0,
        }
    }
    df = build_mappings(
        _minimal_growth_df("Glucose"),
        supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
        custom_mapping=custom,
    )
    other = df.loc[df["name"] == "glucose"].iloc[0]["other_names"]
    assert isinstance(other, str)
    if other:
        assert ";".join(sorted(other.split(";"))) == other


def test_build_mappings_returns_sorted_by_exchange_reaction():
    base = {"glucose": "EX_glc__D_e", "oxygen": "EX_o2_e"}
    df = build_mappings(_minimal_growth_df(), supplement_to_exchange_map=base)
    exchanges = df["exchange_reaction"].tolist()
    assert exchanges == sorted(exchanges)


def test_build_mappings_halt_on_error_raises_on_unknown_exchange():
    """When `halt_on_error=True` and `custom_mapping` declares an exchange not in the base map."""
    custom = {
        "mystery": {
            "exchange_name": "EX_does_not_exist",
            "source": MediumSource.MEDIUM.value,
            "mass_per_litre": 0.0,
        }
    }
    with pytest.raises(ValueError):
        build_mappings(
            _minimal_growth_df(),
            supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
            custom_mapping=custom,
            halt_on_error=True,
        )


def test_build_mappings_halt_on_error_false_does_not_raise_on_unknown_exchange():
    custom = {
        "mystery": {
            "exchange_name": "EX_does_not_exist",
            "source": MediumSource.MEDIUM.value,
            "mass_per_litre": 0.0,
        }
    }
    df = build_mappings(
        _minimal_growth_df(),
        supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
        custom_mapping=custom,
        halt_on_error=False,
    )
    # mystery is skipped, glucose still present
    assert "glucose" in df["name"].values
    assert "mystery" not in df["name"].values


def test_build_mappings_custom_mapping_file_explicit_yaml(tmp_path):
    """Explicit yaml filename (relative) is loaded from `yamls/` directory.

    CURRENT BEHAVIOR: `custom_mapping_file` resolves via `Path(__file__).parent / "yamls" / custom_mapping_file`.
    Yaml entries overlay onto exchanges already present in the base map.
    """
    df = build_mappings(
        _minimal_growth_df(),
        supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
        custom_mapping_file="custom_exchange_mapping.yaml",
    )
    glc = df.loc[df["name"] == "glucose"].iloc[0]
    assert glc["mass_per_litre"] == 4.0


# ---------------------------------------------------------------------------
# build_supplement_flux_dataframe
# ---------------------------------------------------------------------------


def _hand_built_mappings() -> pd.DataFrame:
    """Bypass build_mappings entirely so MW values are deterministic (no PubChem)."""
    mmol_glc = (4.0 / GLUCOSE_MW) * 1000.0
    rows = [
        {
            "name": "glucose", "iupac_name": "", "other_names": "",
            "exchange_reaction": "EX_glc__D_e",
            "mass_per_litre": 4.0, "mmol_concentration": mmol_glc,
            "flux_upper_bound": 0.0,
            "source": MediumSource.SUPPLEMENT.value,
            "record_origin": "updated",
        },
        {
            "name": "oxygen", "iupac_name": "", "other_names": "",
            "exchange_reaction": "EX_o2_e",
            "mass_per_litre": 0.0, "mmol_concentration": 0.0,
            "flux_upper_bound": 10.0,
            "source": MediumSource.FIXED.value,
            "record_origin": "updated",
        },
    ]
    return pd.DataFrame(rows)


def test_build_supplement_flux_columns_sorted_then_mu_max():
    mappings = _hand_built_mappings()
    growth = pd.DataFrame({
        "well": ["A1"],
        "supplements": ["glucose"],
        "mu_max": [0.5],
        "success": [True],
        "max_value": [1.0],
        "max_time": [24.0],
        "total_volume_uL": [200.0],
        "mv_mu_max_value": [0.5],
    })
    flux = build_supplement_flux_dataframe(growth, mappings)
    cols = list(flux.columns)
    assert cols[-1] == "mu_max"
    assert cols[:-1] == sorted(cols[:-1])


def test_build_supplement_flux_value_uses_max_time_and_mv_mu_max_value():
    """Pins formula at amn_mappings.py:563 — flux = mmol / (max_time * od_at_gr)."""
    mappings = _hand_built_mappings()
    growth = pd.DataFrame({
        "well": ["A1"],
        "supplements": ["glucose"],
        "mu_max": [0.5],
        "success": [True],
        "max_value": [1.0],
        "max_time": [24.0],
        "total_volume_uL": [200.0],
        "mv_mu_max_value": [0.5],
    })
    flux = build_supplement_flux_dataframe(growth, mappings)
    expected = ((4.0 / GLUCOSE_MW) * 1000.0) / (24.0 * 0.5)
    assert flux["EX_glc__D_e"].iloc[0] == pytest.approx(expected, rel=1e-9)


def test_build_supplement_flux_skips_unsuccessful_rows():
    mappings = _hand_built_mappings()
    growth = pd.DataFrame({
        "well": ["A1", "A2"],
        "supplements": ["glucose", "glucose"],
        "mu_max": [0.5, 0.3],
        "success": [True, False],
        "max_value": [1.0, 1.0],
        "max_time": [24.0, 24.0],
        "total_volume_uL": [200.0, 200.0],
        "mv_mu_max_value": [0.5, 0.5],
    })
    flux = build_supplement_flux_dataframe(growth, mappings)
    assert len(flux) == 1


def test_build_supplement_flux_exchange_suffix_applied_only_to_exchange_cols():
    mappings = _hand_built_mappings()
    growth = pd.DataFrame({
        "well": ["A1"], "supplements": ["glucose"], "mu_max": [0.5],
        "success": [True], "max_value": [1.0], "max_time": [24.0],
        "total_volume_uL": [200.0], "mv_mu_max_value": [0.5],
    })
    flux = build_supplement_flux_dataframe(growth, mappings, exchange_suffix="_i")
    assert "mu_max" in flux.columns
    assert "mu_max_i" not in flux.columns
    for c in flux.columns:
        if c != "mu_max":
            assert c.endswith("_i")


# ---------------------------------------------------------------------------
# build_AMN_inputs_dataframe
# ---------------------------------------------------------------------------


def test_build_amn_inputs_unsuccessful_rows_kept_or_dropped_per_current_behavior():
    """CURRENT BEHAVIOR: build_AMN_inputs_dataframe also skips rows where success=False."""
    mappings = _hand_built_mappings()
    growth = pd.DataFrame({
        "well": ["A1", "A2"],
        "supplements": ["glucose", "glucose"],
        "mu_max": [0.5, 0.3],
        "success": [True, False],
    })
    out = build_AMN_inputs_dataframe(growth, mappings)
    assert len(out) == 1


def test_build_amn_inputs_dtype_for_exchanges_and_mu_max():
    """CURRENT BEHAVIOR: exchange columns are integer-valued (0/1); mu_max is float."""
    mappings = _hand_built_mappings()
    growth = pd.DataFrame({
        "well": ["A1"], "supplements": ["glucose"], "mu_max": [0.5], "success": [True],
    })
    out = build_AMN_inputs_dataframe(growth, mappings)
    assert pd.api.types.is_float_dtype(out["mu_max"])
    # Exchange columns store ints (0 or 1)
    for c in out.columns:
        if c != "mu_max":
            assert out[c].iloc[0] in (0, 1)


# ---------------------------------------------------------------------------
# build_AMN_levels_dataframe
# ---------------------------------------------------------------------------


def test_build_amn_levels_default_values_for_unknown_exchange():
    """No flux, no custom, no sbml -> default_level=1, default_max_value=1000, ratio=0."""
    mappings = _hand_built_mappings()
    # Build an exchange matrix that includes a NEW exchange not in mappings
    exchange_matrix = pd.DataFrame({
        "EX_unknown_e": [0, 0],
        "mu_max": [0.1, 0.2],
    })
    flux_df = pd.DataFrame({"mu_max": [0.1, 0.2]})
    levels = build_AMN_levels_dataframe(exchange_matrix, mappings, flux_df).set_index("name")
    assert levels.at["level", "EX_unknown_e"] == 1
    assert levels.at["max_value", "EX_unknown_e"] == 1000
    assert levels.at["ratio_drawing", "EX_unknown_e"] == 0


def test_build_amn_levels_sbml_bounds_used_when_no_flux_or_custom():
    mappings = _hand_built_mappings()
    exchange_matrix = pd.DataFrame({"EX_unknown_e": [0], "mu_max": [0.1]})
    flux_df = pd.DataFrame({"mu_max": [0.1]})
    levels = build_AMN_levels_dataframe(
        exchange_matrix, mappings, flux_df,
        sbml_bounds={"EX_unknown_e": (3, 77)},
    ).set_index("name")
    assert levels.at["level", "EX_unknown_e"] == 3
    assert levels.at["max_value", "EX_unknown_e"] == 77


def test_build_amn_levels_fixed_source_uses_level_1():
    mappings = _hand_built_mappings()
    exchange_matrix = pd.DataFrame({"EX_o2_e": [0], "mu_max": [0.1]})
    flux_df = pd.DataFrame({"mu_max": [0.1]})
    levels = build_AMN_levels_dataframe(exchange_matrix, mappings, flux_df).set_index("name")
    assert levels.at["level", "EX_o2_e"] == 1
    assert levels.at["max_value", "EX_o2_e"] == 10  # from flux_upper_bound in mappings


def test_build_amn_levels_ratio_drawing_row_is_zero():
    mappings = _hand_built_mappings()
    exchange_matrix = pd.DataFrame({
        "EX_glc__D_e": [0.0], "EX_o2_e": [10.0], "mu_max": [0.1],
    })
    flux_df = pd.DataFrame({"mu_max": [0.1]})
    levels = build_AMN_levels_dataframe(exchange_matrix, mappings, flux_df).set_index("name")
    assert (levels.loc["ratio_drawing"] == 0).all()


def test_build_amn_levels_first_row_names_are_three_known_labels():
    mappings = _hand_built_mappings()
    exchange_matrix = pd.DataFrame({"EX_glc__D_e": [0.0], "mu_max": [0.1]})
    flux_df = pd.DataFrame({"mu_max": [0.1]})
    levels = build_AMN_levels_dataframe(exchange_matrix, mappings, flux_df)
    assert list(levels["name"]) == ["level", "max_value", "ratio_drawing"]


# ---------------------------------------------------------------------------
# build_AMN_levels_for_different_levels
# ---------------------------------------------------------------------------


def _tiny_amn_levels_template() -> pd.DataFrame:
    return pd.DataFrame({
        "name": ["level", "max_value", "ratio_drawing"],
        "EX_a": [1, 20, 0],     # supplement-like (level==1)
        "EX_b": [100, 1000, 0], # variable (level!=1)
    })


def test_build_AMN_levels_for_different_levels_count():
    out = build_AMN_levels_for_different_levels(
        _tiny_amn_levels_template(),
        fixed_levels=[10, 50, 100],
        fixed_max_values=[5, 20],
    )
    assert len(out) == 6


def test_build_AMN_levels_for_different_levels_none_returns_empty_list():
    out = build_AMN_levels_for_different_levels(_tiny_amn_levels_template(), fixed_levels=None)
    assert out == []


def test_build_AMN_levels_for_different_levels_modifies_only_per_level_branch():
    """level==1 column -> set max_value; else -> set level."""
    out = build_AMN_levels_for_different_levels(
        _tiny_amn_levels_template(),
        fixed_levels=[42],
        fixed_max_values=[7],
    )
    df = out[0]
    # EX_a started with level=1 -> max_value updated to 7, level untouched
    assert df.at[0, "EX_a"] == 1
    assert df.at[1, "EX_a"] == 7
    # EX_b started with level=100 -> level updated to 42, max_value untouched
    assert df.at[0, "EX_b"] == 42
    assert df.at[1, "EX_b"] == 1000


# ---------------------------------------------------------------------------
# SBML XML fallback parsers
# ---------------------------------------------------------------------------


_MINIMAL_SBML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
  <model id="m">
    <listOfReactions>
      <reaction id="EX_glc__D_e" name="D-Glucose exchange"/>
      <reaction id="EX_o2_e" name="Oxygen exchange">
        <kineticLaw>
          <listOfParameters>
            <parameter id="upperFluxBound" value="500"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="BIOMASS_core" name="biomass"/>
      <reaction id="EX_huge_e" name="huge bound">
        <kineticLaw>
          <listOfParameters>
            <parameter id="upperFluxBound" value="1000000000"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
"""


def test_parse_sbml_exchanges_fallback_only_returns_EX_reactions(tmp_path):
    sbml_path = tmp_path / "toy.sbml"
    sbml_path.write_text(_MINIMAL_SBML)
    mapping = _parse_sbml_exchanges_fallback(sbml_path)
    # Only EX_ reactions appear; BIOMASS_core is excluded.
    keys = set(mapping.values())
    assert "EX_glc__D_e" in keys
    assert "EX_o2_e" in keys
    assert "BIOMASS_core" not in keys


def test_parse_sbml_exchange_bounds_fallback_caps_large_value_at_1000(tmp_path):
    sbml_path = tmp_path / "toy.sbml"
    sbml_path.write_text(_MINIMAL_SBML)
    bounds = _parse_sbml_exchange_bounds_fallback(sbml_path, default_level=1)
    # EX_huge_e had value=1e9 -> capped at 1000
    assert bounds["EX_huge_e"] == (1, 1000)


def test_parse_sbml_exchange_bounds_fallback_default_level_propagates(tmp_path):
    sbml_path = tmp_path / "toy.sbml"
    sbml_path.write_text(_MINIMAL_SBML)
    bounds = _parse_sbml_exchange_bounds_fallback(sbml_path, default_level=7)
    for level, _maxv in bounds.values():
        assert level == 7


# ---------------------------------------------------------------------------
# iML1515 defaults
# ---------------------------------------------------------------------------


def test_load_default_iml1515_mapping_includes_expected_entries():
    m = load_default_iml1515_mapping()
    assert isinstance(m, dict)
    # Pin a handful of entries we rely on
    assert m["glucose"] == "EX_glc__D_e"
    assert m["fructose"] == "EX_fru_e_i"
    assert m["adenine"] == "EX_ade_e"
    assert m["uracil"] == "EX_ura_e"
    # Pin total count snapshot (current behavior).
    assert len(m) == 25


def test_load_minimal_media_exchanges_returns_known_list():
    items = load_minimal_media_exchanges()
    # All entries are EX_*_e_i
    assert all(s.startswith("EX_") for s in items)
    # Pin count and a few specific membership facts
    assert len(items) == 23
    must_have = {"EX_o2_e_i", "EX_h2o_e_i", "EX_pi_e_i", "EX_nh4_e_i", "EX_so4_e_i"}
    assert must_have.issubset(set(items))
