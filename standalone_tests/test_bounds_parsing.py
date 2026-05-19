"""SBML bounds parsing tests for the current API."""

from pathlib import Path

import pandas as pd

from labUtils.amn_mappings import (
    _parse_sbml_exchange_bounds_fallback,
    build_AMN_levels_dataframe,
    build_mappings,
)


def test_parse_sbml_exchange_bounds_fallback_reads_upper_bounds(tmp_path):
    sbml_path = Path(tmp_path) / "toy_model.xml"
    sbml_path.write_text(
        """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\">
  <model>
    <listOfReactions>
      <reaction id=\"EX_glc__D_e\" name=\"glucose exchange\">
        <kineticLaw>
          <listOfParameters>
            <parameter id=\"UPPER_BOUND\" value=\"25\"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id=\"EX_o2_e_i\" name=\"oxygen exchange\"/>
      <reaction id=\"R_NOT_EX\" name=\"internal reaction\"/>
    </listOfReactions>
  </model>
</sbml>
""",
        encoding="utf-8",
    )

    bounds = _parse_sbml_exchange_bounds_fallback(sbml_path, default_level=2)

    assert bounds["EX_glc__D_e"] == (2, 25)
    assert bounds["EX_o2_e_i"] == (2, 1000)
    assert "R_NOT_EX" not in bounds


def test_build_amn_levels_custom_bounds_override_mapping_values():
    growth_df = pd.DataFrame(
        {
            "well": ["A1"],
            "supplements": ["Glucose"],
            "mu_max": [0.2],
        }
    )
    mappings_df = build_mappings(
        growth_df,
        supplement_to_exchange_map={"glucose": "EX_glc__D_e"},
    )
    exchange_matrix = pd.DataFrame({"EX_glc__D_e": [1], "mu_max": [0.2]})
    flux_df = pd.DataFrame({"EX_glc__D_e": [5.0], "mu_max": [0.2]})

    levels_df = build_AMN_levels_dataframe(
        exchange_matrix,
        mappings_df,
        flux_df,
        custom_bounds={"EX_glc__D_e": (9, 90)},
    ).set_index("name")

    assert levels_df.at["level", "EX_glc__D_e"] == 9
    assert levels_df.at["max_value", "EX_glc__D_e"] == 90
