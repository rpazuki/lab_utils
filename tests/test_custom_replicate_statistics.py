import pandas as pd
import pytest

from labUtils.media_bot import calculate_replicate_statistics_by_custom


def _build_variable_pattern_group_data() -> pd.DataFrame:
    records = []
    strain_wells = {
        "EGMB_61": ["A1", "B1"],
        "EGMB_62": ["A2", "B2"],
        "EGMB_63": ["A3", "B3"],
        "EGMB_64": ["A4", "B4"],
        "EGMB_65": ["B5"],
        "control": ["C1", "C2"],
    }

    for strain, wells in strain_wells.items():
        for well in wells:
            row = well[0]
            col = int(well[1:])
            for time_h in (0, 1):
                records.append(
                    {
                        "well": well,
                        "well_row": row,
                        "well_col": col,
                        "strain": strain,
                        "time_h": time_h,
                        "time_label": f"{time_h} h",
                        "od": float(col + time_h),
                    }
                )

    return pd.DataFrame(records)


def _build_fixed_pattern_group_data() -> pd.DataFrame:
    records = []
    strain_wells = {
        "EGMB_61": ["A1", "A2"],
        "EGMB_62": ["B1", "B2"],
        "control": ["C1", "C2"],
    }

    for strain, wells in strain_wells.items():
        for well in wells:
            row = well[0]
            col = int(well[1:])
            for time_h in (0, 1):
                records.append(
                    {
                        "well": well,
                        "well_row": row,
                        "well_col": col,
                        "strain": strain,
                        "time_h": time_h,
                        "time_label": f"{time_h} h",
                        "od": float(col + time_h),
                    }
                )

    return pd.DataFrame(records)


def test_pattern_rule_without_sample_size_groups_full_directional_sets():
    df = _build_variable_pattern_group_data()

    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=None,
        custom_rules={
            "EGMB": {
                "pattern": r"EGMB_\d+",
                "direction": "alphabetical",
            }
        },
    )

    assert sorted(result["strain"].unique()) == ["EGMB"]
    assert sorted(result["group_id"].unique()) == ["A_1-4", "B_1-5"]
    assert set(result.loc[result["group_id"] == "A_1-4", "n_replicates"]) == {4}
    assert set(result.loc[result["group_id"] == "B_1-5", "n_replicates"]) == {5}

    row_a_time0 = result[(result["group_id"] == "A_1-4") & (result["time_h"] == 0)].iloc[0]
    row_b_time1 = result[(result["group_id"] == "B_1-5") & (result["time_h"] == 1)].iloc[0]

    assert row_a_time0["od_mean"] == pytest.approx(2.5)
    assert row_b_time1["od_mean"] == pytest.approx(4.0)


def test_pattern_rule_with_sample_size_keeps_per_strain_behavior():
    df = _build_fixed_pattern_group_data()

    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=None,
        custom_rules={
            "EGMB_rule": {
                "pattern": r"EGMB_\d+",
                "direction": "alphabetical",
                "sample_size": 2,
            }
        },
    )

    assert sorted(result["strain"].unique()) == ["EGMB_61", "EGMB_62"]
    assert sorted(result["group_id"].unique()) == ["A_1-2", "B_1-2"]
    assert set(result["n_replicates"]) == {2}


@pytest.mark.parametrize("direction", ["Alphabetica", "numerical", "row", "column"])
def test_custom_rule_direction_aliases_are_supported(direction: str):
    df = _build_fixed_pattern_group_data()

    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=None,
        custom_rules={
            "EGMB_rule": {
                "pattern": r"EGMB_\d+",
                "direction": direction,
                "sample_size": 2,
            }
        },
    )

    assert not result.empty


def test_pattern_and_sample_size_explicit_none_are_supported():
    df = _build_variable_pattern_group_data()

    result = calculate_replicate_statistics_by_custom(
        df,
        strain_pattern=None,
        custom_rules={
            "EGMB": {
                "pattern": r"EGMB_\d+",
                "direction": "alphabetical",
                "sample_size": None,
            },
            "control": {
                "pattern": None,
                "direction": "row",
                "sample_size": 2,
            },
        },
    )

    assert not result.empty
    assert "EGMB" in set(result["strain"].unique())
    assert "control" in set(result["strain"].unique())
