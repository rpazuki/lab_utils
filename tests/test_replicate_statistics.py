"""Tests for replicate statistics helpers."""
from pathlib import Path

import pytest

from labUtils.media_bot import calculate_replicate_statistics_by_well, parse


@pytest.fixture
def parsed_df():
    tests_dir = Path(__file__).parent
    return parse(tests_dir / "rw_data.csv", tests_dir / "protocol_metadata.csv")


def test_replicate_statistics_alphabetical(parsed_df):
    stats_alpha = calculate_replicate_statistics_by_well(
        parsed_df,
        direction="alphabetical",
        sample_size=2,
        ddof=1,
    )

    assert set(stats_alpha["group_id"].unique()) == {"A_1-2", "B_1-2", "C_1-2"}
    assert {"od_mean", "od_std", "n_replicates"}.issubset(stats_alpha.columns)
    assert (stats_alpha["n_replicates"] == 2).all()


def test_replicate_statistics_numerical(parsed_df):
    stats_num = calculate_replicate_statistics_by_well(
        parsed_df,
        direction="numerical",
        sample_size=3,
        ddof=1,
    )

    assert set(stats_num["group_id"].unique()) == {"ABC_1", "ABC_2"}
    assert (stats_num["n_replicates"] == 3).all()


def test_replicate_statistics_ddof_affects_std(parsed_df):
    stats_unbiased = calculate_replicate_statistics_by_well(parsed_df, direction="numerical", sample_size=3, ddof=1)
    stats_biased = calculate_replicate_statistics_by_well(parsed_df, direction="numerical", sample_size=3, ddof=0)

    comparison = stats_unbiased[stats_unbiased["time_h"] == 0.0][["group_id", "od_std"]].merge(
        stats_biased[stats_biased["time_h"] == 0.0][["group_id", "od_std"]],
        on="group_id",
        suffixes=("_unbiased", "_biased"),
    )
    assert (comparison["od_std_unbiased"] >= comparison["od_std_biased"]).all()


def test_replicate_statistics_custom_value_column():
    tests_dir = Path(__file__).parent
    df_custom = parse(tests_dir / "rw_data.csv", tests_dir / "protocol_metadata.csv", value_column_name="absorbance")
    stats_custom = calculate_replicate_statistics_by_well(
        df_custom,
        direction="numerical",
        sample_size=3,
        value_column_name="absorbance",
    )

    assert {"absorbance_mean", "absorbance_std"}.issubset(stats_custom.columns)
    assert not stats_custom.empty
