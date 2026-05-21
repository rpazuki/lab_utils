"""Characterization tests for `labUtils.media_bot`.

Pin current behavior before refactor. No network.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from labUtils.media_bot import (
    calculate_replicate_statistics_by_custom,
    calculate_replicate_statistics_by_well,
    parse,
    parse_protocol_metadata,
    parse_raw_CLARIOstar_export,
    parse_time_label,
    report,
)


# ---------------------------------------------------------------------------
# parse_time_label — edge cases
# ---------------------------------------------------------------------------


def test_parse_time_label_unicode_whitespace_handled():
    #   (thin space) and \xa0 (nbsp) are normalized to regular space.
    assert parse_time_label("3 h\xa015 min") == (3.25, 3, 15)


def test_parse_time_label_decimal_hours_not_supported():
    """CURRENT BEHAVIOR: regex matches only integer minutes/hours; '1.5 h' is rejected."""
    assert parse_time_label("1.5 h") == (None, None, None)


def test_parse_time_label_none_returns_zero_tuple():
    assert parse_time_label(None) == (0.0, 0, 0)


def test_parse_time_label_garbage_returns_none_tuple():
    assert parse_time_label("not a time at all") == (None, None, None)


# ---------------------------------------------------------------------------
# parse_raw_CLARIOstar_export — error paths
# ---------------------------------------------------------------------------


def test_parse_raw_clariostar_export_missing_header_raises(tmp_path):
    bad = tmp_path / "no_header.csv"
    bad.write_text("Foo,Bar,Baz\n1,2,3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Could not find the data header row"):
        parse_raw_CLARIOstar_export(bad, value_column_name="od600")


def test_parse_raw_clariostar_export_well_col_dtype_is_Int64(tmp_path):
    raw_text = (
        "Some preamble\n"
        "Well Row,Well Col,Content,Raw Data,Raw Data\n"
        ",,Time,0 h 0 min,0 h 30 min\n"
        "A,1,Sample,0.1,0.11\n"
        "B,2,Sample,0.2,0.21\n"
    )
    p = tmp_path / "rw.csv"
    p.write_text(raw_text, encoding="utf-8")
    df = parse_raw_CLARIOstar_export(p, value_column_name="od600")
    assert df["well_col"].dtype == "Int64"


# ---------------------------------------------------------------------------
# parse_protocol_metadata — branches
# ---------------------------------------------------------------------------


def test_parse_protocol_metadata_missing_well_column_raises(tmp_path):
    meta = """=== Experiment Data ===
Strain,Media Type
SLAB1,2x_M9
"""
    p = tmp_path / "meta_no_well.csv"
    p.write_text(meta, encoding="utf-8")
    with pytest.raises(ValueError, match='missing a "Well" column'):
        parse_protocol_metadata(p)


def test_parse_protocol_metadata_is_blank_inference(tmp_path):
    meta = """=== Experiment Data ===
Well,Strain,Strain Well,Media Type,Supplements
A1,Blank,A1,2x_M9,None
B1,SLAB1,B1,2x_M9,Glucose
C1,BLANK,C1,2x_M9,None
"""
    p = tmp_path / "meta_blank.csv"
    p.write_text(meta, encoding="utf-8")
    df = parse_protocol_metadata(p)
    df = df.set_index("well")
    assert bool(df.at["A1", "is_blank"]) is True
    assert bool(df.at["B1", "is_blank"]) is False
    assert bool(df.at["C1", "is_blank"]) is True  # case-insensitive


def test_parse_protocol_metadata_strain_well_kept_as_separate_column(tmp_path):
    """CURRENT BEHAVIOR: `Strain Well` becomes `strain_well` (kept), `Well` becomes `well`."""
    meta = """=== Experiment Data ===
Well,Strain,Strain Well,Media Type,Supplements
A1,SLAB1,X9,2x_M9,Glucose
"""
    p = tmp_path / "meta.csv"
    p.write_text(meta, encoding="utf-8")
    df = parse_protocol_metadata(p)
    assert "well" in df.columns
    assert "strain_well" in df.columns
    assert df.iloc[0]["well"] == "A1"
    assert df.iloc[0]["strain_well"] == "X9"


def test_parse_protocol_metadata_drops_rows_with_invalid_well_format(tmp_path):
    """CURRENT BEHAVIOR: Well must fullmatch `[A-Za-z]+\\d+`; other rows are filtered."""
    meta = """=== Experiment Data ===
Well,Strain,Strain Well,Media Type,Supplements
A1,SLAB1,A1,2x_M9,Glucose
not_a_well,SLAB2,B1,2x_M9,Sucrose
B2,SLAB3,B2,2x_M9,Acetate
"""
    p = tmp_path / "meta_bad_well.csv"
    p.write_text(meta, encoding="utf-8")
    df = parse_protocol_metadata(p)
    assert set(df["well"]) == {"A1", "B2"}


# ---------------------------------------------------------------------------
# report — issue types
# ---------------------------------------------------------------------------


@pytest.fixture
def real_raw():
    return parse_raw_CLARIOstar_export(Path(__file__).parent / "rw_data.csv", value_column_name="od600")


@pytest.fixture
def real_meta():
    return parse_protocol_metadata(Path(__file__).parent / "protocol_metadata.csv")


def test_report_columns_when_issues_present(real_raw, real_meta):
    """CURRENT BEHAVIOR: columns are at minimum {issue, count}; `wells` is only added for some issues."""
    reduced = real_meta[~real_meta["well"].isin(["A2", "B2", "C2"])].copy()
    rep = report(real_raw, reduced, value_column_name="od600")
    assert "issue" in rep.columns
    assert "count" in rep.columns


def test_report_wells_missing_in_raw(real_raw, real_meta):
    # Drop A1, B1 from raw -> issue should be `wells_missing_in_raw`
    reduced_raw = real_raw[~real_raw["well"].isin(["A1", "B1"])].copy()
    rep = report(reduced_raw, real_meta, value_column_name="od600")
    miss_raw = rep[rep["issue"] == "wells_missing_in_raw"]
    assert not miss_raw.empty
    assert miss_raw["count"].iloc[0] == 2
    assert set(miss_raw["wells"].iloc[0].split(";")) == {"A1", "B1"}


def test_report_duplicate_wells_in_metadata(real_raw, real_meta):
    dup_meta = pd.concat([real_meta, real_meta.iloc[[0]]], ignore_index=True)
    rep = report(real_raw, dup_meta, value_column_name="od600")
    dup_row = rep[rep["issue"] == "duplicate_wells_in_metadata"]
    assert not dup_row.empty
    assert dup_row["count"].iloc[0] == 1


def test_report_non_numeric_values(real_raw, real_meta):
    bad_raw = real_raw.copy()
    bad_raw.loc[bad_raw.index[:5], "od600"] = pd.NA
    rep = report(bad_raw, real_meta, value_column_name="od600")
    non_num = rep[rep["issue"] == "non_numeric_od600_values"]
    assert not non_num.empty
    assert non_num["count"].iloc[0] == 5


# ---------------------------------------------------------------------------
# calculate_replicate_statistics_by_well
# ---------------------------------------------------------------------------


@pytest.fixture
def parsed_real(real_raw, real_meta):
    return parse(real_raw, real_meta, value_column_name="od600")


def test_replicate_stats_group_id_format_alphabetical(parsed_real):
    """`A_1-2` style."""
    import re
    stats = calculate_replicate_statistics_by_well(
        parsed_real, direction="alphabetical", sample_size=2, value_column_name="od600",
    )
    assert all(re.fullmatch(r"[A-Z]_\d+(-\d+)?", g) for g in stats["group_id"].unique())


def test_replicate_stats_group_id_format_numerical(parsed_real):
    """`ABC_1` style."""
    import re
    stats = calculate_replicate_statistics_by_well(
        parsed_real, direction="numerical", sample_size=3, value_column_name="od600",
    )
    assert all(re.fullmatch(r"[A-Z]+_\d+", g) for g in stats["group_id"].unique())


def test_replicate_stats_wells_column_is_comma_joined(parsed_real):
    """CURRENT BEHAVIOR: `wells` joins entries with `,` (not `;`)."""
    stats = calculate_replicate_statistics_by_well(
        parsed_real, direction="alphabetical", sample_size=2, value_column_name="od600",
    )
    sample = stats["wells"].iloc[0]
    assert "," in sample
    # Splitting on comma should yield exactly `sample_size` items
    assert len(sample.split(",")) == 2


def test_replicate_stats_value_column_rename_in_output(parsed_real):
    stats = calculate_replicate_statistics_by_well(
        parsed_real, direction="numerical", sample_size=3, value_column_name="od600",
    )
    assert "od600_mean" in stats.columns
    assert "od600_std" in stats.columns
    # n_replicates emitted
    assert "n_replicates" in stats.columns


def test_replicate_stats_ddof_zero_less_or_equal_one(parsed_real):
    stats_one = calculate_replicate_statistics_by_well(
        parsed_real, direction="numerical", sample_size=3, value_column_name="od600", ddof=1,
    )
    stats_zero = calculate_replicate_statistics_by_well(
        parsed_real, direction="numerical", sample_size=3, value_column_name="od600", ddof=0,
    )
    merged = stats_one[["group_id", "time_h", "od600_std"]].merge(
        stats_zero[["group_id", "time_h", "od600_std"]],
        on=["group_id", "time_h"], suffixes=("_one", "_zero"),
    )
    assert (merged["od600_std_one"] >= merged["od600_std_zero"]).all()


def test_replicate_stats_invalid_direction_raises(parsed_real):
    with pytest.raises(ValueError, match="Invalid direction"):
        calculate_replicate_statistics_by_well(
            parsed_real, direction="diagonal", sample_size=2, value_column_name="od600",
        )


def test_replicate_stats_indivisible_size_raises(parsed_real):
    with pytest.raises(ValueError, match="not divisible"):
        calculate_replicate_statistics_by_well(
            parsed_real, direction="alphabetical", sample_size=5, value_column_name="od600",
        )


def test_replicate_stats_sample_size_zero_raises(parsed_real):
    with pytest.raises(ValueError, match="at least 1"):
        calculate_replicate_statistics_by_well(
            parsed_real, direction="alphabetical", sample_size=0, value_column_name="od600",
        )


def test_replicate_stats_missing_required_column_raises():
    df = pd.DataFrame({"well": ["A1"], "well_row": ["A"]})  # no well_col, no value col
    with pytest.raises(ValueError, match="missing required columns"):
        calculate_replicate_statistics_by_well(df, value_column_name="od")


# ---------------------------------------------------------------------------
# calculate_replicate_statistics_by_custom
# ---------------------------------------------------------------------------


def test_replicate_stats_custom_rule_by_exact_strain(parsed_real):
    """Pass a custom rule keyed by exact strain name."""
    rules = {"SLAB1": {"direction": "numerical", "sample_size": 1}}
    out = calculate_replicate_statistics_by_custom(
        parsed_real,
        strain_pattern=None,
        custom_rules=rules,
        value_column_name="od600",
    )
    assert not out.empty
    # Only SLAB1 group(s) should appear
    assert set(out["strain"].unique()).issubset({"SLAB1"})


def test_replicate_stats_custom_no_rules_returns_empty():
    """CURRENT BEHAVIOR: empty custom_rules -> empty result."""
    parsed = parse(
        Path(__file__).parent / "rw_data.csv",
        Path(__file__).parent / "protocol_metadata.csv",
        value_column_name="od600",
    )
    out = calculate_replicate_statistics_by_custom(
        parsed, custom_rules=None, value_column_name="od600",
    )
    assert isinstance(out, pd.DataFrame)
    assert out.empty


def test_replicate_stats_custom_default_strain_pattern_is_letters_only():
    """CURRENT BEHAVIOR: default strain_pattern=r'[A-Za-z]+' groups SLAB1/SLAB2 as `SLAB`."""
    parsed = parse(
        Path(__file__).parent / "rw_data.csv",
        Path(__file__).parent / "protocol_metadata.csv",
        value_column_name="od600",
    )
    rules = {"SLAB": {"direction": "numerical", "sample_size": 3}}
    out = calculate_replicate_statistics_by_custom(
        parsed, custom_rules=rules, value_column_name="od600",
    )
    assert not out.empty
    assert set(out["strain"].unique()).issubset({"SLAB"})
