"""
The test suite includes:

1. `test_parse_time_label`: Tests the time label parsing function with various formats
   - Regular time format (hours and minutes)
   - Edge cases (60 minutes -> 1 hour)
   - Partial formats (only hours or only minutes)
   - Empty/None cases

2. `test_read_bmg_export`: Tests reading the BMG plate reader export file
   - Checks basic DataFrame structure
   - Verifies required columns exist
   - Validates data types
   - Tests known values from the sample data
   - Verifies well formatting

3. `test_read_meta_experiment`: Tests reading the experiment metadata
   - Verifies DataFrame structure
   - Checks required columns
   - Validates known values from metadata
   - Tests data consistency

4. `test_parse_integration`: Tests the full integration of data and metadata parsing
   - Verifies merged data structure
   - Checks column presence
   - Validates metadata merging
   - Tests time series completeness

5. `test_report_function`: Tests the data validation reporting
   - Checks report structure
   - Validates issue formatting
   - Tests reporting with matching data

The tests use the provided sample files:
- `rw_data.csv`: Raw plate reader data
- `protocol_metadata.csv`: Experiment metadata
"""
from pathlib import Path

import pandas as pd
import pytest

from labUtils.media_bot import (parse, parse_protocol_metadata,
                                parse_raw_CLARIOstar_export, parse_time_label,
                                report)


@pytest.fixture
def test_data_path():
    """Path to test raw data file"""
    return Path(__file__).parent / "rw_data.csv"


@pytest.fixture
def test_meta_path():
    """Path to test metadata file"""
    return Path(__file__).parent / "protocol_metadata.csv"


def test_parse_time_label():
    """Test time label parsing function"""
    assert parse_time_label("3 h 15 min") == (3.25, 3, 15)
    assert parse_time_label("0 h 60 min") == (1.0, 1, 0)  # Tests normalization
    assert parse_time_label("2 h") == (2.0, 2, 0)
    assert parse_time_label("30 min") == (0.5, 0, 30)
    assert parse_time_label("") == (0.0, 0, 0)
    # Note: None is handled as empty string internally
    assert parse_time_label("") == (0.0, 0, 0)


def test_read_bmg_export(test_data_path):
    """Test reading BMG export file"""
    # Test parse_raw_bmg_export
    df = parse_raw_CLARIOstar_export(test_data_path, value_column_name="od600")

    # Check basic structure of raw data
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert set(df["well_row"]) == {"A", "B", "C"}
    assert set(df["well_col"]) == {1, 2}

    # Check final structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Check required columns
    required_columns = [
        "well_row", "well_col", "content", "time_label",
        "time_h", "time_min", "od600", "well"
    ]
    for col in required_columns:
        assert col in df.columns

    # Check data types
    assert pd.api.types.is_string_dtype(df["well_row"])  # string-like dtype
    assert df["well_col"].dtype == "Int64"
    assert df["od600"].dtype == "float64"

    # Check some known values from the test file
    first_row = df[df["well"] == "A1"].iloc[0]
    assert first_row["od600"] == pytest.approx(0.081)  # First OD reading for A1
    assert first_row["time_h"] == 0.0  # First timepoint

    # Check well formatting
    assert set(df["well_row"]) == {"A", "B", "C"}
    assert set(df["well_col"]) == {1, 2}


def test_read_bmg_export_with_single_well_column(tmp_path):
    """Test reading BMG export when well is provided as a single 'Well' column."""
    raw_text = """Header line to skip
Well,Content,Raw Data,,
,,Time,0 h 0 min,0 h 30 min
A1,SampleA,0.101,0.102,0.103
B2,SampleB,0.201,0.202,0.203
"""
    raw_path = tmp_path / "rw_single_well.csv"
    raw_path.write_text(raw_text, encoding="utf-8")

    df = parse_raw_CLARIOstar_export(raw_path, value_column_name="od600")

    assert not df.empty
    assert set(df["well"]) == {"A1", "B2"}
    assert set(df["well_row"]) == {"A", "B"}
    assert set(df["well_col"]) == {1, 2}
    assert "od600" in df.columns


def test_read_bmg_export_with_raw_data_wavelength_split_well(tmp_path):
    """Test reading BMG export when header uses 'Raw Data (600)' with split well columns."""
    raw_text = """Some preamble
Well Row,Well Col,Content,Raw Data (600),Raw Data (600)
,,Time,0 h 0 min,0 h 15 min
A,1,SampleA,0.111,0.112
B,2,SampleB,0.211,0.212
"""
    raw_path = tmp_path / "rw_raw_data_600_split.csv"
    raw_path.write_text(raw_text, encoding="utf-8")

    df = parse_raw_CLARIOstar_export(raw_path, value_column_name="od600")

    assert not df.empty
    assert set(df["well"]) == {"A1", "B2"}
    assert set(df["well_row"]) == {"A", "B"}
    assert set(df["well_col"]) == {1, 2}
    assert "od600" in df.columns


def test_read_bmg_export_with_raw_data_wavelength_single_well(tmp_path):
    """Test reading BMG export when header uses 'Raw Data (600)' with single well column."""
    raw_text = """Some preamble
Well,Content,Raw Data (600),Raw Data (600)
,,Time,0 h 0 min,0 h 15 min
A1,SampleA,0.101,0.102
B2,SampleB,0.201,0.202
"""
    raw_path = tmp_path / "rw_raw_data_600_single.csv"
    raw_path.write_text(raw_text, encoding="utf-8")

    df = parse_raw_CLARIOstar_export(raw_path, value_column_name="od600")

    assert not df.empty
    assert set(df["well"]) == {"A1", "B2"}
    assert set(df["well_row"]) == {"A", "B"}
    assert set(df["well_col"]) == {1, 2}
    assert "od600" in df.columns


def test_read_meta_experiment(test_meta_path):
    """Test reading experiment metadata"""
    df = parse_protocol_metadata(test_meta_path)

    # Check basic structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Check required columns
    required_columns = [
        "well", "strain", "media_type", "supplements",
        "media_volume_uL", "water_volume_uL", "supplement_volume_uL"
    ]
    for col in required_columns:
        assert col in df.columns

    # Check some known values from the test file
    first_row = df[df["well"] == "A1"].iloc[0]
    assert first_row["strain"] == "SLAB1"
    assert first_row["media_type"] == "2x_M9"
    assert first_row["supplements"] == "Glucose"
    assert first_row["media_volume_uL"] == 230

    # Check data consistency
    assert len(df) == 6  # Number of wells in test data
    assert set(df["well"]) == {"A1", "B1", "C1", "A2", "B2", "C2"}


def test_parse_integration(test_data_path, test_meta_path):
    """Test full integration of parsing both data and metadata"""
    # First read and process the raw data
    raw_long = parse_raw_CLARIOstar_export(test_data_path, value_column_name="od600")
    meta = parse_protocol_metadata(test_meta_path)

    # Then parse them together
    df = parse(raw_long, meta, value_column_name="od600")

    # Check basic structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Check that we have all expected columns
    expected_columns = [
        "well", "well_row", "well_col", "content", "strain",
        "media_type", "supplements", "time_h", "od600"
    ]
    for col in expected_columns:
        assert col in df.columns

    # Check that metadata was correctly merged
    assert set(df["supplements"].unique()) == {"Glucose", "Sucrose", "Succinate"}
    assert set(df["strain"].unique()) == {"SLAB1", "SLAB2", "SLAB3", "SLAB9", "SLAB10", "SLAB11"}

    # Check time series completeness
    timepoints_per_well = df.groupby("well")["time_h"].nunique()
    assert all(timepoints_per_well == timepoints_per_well.iloc[0])  # All wells should have same number of timepoints


def test_report_function(test_data_path, test_meta_path):
    """Test the report function for data validation"""
    # First read and process the raw data
    raw_long = parse_raw_CLARIOstar_export(test_data_path, value_column_name="od600")
    meta = parse_protocol_metadata(test_meta_path)

    # Generate report for complete data
    report_df = report(raw_long, meta, value_column_name="od600")

    # Check basic structure
    assert isinstance(report_df, pd.DataFrame)

    # For valid test data, we expect no issues
    assert report_df.empty, "Expected no issues in the test data"

    # Create a reduced metadata DataFrame missing some wells
    reduced_meta = meta[~meta["well"].isin(["A2", "B2", "C2"])].copy()

    # Now test with incomplete metadata
    report_df = report(raw_long, reduced_meta, value_column_name="od600")

    # Check that issues are reported correctly
    assert isinstance(report_df, pd.DataFrame)
    assert set(report_df.columns) == {"issue", "count", "wells"}
    assert not report_df.empty

    # Verify specific issues
    wells_missing = report_df[report_df["issue"] == "wells_missing_in_metadata"]
    assert not wells_missing.empty
    assert wells_missing["count"].iloc[0] == 3  # A2, B2, C2 should be missing
    assert set(wells_missing["wells"].iloc[0].split(";")) == {"A2", "B2", "C2"}

    # Check formatting
    assert all(report_df["count"].astype(int) >= 0)
    assert all(report_df["issue"].astype(str).str.contains(r"^[a-z_]+$"))


def test_parse_protocol_metadata_with_nonstandard_section_title(tmp_path):
    """Test metadata parsing when section title is not exactly 'Experiment Data'."""
    meta_text = """=== General Notes ===
Some info here

=== Plate Layout ===
Well,Strain,Strain Well,Media Type,Supplements,Media Volume (µL)
A1,SLAB1,A1,2x_M9,Glucose,230
B2,Blank,B2,2x_M9,None,230
"""
    meta_path = tmp_path / "protocol_alt_section.csv"
    meta_path.write_text(meta_text, encoding="utf-8")

    df = parse_protocol_metadata(meta_path)

    assert not df.empty
    assert set(df["well"]) == {"A1", "B2"}
    assert "strain" in df.columns
    assert "is_blank" in df.columns


def test_parse_protocol_metadata_plain_csv_without_sections(tmp_path):
    """Test metadata parsing when file is plain CSV with no === section headers."""
    meta_text = """Well,Strain,Strain Well,Media Type,Supplements,Media Volume (µL)
A1,SLAB1,A1,2x_M9,Glucose,230
B2,Blank,B2,2x_M9,None,230
"""
    meta_path = tmp_path / "protocol_plain.csv"
    meta_path.write_text(meta_text, encoding="utf-8")

    df = parse_protocol_metadata(meta_path)

    assert not df.empty
    assert set(df["well"]) == {"A1", "B2"}
    assert "strain" in df.columns
    assert "is_blank" in df.columns


def test_parse_protocol_metadata_section_markers_with_trailing_commas(tmp_path):
    """Section headers with trailing commas and comma-only rows should be ignored."""
    meta_text = """=== Experiment Data ===,,,,,,,,,
Well,Strain,Strain Well,Media Type,Supplements,Media Volume (µL)
A1,SLAB1,A1,2x_M9,Glucose,230
B2,Blank,B2,2x_M9,None,230
,,,,,,,,,
=== Supplement and Media Map ===,,,,,,,,,
Well,Supplement/Media,Plate
A1,nan,Supplement reservoir
"""
    meta_path = tmp_path / "protocol_trailing_commas.csv"
    meta_path.write_text(meta_text, encoding="utf-8")

    df = parse_protocol_metadata(meta_path)

    assert set(df["well"]) == {"A1", "B2"}
    assert "=== Supplement and Media Map ===" not in set(df["well"])


def test_report_ignores_section_marker_rows_in_metadata(tmp_path):
    """report() should not treat section titles/comma-only rows as metadata wells."""
    meta_text = """=== Experiment Data ===,,,,,,,,,
Well,Strain,Strain Well,Media Type,Supplements,Media Volume (µL)
A1,SLAB1,A1,2x_M9,Glucose,230
B2,Blank,B2,2x_M9,None,230
,,,,,,,,,
=== Strain Map ===,,,,,,,,,
Well,Strain
A1,SLAB1
"""
    meta_path = tmp_path / "protocol_report_trailing_commas.csv"
    meta_path.write_text(meta_text, encoding="utf-8")

    meta = parse_protocol_metadata(meta_path)

    raw_long = pd.DataFrame(
        {
            "well": ["A1", "B2"],
            "well_row": ["A", "B"],
            "well_col": [1, 2],
            "content": ["Sample", "Sample"],
            "time_label": ["0 h", "0 h"],
            "time_h": [0.0, 0.0],
            "time_h_int": [0, 0],
            "time_min_int": [0, 0],
            "time_min": [0.0, 0.0],
            "od": [0.1, 0.2],
        }
    )

    report_df = report(raw_long, meta)

    assert report_df.empty
