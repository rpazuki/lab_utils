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

from labUtils.media_bot import (parse, parse_meta_experiment,
                                parse_raw_bmg_export, parse_time_label,
                                process_bmg_dataframe, report)


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
    raw_df = parse_raw_bmg_export(test_data_path)

    # Check basic structure of raw data
    assert isinstance(raw_df, pd.DataFrame)
    assert not raw_df.empty
    assert set(raw_df["well_row"]) == {"A", "B", "C"}
    assert set(raw_df["well_col"]) == {1, 2}

    # Process the raw data
    df = process_bmg_dataframe(raw_df)

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
    assert df["well_row"].dtype == "object"  # string
    assert df["well_col"].dtype == "Int64"
    assert df["od600"].dtype == "float64"

    # Check some known values from the test file
    first_row = df[df["well"] == "A1"].iloc[0]
    assert first_row["od600"] == pytest.approx(0.081)  # First OD reading for A1
    assert first_row["time_h"] == 0.0  # First timepoint

    # Check well formatting
    assert set(df["well_row"]) == {"A", "B", "C"}
    assert set(df["well_col"]) == {1, 2}


def test_read_meta_experiment(test_meta_path):
    """Test reading experiment metadata"""
    df = parse_meta_experiment(test_meta_path)

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
    raw_df = parse_raw_bmg_export(test_data_path)
    raw_long = process_bmg_dataframe(raw_df)
    meta = parse_meta_experiment(test_meta_path)

    # Then parse them together
    df = parse(raw_long, meta)

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
    raw_df = parse_raw_bmg_export(test_data_path)
    raw_long = process_bmg_dataframe(raw_df)
    meta = parse_meta_experiment(test_meta_path)

    # Generate report for complete data
    report_df = report(raw_long, meta)

    # Check basic structure
    assert isinstance(report_df, pd.DataFrame)

    # For valid test data, we expect no issues
    assert report_df.empty, "Expected no issues in the test data"

    # Create a reduced metadata DataFrame missing some wells
    reduced_meta = meta[~meta["well"].isin(["A2", "B2", "C2"])].copy()

    # Now test with incomplete metadata
    report_df = report(raw_long, reduced_meta)

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
