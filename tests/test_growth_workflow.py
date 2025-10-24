"""
The test does the following:

1. **Data Parsing**:
   - Uses `media_bot.parse()` to combine raw data with metadata
   - Verifies that all wells are present
   - Checks for required columns

2. **Growth Curve Fitting**:
   - Applies the Gompertz model to each well
   - Verifies that all wells were processed
   - Checks for required output columns

3. **Quality Control**:
   - Focuses on glucose conditions (known growing cultures)
   - Verifies high R² values for growing cultures
   - Checks that fitted parameters are biologically meaningful:
     - Initial OD (y0) around 0.1
     - Reasonable growth amplitude (A)
     - Positive but reasonable maximum growth rate
     - Sensible lag times
   - Verifies prediction quality using RMSE

This integration test demonstrates the complete workflow from raw data to fitted growth curves and verifies that:

1. Raw data parsing works correctly
   - Combines plate reader data with metadata
   - Preserves all wells and time points
   - Creates proper time and OD columns

2. Growth curve fitting works as expected
   - Successfully fits all wells
   - Produces high-quality fits for growing conditions (glucose)
   - Generates parameters in biologically reasonable ranges:
     - Initial OD (y0): 0.05 - 0.2
     - Growth amplitude (A): 0.7 - 2.0 OD units
     - Maximum growth rate (μmax): 0.1 - 1.0 h⁻¹
     - Lag time (λ): 1.0 - 12.0 h

3. Output data structure is correct
   - Includes fitted curves and residuals
   - Preserves all input data
   - Good fit quality (RMSE < 0.1)

This test ensures that the entire data processing pipeline works correctly and produces biologically meaningful results.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from labUtils.growth_rates import fit_modified_gompertz_per_series
from labUtils.media_bot import (parse, parse_protocol_metadata,
                                parse_raw_CLARIOstar_export)


@pytest.fixture
def test_data_path():
    """fixture for test raw data file"""
    return Path(__file__).parent / "rw_data.csv"


@pytest.fixture
def test_meta_path():
    """fixture for test metadata file"""
    return Path(__file__).parent / "protocol_metadata.csv"


def test_complete_growth_workflow(test_data_path, test_meta_path):
    """Test the complete workflow from raw data to growth curve fitting"""
    # Step 1: Parse raw data and metadata
    raw_long = parse_raw_CLARIOstar_export(test_data_path)
    meta = parse_protocol_metadata(test_meta_path)
    df = parse(raw_long, meta)

    # Verify basic data parsing
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert set(df["well"]) == {"A1", "B1", "C1", "A2", "B2", "C2"}
    assert "time_h" in df.columns
    assert "od600" in df.columns

    # Step 2: Fit growth curves for each well
    params_df, preds_df = fit_modified_gompertz_per_series(
        df,
        time_col="time_h",
        value_col="od600",
        group_cols=["well"],
        min_points=5
    )

    # Verify fitting results
    assert isinstance(params_df, pd.DataFrame)
    assert isinstance(preds_df, pd.DataFrame)

    # Check that we have results for all wells
    assert len(params_df) == 6  # One row per well
    assert set(params_df["well"]) == {"A1", "B1", "C1", "A2", "B2", "C2"}

    # Verify prediction structure
    assert "od600_fit" in preds_df.columns
    assert "residual" in preds_df.columns

    # Check fitting quality for growing cultures (glucose conditions)
    glucose_wells = df[df["supplements"] == "Glucose"]["well"].unique()
    glucose_fits = params_df[params_df["well"].isin(glucose_wells)]

    # Verify that glucose conditions show growth (good fits)
    assert all(glucose_fits["success"])  # All fits should succeed
    assert all(glucose_fits["r2"] > 0.95)  # Good R² for growing cultures

    # Check biological meaning of parameters for glucose conditions
    for _, row in glucose_fits.iterrows():
        # Initial OD should be close to observed
        assert 0.05 <= row["y0"] <= 0.2

        # Growth amplitude should be reasonable for OD600
        assert 0.7 <= row["A"] <= 2.0

        # Maximum growth rate should be positive but not unreasonable
        assert 0.1 <= row["mu_max"] <= 1.0

        # Lag time should be non-negative and reasonable
        assert 1.0 <= row["lambda"] <= 12.0    # Verify predictions
    assert not preds_df["od600_fit"].isna().any()  # No missing predictions
    assert len(preds_df) == len(df)  # Same number of rows as input

    # Calculate overall fit quality
    growing_mask = preds_df["well"].isin(glucose_wells)
    growing_residuals = preds_df.loc[growing_mask, "residual"]
    rmse = np.sqrt(np.mean(growing_residuals**2))
    assert rmse < 0.1  # RMSE should be small for growing cultures
