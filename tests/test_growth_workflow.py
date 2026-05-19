"""Integration test for the current growth workflow."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from labUtils.growth_rates import (
   fit_modified_gompertz_per_series,
   predict_modified_gompertz_per_series,
   transform_to_log_n_n0,
)
from labUtils.media_bot import parse, parse_protocol_metadata, parse_raw_CLARIOstar_export


@pytest.fixture
def test_data_path():
   """fixture for test raw data file"""
   return Path(__file__).parent / "rw_data.csv"


@pytest.fixture
def test_meta_path():
   """fixture for test metadata file"""
   return Path(__file__).parent / "protocol_metadata.csv"


def test_complete_growth_workflow(test_data_path, test_meta_path):
   """Test the complete workflow from raw data to transformed growth-curve fitting."""
   raw_long = parse_raw_CLARIOstar_export(test_data_path, value_column_name="od600")
   meta = parse_protocol_metadata(test_meta_path)
   df = parse(raw_long, meta, value_column_name="od600")

   assert isinstance(df, pd.DataFrame)
   assert not df.empty
   assert set(df["well"]) == {"A1", "B1", "C1", "A2", "B2", "C2"}
   assert "time_h" in df.columns
   assert "od600" in df.columns

   df_transformed = transform_to_log_n_n0(
      df,
      value_col="od600",
      transformed_col="log_n_n0",
      group_cols=["well"],
   )
   params_df = fit_modified_gompertz_per_series(
      df_transformed,
      time_col="time_h",
      value_col="log_n_n0",
      group_cols=["well"],
      min_points=5,
   )
   preds_df = predict_modified_gompertz_per_series(
      df_transformed,
      params_df,
      time_col="time_h",
      value_col="log_n_n0",
      group_cols=["well"],
   )

   assert isinstance(params_df, pd.DataFrame)
   assert isinstance(preds_df, pd.DataFrame)
   assert len(params_df) == 6
   assert set(params_df["well"]) == {"A1", "B1", "C1", "A2", "B2", "C2"}
   assert "od600_fit" in preds_df.columns
   assert "residual" in preds_df.columns

   glucose_wells = df_transformed[df_transformed["supplements"] == "Glucose"]["well"].unique()
   glucose_fits = params_df[params_df["well"].isin(glucose_wells)]

   assert all(glucose_fits["success"])
   assert all(glucose_fits["r2"] > 0.95)
   assert all(glucose_fits["mu_max"] > 0)
   successful_wells = params_df.loc[params_df["success"], "well"]
   successful_preds = preds_df[preds_df["well"].isin(successful_wells)]
   assert not successful_preds["od600_fit"].isna().any()
   assert len(preds_df) == len(df_transformed)

   growing_mask = preds_df["well"].isin(glucose_wells)
   growing_residuals = preds_df.loc[growing_mask, "residual"]
   rmse = np.sqrt(np.mean(growing_residuals**2))
   assert rmse < 0.15
