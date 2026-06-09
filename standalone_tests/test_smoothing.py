"""Tests for smoothing support in max growth rate fitting."""

import numpy as np
import pandas as pd

from labUtils.growth_rates import fit_max_growth_rate_per_series, transform_to_log_n_n0


def _growth_df(num_points=50):
    np.random.seed(42)
    time = np.linspace(0, 24, num_points)
    y_true = np.exp(0.3 * time)
    y_noisy = y_true + np.random.normal(0, 0.1, len(time))
    df = pd.DataFrame({"time": time, "OD": y_noisy, "well": ["A1"] * len(time)})
    return transform_to_log_n_n0(df, value_col="OD", transformed_col="log_n_n0", group_cols=["well"])


def test_fit_max_growth_rate_with_smoothing_runs():
    df = _growth_df()
    baseline = fit_max_growth_rate_per_series(
        df,
        time_col="time",
        value_col="log_n_n0",
        group_cols=["well"],
        moving_window_size=5,
        smoothing_iterations=0,
    )
    smoothed = fit_max_growth_rate_per_series(
        df,
        time_col="time",
        value_col="log_n_n0",
        group_cols=["well"],
        moving_window_size=5,
        smoothing_iterations=1,
        smooth_window_size=3,
    )

    assert baseline["success"].iloc[0]
    assert smoothed["success"].iloc[0]
    assert baseline["mv_mu_max"].iloc[0] > 0
    assert smoothed["mv_mu_max"].iloc[0] > 0


def test_fit_max_growth_rate_small_dataset_still_returns_row():
    df_small = _growth_df(num_points=10)
    result = fit_max_growth_rate_per_series(
        df_small,
        time_col="time",
        value_col="log_n_n0",
        group_cols=["well"],
        moving_window_size=5,
        smoothing_iterations=1,
        smooth_window_size=3,
    )

    assert len(result) == 1
    assert "success" in result.columns
