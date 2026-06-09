"""Tests for the current Gompertz fitting API."""

import numpy as np
import pandas as pd

from labUtils.growth_rates import (fit_modified_gompertz_per_series, gompertz,
                                   transform_to_log_n_n0)


def _synthetic_df():
    np.random.seed(42)
    time = np.linspace(0, 24, 25)
    y_true = gompertz(time, 0.05, 2.0, 0.4, 3.0, clip_exp=100)
    y_noisy = y_true + np.random.normal(0, 0.01, len(time))
    df = pd.DataFrame({"time": time, "OD": y_noisy, "well": ["A1"] * len(time)})
    return transform_to_log_n_n0(df, value_col="OD", transformed_col="log_n_n0", group_cols=["well"])


def test_fit_modified_gompertz_supports_fixed_params():
    df = _synthetic_df()
    result = fit_modified_gompertz_per_series(df, time_col="time", value_col="log_n_n0", group_cols=["well"])
    assert result["success"].iloc[0]
    assert result["r2"].iloc[0] > 0.95

    fixed_y0 = fit_modified_gompertz_per_series(
        df,
        time_col="time",
        value_col="log_n_n0",
        group_cols=["well"],
        fixed_params={"y0": 0.0},
    )
    assert fixed_y0["success"].iloc[0]
    assert fixed_y0["y0"].iloc[0] == 0.0

    fixed_multiple = fit_modified_gompertz_per_series(
        df,
        time_col="time",
        value_col="log_n_n0",
        group_cols=["well"],
        fixed_params={"y0": 0.0, "lambda": 3.0},
    )
    assert fixed_multiple["success"].iloc[0]
    assert fixed_multiple["lambda"].iloc[0] == 3.0
