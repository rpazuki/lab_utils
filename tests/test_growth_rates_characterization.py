"""Characterization tests for `labUtils.growth_rates`.

Pin current behavior before refactor. No display / network.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # noqa: E402 — must be set before any matplotlib import

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from labUtils.growth_rates import (
    fit_max_growth_rate_per_series,
    fit_modified_gompertz_per_series,
    gompertz,
    plot_and_save,
    plot_single_series,
    predict_modified_gompertz_per_series,
    transform_to_log_n_n0,
)


# ---------------------------------------------------------------------------
# gompertz
# ---------------------------------------------------------------------------


def test_gompertz_a_zero_returns_y0_filled_array():
    t = np.array([0.0, 1.0, 5.0])
    result = gompertz(t, y0=0.7, A_0=0.0, mu_max=0.5, lam=2.0, clip_exp=50.0)
    assert np.allclose(result, 0.7)


def test_gompertz_clip_exp_keeps_values_finite():
    """Very large mu_max would overflow exp(exp(...)); clip_exp keeps finite."""
    t = np.linspace(0, 10, 50)
    result = gompertz(t, y0=0.0, A_0=1.0, mu_max=1e6, lam=2.0, clip_exp=50.0)
    assert np.all(np.isfinite(result))


def test_gompertz_returns_1d_ndarray_for_1d_input():
    t = np.linspace(0, 5, 10)
    result = gompertz(t, 0.1, 1.0, 0.5, 2.0, 50.0)
    assert isinstance(result, np.ndarray)
    assert result.shape == t.shape


# ---------------------------------------------------------------------------
# transform_to_log_n_n0
# ---------------------------------------------------------------------------


def test_transform_to_log_n_n0_adds_n0_and_log_columns():
    df = pd.DataFrame({
        "well": ["A1"] * 5,
        "time_h": [0.0, 1.0, 2.0, 3.0, 4.0],
        "od600": [0.1, 0.15, 0.2, 0.3, 0.4],
    })
    out = transform_to_log_n_n0(df, value_col="od600", transformed_col="log_n_n0", group_cols=["well"])
    assert "log_n_n0" in out.columns
    assert "n0" in out.columns
    # First-point baseline
    assert out["n0"].iloc[0] == 0.1
    # First log_n_n0 is log(0.1/0.1) = 0
    assert out["log_n_n0"].iloc[0] == pytest.approx(0.0)


def test_transform_to_log_n_n0_window_averaging():
    """`OD_0_averaging_window=3` -> n0 = mean of first 3 values."""
    df = pd.DataFrame({
        "well": ["A1"] * 5,
        "time_h": [0.0, 1.0, 2.0, 3.0, 4.0],
        "od600": [0.10, 0.20, 0.30, 0.40, 0.50],
    })
    out = transform_to_log_n_n0(
        df, value_col="od600", group_cols=["well"], OD_0_averaging_window=3,
    )
    assert out["n0"].iloc[0] == pytest.approx(0.20)  # mean(0.1, 0.2, 0.3)


def test_transform_to_log_n_n0_per_group_baseline():
    df = pd.DataFrame({
        "well": ["A1", "A1", "A2", "A2"],
        "time_h": [0.0, 1.0, 0.0, 1.0],
        "od600": [0.10, 0.20, 0.50, 0.60],
    })
    out = transform_to_log_n_n0(df, value_col="od600", group_cols=["well"])
    a1_n0 = out.loc[out["well"] == "A1", "n0"].unique()
    a2_n0 = out.loc[out["well"] == "A2", "n0"].unique()
    assert a1_n0.tolist() == [0.10]
    assert a2_n0.tolist() == [0.50]


def test_transform_to_log_n_n0_zero_baseline_produces_nan():
    """CURRENT BEHAVIOR: n0=0 -> log_n_n0 set to NaN (np.where with n0>0 branch)."""
    df = pd.DataFrame({
        "well": ["A1"] * 3,
        "time_h": [0.0, 1.0, 2.0],
        "od600": [0.0, 0.1, 0.2],
    })
    out = transform_to_log_n_n0(df, value_col="od600", group_cols=["well"])
    assert out["log_n_n0"].isna().all()


# ---------------------------------------------------------------------------
# fit_max_growth_rate_per_series
# ---------------------------------------------------------------------------


def _linear_log_n_n0_df(slope: float = 0.4, n_points: int = 20, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 10, n_points)
    log_vals = slope * t + rng.normal(0, 0.01, size=n_points)
    return pd.DataFrame({
        "well": ["A1"] * n_points,
        "time_h": t,
        "log_n_n0": log_vals,
        "n0": [0.1] * n_points,
    })


def test_fit_max_growth_rate_returns_expected_columns():
    df = _linear_log_n_n0_df()
    params = fit_max_growth_rate_per_series(df, value_col="log_n_n0")
    # CURRENT BEHAVIOR: prefix is `mv_` (moving-window) for r2/rmse but NOT for mu_max.
    expected_cols = {
        "well", "n", "success", "message",
        "mv_mu_max", "mv_r2", "mv_rmse",
        "mv_mu_max_time", "mv_mu_max_value",
        "max_value", "max_time",
    }
    assert expected_cols.issubset(set(params.columns))


def test_fit_max_growth_rate_recovers_known_slope():
    df = _linear_log_n_n0_df(slope=0.4)
    params = fit_max_growth_rate_per_series(df, value_col="log_n_n0", moving_window_size=5)
    assert params["success"].iloc[0]
    assert params["mv_mu_max"].iloc[0] == pytest.approx(0.4, abs=0.05)


def test_fit_max_growth_rate_marks_flat_as_unsuccessful():
    """`mv_mu_max` should report 0/NaN-ish behavior with truly flat data."""
    df = pd.DataFrame({
        "well": ["A1"] * 10,
        "time_h": np.linspace(0, 10, 10),
        "log_n_n0": np.zeros(10),
        "n0": [0.1] * 10,
    })
    params = fit_max_growth_rate_per_series(df, value_col="log_n_n0")
    # CURRENT BEHAVIOR: insufficient/flat path returns success=False
    assert not bool(params["success"].iloc[0])


def test_fit_max_growth_rate_insufficient_data_returns_success_false():
    df = pd.DataFrame({
        "well": ["A1"] * 2,
        "time_h": [0.0, 1.0],
        "log_n_n0": [0.0, 0.1],
        "n0": [0.1, 0.1],
    })
    params = fit_max_growth_rate_per_series(df, value_col="log_n_n0", min_points=5)
    assert not bool(params["success"].iloc[0])


# ---------------------------------------------------------------------------
# fit_modified_gompertz_per_series
# ---------------------------------------------------------------------------


def _gompertz_dataset(seed: int = 42, n: int = 30) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 10, n)
    y = gompertz(t, 0.1, 1.0, 0.5, 2.0, 50.0) + rng.normal(0, 0.02, size=n)
    df = pd.DataFrame({
        "well": ["A1"] * n,
        "strain": ["SLAB"] * n,  # plot_single_series reads `strain`
        "supplements": ["Glucose"] * n,
        "time_h": t,
        "od600": y,
    })
    return transform_to_log_n_n0(df, value_col="od600", group_cols=["well"])


def test_fit_gompertz_returns_full_column_set():
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    expected_cols = {
        "well", "n", "success", "message",
        "y0", "A", "mu_max", "lambda",
        "r2", "rmse",
        "mu_max_time", "mu_max_value",
        "max_value", "max_time",
    }
    assert expected_cols.issubset(set(params.columns))


def test_fit_gompertz_fixed_params_freezes_value():
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(
        df, value_col="log_n_n0", fixed_params={"lambda": 2.0},
    )
    assert params["lambda"].iloc[0] == 2.0


def test_fit_gompertz_message_contains_ok_on_success():
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    assert params["success"].iloc[0]
    assert params["message"].iloc[0] == "ok"


def test_fit_gompertz_runs_deterministic_across_repeated_calls():
    """Same input -> same output, no random state leakage."""
    df = _gompertz_dataset()
    p1 = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    p2 = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    pd.testing.assert_frame_equal(p1, p2)


# ---------------------------------------------------------------------------
# predict_modified_gompertz_per_series
# ---------------------------------------------------------------------------


def test_predict_residual_is_value_minus_fit():
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    preds = predict_modified_gompertz_per_series(df, params, value_col="log_n_n0")
    # residual = log_n_n0 - od600_fit
    assert np.allclose(
        preds["residual"].values,
        preds["log_n_n0"].values - preds["od600_fit"].values,
        equal_nan=True,
    )


def test_predict_input_df_not_modified():
    df = _gompertz_dataset()
    df_copy = df.copy()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    _ = predict_modified_gompertz_per_series(df, params, value_col="log_n_n0")
    pd.testing.assert_frame_equal(df, df_copy)


def test_predict_save_plot_data_false_writes_nothing(tmp_path):
    """CURRENT BEHAVIOR: default `save_plot_data=False` -> no files written even when output_dir provided."""
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    out = tmp_path / "plot_out"
    out.mkdir()
    _ = predict_modified_gompertz_per_series(
        df, params, value_col="log_n_n0", output_dir=out, save_plot_data=False,
    )
    # No PNGs created
    assert list(out.rglob("*.png")) == []


def test_predict_save_plot_data_true_writes_png(tmp_path):
    """CURRENT BEHAVIOR: `save_plot_data=True` with `output_dir` writes PNGs."""
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    out = tmp_path / "plot_out"
    out.mkdir()
    _ = predict_modified_gompertz_per_series(
        df, params, value_col="log_n_n0", output_dir=out, save_plot_data=True,
    )
    pngs = list(out.rglob("*.png"))
    assert len(pngs) >= 1
    for p in pngs:
        assert p.stat().st_size > 0


# ---------------------------------------------------------------------------
# plot_single_series / plot_and_save
# ---------------------------------------------------------------------------


def test_plot_single_series_writes_one_png(tmp_path):
    df = _gompertz_dataset()
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    preds = predict_modified_gompertz_per_series(df, params, value_col="log_n_n0")
    # plot_single_series expects the fit params to be merged into the series df.
    merged = preds.merge(
        params[["well", "y0", "A", "mu_max", "lambda", "r2", "rmse", "success", "message"]],
        on="well", how="left",
    )
    series = merged[merged["well"] == "A1"]
    plot_single_series(
        series, ("A1",), time_col="time_h", value_col="log_n_n0",
        group_cols=["well"], output_dir=tmp_path, save_plot=True,
    )
    plots_dir = tmp_path / "plots"
    pngs = list(plots_dir.rglob("*.png"))
    assert len(pngs) == 1
    assert pngs[0].stat().st_size > 0


def test_plot_and_save_iterates_all_groups(tmp_path):
    """Two wells -> two PNGs."""
    rng = np.random.default_rng(0)
    t = np.linspace(0, 10, 30)
    rows = []
    for well, params_t in [("A1", (0.1, 1.0, 0.5, 2.0)), ("A2", (0.2, 0.8, 0.6, 1.5))]:
        y = gompertz(t, *params_t, 50.0) + rng.normal(0, 0.02, size=len(t))
        for ti, yi in zip(t, y, strict=True):
            rows.append({"well": well, "strain": "SLAB", "supplements": "Glucose", "time_h": ti, "od600": yi})
    df = pd.DataFrame(rows)
    df = transform_to_log_n_n0(df, value_col="od600", group_cols=["well"])
    params = fit_modified_gompertz_per_series(df, value_col="log_n_n0")
    preds = predict_modified_gompertz_per_series(df, params, value_col="log_n_n0")
    plot_and_save(
        preds, params, time_col="time_h", value_col="log_n_n0",
        group_cols=["well"], output_dir=tmp_path,
    )
    pngs = list((tmp_path / "plots").rglob("*.png"))
    assert len(pngs) == 2
