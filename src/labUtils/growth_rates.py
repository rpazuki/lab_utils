##############################################################
#  labUtils: A python package for fitting modified Gompertz  #
#  model and getting growth rates from microbial growth      #
#  and automation.                                           #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2025                          #
#  License: MIT                                              #
##############################################################

from math import e
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def gompertz(t, y0, A_0, mu_max, lam, clip_exp):
    """Zwietering-modified Gompertz model for microbial growth"""
    # safe version to avoid exp overflow
    if A_0 == 0:
        return np.full_like(t, y0)

    z = (mu_max * e / A_0) * (lam - t) + 1.0
    z = np.clip(z, -clip_exp, clip_exp)
    return y0 + A_0 * np.exp(-np.exp(z))


def transform_to_log_n_n0(
    df: pd.DataFrame,
    value_col: str = "od600",
    transformed_col: str = "log_n_n0",
    group_cols: list[str] = ["well"],
    OD_0_averaging_window: int = 1,  # number of initial points to average for OD_0
) -> pd.DataFrame:
    """Transforms the value_col to log(n/n0) per group defined by group_cols."""
    df_transformed = df.copy()
    df_transformed[transformed_col] = np.nan  # Initialize the new column with NaNs
    for keys, g in df_transformed.groupby(group_cols, sort=False):
        keys = (keys,) if not isinstance(keys, tuple) else keys
        y = g[value_col].to_numpy(dtype=float)
        m = np.isfinite(y)
        y_valid = y[m]
        if len(y_valid) == 0:
            continue
        n0 = float(np.nanmedian(np.mean(y_valid[0:OD_0_averaging_window])))  # baseline from the first valid point as N0
        n = y_valid
        with np.errstate(divide="ignore", invalid="ignore"):
            log_n_n0 = np.where(n0 > 0, np.log(n / n0), np.nan)
        df_transformed.loc[g.index[m], transformed_col] = log_n_n0

    return df_transformed


def fit_modified_gompertz_per_series(
    df: pd.DataFrame,
    time_col: str = "time_h",  # or "time_min" (then your mu_max units change accordingly)
    value_col: str = "od600",
    group_cols: list[str] = ["well"],
    clip_exp: float = 50.0,  # numerical safety against overflow in exp(exp(.))
    min_points: int = 5,  # require at least this many points to fit
    fixed_params: dict[str, float] | None = None,
) -> pd.DataFrame:
    """
    Fits the modified Gompertz model to each series in `df` and returns the fit parameters.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with time series data
    time_col : str
        Column name for time values (default: "time_h")
    value_col : str
        Column name for measured values (default: "od600")
    group_cols : List[str]
        Columns to group by (default: ["well"])
    clip_exp : float
        Numerical safety parameter against overflow in exp(exp(.)) (default: 50.0)
    min_points : int
        Minimum number of points required to fit (default: 5)
    fixed_params : Optional[Dict[str, float]]
        Dictionary specifying which parameters to fix and their values (default: None).
        Valid keys: 'y0', 'A', 'mu_max', 'lambda'
        Example: {'y0': 0.0, 'lambda': 2.0} fixes y0 to 0 and lambda to 2,
        while A and mu_max are fitted.

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per series (group), containing:
        - Group columns (e.g., well)
        - y0, A, mu_max, lambda: Fit parameters (fixed or fitted)
        - r2, rmse: Goodness of fit metrics
        - n: Number of data points
        - success: Whether fit was successful
        - message: Status or error message
        - Metadata columns from original data
    """

    # Normalize fixed_params
    if fixed_params is None:
        fixed_params = {}

    # Define parameter names and their order
    param_names = ["y0", "A", "mu_max", "lambda"]

    # Determine which parameters are free (to be fitted)
    free_param_names = [p for p in param_names if p not in fixed_params]

    def fit_one(gdf: pd.DataFrame) -> dict[str, Any]:
        """Unified fitting function that handles any combination of fixed parameters."""
        t = gdf[time_col].to_numpy(dtype=float)
        y = gdf[value_col].to_numpy(dtype=float)
        m = np.isfinite(t) & np.isfinite(y)
        t, y = t[m], y[m]
        out: dict[str, Any] = {"n": int(len(y))}

        # np.exp(0.1) is used as the upper limit, by assuming the data are log(n/n0) transformed
        if len(y) < min_points or (float(np.nanmax(y)) - float(np.nanmin(y)) <= np.exp(0.1)):
            out.update({"success": False, "message": "insufficient or flat data"})
            return out

        # Heuristic initial guesses for all parameters
        y0_guess = float(np.nanpercentile(y, 5))
        ymax = float(np.nanpercentile(y, 95))
        a_guess = max(ymax - y0_guess, 1e-3)

        if len(t) >= 3:
            dt = np.diff(t)
            dy = np.diff(y)
            slopes = dy / np.where(dt == 0, np.nan, dt)
            mu_guess = float(np.nanmax(slopes)) if np.isfinite(slopes).any() else 0.1
            mu_guess = max(mu_guess, 1e-3)
            k = int(np.nanargmax(slopes)) if np.isfinite(slopes).any() else 0
            lam_guess = float(t[min(k, len(t) - 1)])
        else:
            mu_guess = 0.2
            lam_guess = float(np.nanmedian(t))

        # Build initial guesses and bounds for free parameters only
        all_guesses = {"y0": y0_guess, "A": a_guess, "mu_max": mu_guess, "lambda": lam_guess}
        all_lower_bounds = {"y0": 0.0, "A": 1e-5, "mu_max": 1e-5, "lambda": 0.0}
        all_upper_bounds = {
            "y0": max(2.0, float(np.nanmax(y)) * 2),
            "A": max(1e-3, float(np.nanmax(y)) * 2),
            "mu_max": 10.0,
            "lambda": float(np.nanmax(t)) * 1.5 + 1,
        }

        p0 = [all_guesses[p] for p in free_param_names]
        lb = [all_lower_bounds[p] for p in free_param_names]
        ub = [all_upper_bounds[p] for p in free_param_names]

        # Create wrapper function that inserts fixed parameters
        def gompertz_wrapper(t_arr, *free_params):
            # Reconstruct full parameter set
            params = {}
            free_idx = 0
            for p in param_names:
                if p in fixed_params:
                    params[p] = fixed_params[p]
                else:
                    params[p] = free_params[free_idx]
                    free_idx += 1
            return gompertz(t_arr, params["y0"], params["A"], params["mu_max"], params["lambda"], clip_exp)

        try:
            result = curve_fit(gompertz_wrapper, t, y, p0=p0, bounds=(lb, ub), maxfev=20000)
            popt = result[0]

            # Reconstruct full parameter set from fitted and fixed values
            fitted_params = {}
            free_idx = 0
            for p in param_names:
                if p in fixed_params:
                    fitted_params[p] = fixed_params[p]
                else:
                    fitted_params[p] = popt[free_idx]
                    free_idx += 1

            yhat = gompertz(
                t, fitted_params["y0"], fitted_params["A"], fitted_params["mu_max"], fitted_params["lambda"], clip_exp
            )
            resid = y - yhat
            ss_res = float(np.sum(resid**2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            rmse = float(np.sqrt(ss_res / max(1, len(y) - len(popt))))

            out.update(
                {
                    "success": True,
                    "message": "ok",
                    "y0": fitted_params["y0"],
                    "A": fitted_params["A"],
                    "mu_max": fitted_params["mu_max"],
                    "lambda": fitted_params["lambda"],
                    "r2": r2,
                    "rmse": rmse,
                }
            )

            # predictions aligned to original (including any rows filtered by NaN check)
            t_full = gdf[time_col].to_numpy(dtype=float)
            out["pred"] = gompertz(
                t_full,
                fitted_params["y0"],
                fitted_params["A"],
                fitted_params["mu_max"],
                fitted_params["lambda"],
                clip_exp,
            )
        except (RuntimeError, ValueError) as ex:
            out.update({"success": False, "message": str(ex)[:200]})

        return out

    # Fit each series
    param_rows = []

    # Identify time-related columns that should be excluded from metadata
    time_related_cols = ["time_label", "time_h", "time_min", "time_h_int", "time_min_int"]

    # Identify metadata columns (all columns except group_cols, time_cols, value_col, well position cols)
    exclude_cols = set(group_cols + time_related_cols + [time_col, value_col])
    metadata_cols = [col for col in df.columns if col not in exclude_cols]

    for keys, g in df.groupby(group_cols, sort=False):
        keys = (keys,) if not isinstance(keys, tuple) else keys
        result = fit_one(g)
        row = {col: val for col, val in zip(group_cols, keys)}
        row.update({k: v for k, v in result.items() if k != "pred"})

        # Add metadata from the first row of the group
        first_row = g.iloc[0]
        for col in metadata_cols:
            if col in first_row.index:
                row[col] = first_row[col]

        param_rows.append(row)

    params_df = pd.DataFrame(param_rows)

    return params_df


def predict_modified_gompertz_per_series(
    df: pd.DataFrame,
    params_df: pd.DataFrame,
    time_col: str = "time_h",
    value_col: str = "od600",
    group_cols: list[str] = ["well"],
    clip_exp: float = 50.0,
    save_plot_data: bool = False,
    output_dir: str | Path | None = None,
    standard_deviation_column: str | None = None,
) -> pd.DataFrame:
    """
    Calculate predictions and residuals using fitted Gompertz parameters.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with time series data (same as used for fitting)
    params_df : pd.DataFrame
        DataFrame with fitted parameters from fit_modified_gompertz_per_series
    time_col : str
        Column name for time values (default: "time_h")
    value_col : str
        Column name for measured values (default: "od600")
    group_cols : List[str]
        Columns to group by (default: ["well"])
    clip_exp : float
        Numerical safety parameter against overflow (default: 50.0)
    save_plot_data : bool
        Whether to save plots of the data and fit (default: False)
    output_dir : Optional[str | Path]
        Directory to save plots (default: None)
    standard_deviation_column : Optional[str]
        Column name containing standard deviation values for error bars in plots (default: None)

    Returns
    -------
    pd.DataFrame
        DataFrame with same rows as input df plus columns:
        - od600_fit: Predicted values
        - residual: Difference between measured and predicted values
    """

    def gompertz_internal(t, y0, A, mu_max, lam):
        return gompertz(t, y0, A, mu_max, lam, clip_exp)

    preds_accum = []

    for keys, g in df.groupby(group_cols, sort=False):
        keys = (keys,) if not isinstance(keys, tuple) else keys

        # Get parameters for this group
        param_filter = params_df
        for col, val in zip(group_cols, keys):
            param_filter = param_filter[param_filter[col] == val]

        if len(param_filter) == 0:
            # No parameters found for this group, skip
            h = g.copy()
            h["od600_fit"] = np.nan
            h["residual"] = np.nan
            preds_accum.append(h)
            continue

        param_row = param_filter.iloc[0]

        # Generate predictions
        t = g[time_col].to_numpy(dtype=float)
        if param_row["success"]:
            pred = gompertz_internal(t, param_row["y0"], param_row["A"], param_row["mu_max"], param_row["lambda"])
        else:
            pred = np.full(len(g), np.nan)

        h = g.copy()
        h["od600_fit"] = pred
        h["residual"] = h[value_col] - h["od600_fit"]
        preds_accum.append(h)

    preds_df = pd.concat(preds_accum, ignore_index=True)

    if save_plot_data:
        plot_and_save(preds_df, params_df, time_col, value_col, group_cols, output_dir, standard_deviation_column)

    return preds_df


def plot_and_save(
    preds_df: pd.DataFrame,
    params_df: pd.DataFrame,
    time_col: str = "time_h",  # or "time_min" (then your mu_max units change accordingly)
    value_col: str = "od600",
    group_cols: list[str] = ["well"],
    output_dir: str | Path | None = None,
    standard_deviation_column: str | None = None,
):
    """Plots and saves growth curves with fitted model for each series.

    Parameters
    ----------
    standard_deviation_column : Optional[str]
        Column name containing standard deviation values for error bars.
        If provided, error bars will be added to the scatter plot.
    """
    # Only merge the fit parameters, not metadata columns that already exist in preds_df
    fit_cols = group_cols + ["y0", "A", "mu_max", "lambda", "r2", "rmse", "n", "success", "message"]
    params_to_merge = params_df[[col for col in fit_cols if col in params_df.columns]]
    df_combined = preds_df.merge(params_to_merge, on=group_cols, how="left")

    pred_col: str = "od600_fit"

    # Get unique groups to avoid duplicates
    unique_groups = df_combined[group_cols].drop_duplicates()

    for _, (_, group_row) in enumerate(unique_groups.iterrows()):
        # Filter data for this specific group
        mask = True
        for col in group_cols:
            mask = mask & (df_combined[col] == group_row[col])
        g = df_combined[mask]

        if g["success"].iloc[0] is False or pd.isna(g["mu_max"].iloc[0]):
            continue

        keys = tuple(group_row[col] for col in group_cols)
        if len(keys) == 1:
            keys = keys[0]

        t = g[time_col].to_numpy(dtype=float)
        y = g[value_col].to_numpy(dtype=float)
        y_hat = g[pred_col].to_numpy(dtype=float)

        # Get standard deviation if provided
        yerr = None
        if standard_deviation_column is not None and standard_deviation_column in g.columns:
            yerr = g[standard_deviation_column].to_numpy(dtype=float)

        plt.figure(figsize=(8, 5))
        if yerr is not None:
            plt.errorbar(
                t,
                y,
                yerr=yerr,
                fmt=".",
                label=f"{group_cols[0]}={keys}"
                if not isinstance(keys, tuple)
                else ", ".join(f"{col}={val}" for col, val in zip(group_cols, keys)),
                capsize=3,
                capthick=1,
                elinewidth=1,
            )
        else:
            plt.scatter(
                t,
                y,
                label=f"{group_cols[0]}={keys}"
                if not isinstance(keys, tuple)
                else ", ".join(f"{col}={val}" for col, val in zip(group_cols, keys)),
                marker=".",
            )
        plt.plot(t, y_hat, label="Predicted", linestyle="--", color="orange")
        plt.xlabel("Time (h)")
        plt.ylabel(r"$\ln(OD/OD_0)$")

        # Get the first key value for filename and title
        first_key = keys if not isinstance(keys, tuple) else keys[0]
        plt.title(
            f"Growth Curve for  well {first_key} growth_rate:{g['mu_max'].iloc[0]:.4f} 1/h \n"
            f"y0={g['y0'].iloc[0]:.3f}, A={g['A'].iloc[0]:.3f}, lambda={g['lambda'].iloc[0]:.3f}"
        )
        plt.legend()
        # plt.ylim(bottom=0)
        plt.grid(True)

        # Save figure instead of showing
        filename = f"growth_curve_well_{first_key}.png"
        if output_dir is not None:
            output_dir_plots = Path(output_dir) / "plots"
            output_dir_plots.mkdir(parents=True, exist_ok=True)
            filepath = output_dir_plots / filename
        else:
            filepath = filename
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
        print(f"Saved plot: {filepath}")
        plt.close()  # Close figure to free memory
