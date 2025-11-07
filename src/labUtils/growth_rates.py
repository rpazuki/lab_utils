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
from typing import Any, Dict, List, Optional, Tuple

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

def transform_to_log_n_n0(df: pd.DataFrame,
                           value_col: str = "od600",
                           group_cols: List[str] = ["well"],
                           OD_0_averaging_window: int = 1,  # number of initial points to average for OD_0
                           ) -> pd.DataFrame:
    """Transforms the value_col to log(n/n0) per group defined by group_cols."""
    df_transformed = df.copy()
    df_transformed["log_n_n0"] = np.nan # Initialize the new column with NaNs
    for keys, g in df_transformed.groupby(group_cols, sort=False):
        keys = (keys,) if not isinstance(keys, tuple) else keys
        y = g[value_col].to_numpy(dtype=float)
        m = np.isfinite(y)
        y_valid = y[m]
        if len(y_valid) == 0:
            continue
        n0 = float(np.nanmedian(np.mean(y_valid[0:OD_0_averaging_window])))  # baseline from the first valid point as N0
        n = y_valid
        with np.errstate(divide='ignore', invalid='ignore'):
            log_n_n0 = np.where(n0 > 0, np.log(n / n0), np.nan)
        df_transformed.loc[g.index[m], "log_n_n0"] = log_n_n0

    return df_transformed

def fit_modified_gompertz_per_series(
    df: pd.DataFrame,
    time_col: str = "time_h",        # or "time_min" (then your mu_max units change accordingly)
    value_col: str = "od600",
    group_cols: List[str] = ["well"],
    clip_exp: float = 50.0,           # numerical safety against overflow in exp(exp(.))
    min_points: int = 5,             # require at least this many points to fit
    save_plot_data: bool = False,    # whether to save plot of the data and fit
    output_dir: Optional[str | Path] = None,
    y_0_fixed: Optional[float] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fits the modified Gompertz model to each series in `df` and returns:
      - params_df: one row per series (group), with y0, A, mu_max, lambda, r2, rmse, n, success, message
      - preds_df: same rows as input df + columns `od600_fit` and `residual`

    The time column should be numeric (e.g., hours). If you use minutes, adjust interpretation of mu_max.
    """

    def gompertz_internal(t, y0, A, mu_max, lam):
        return gompertz(t, y0, A, mu_max, lam, clip_exp)


    def fit_one(gdf: pd.DataFrame) -> Dict[str, Any]:
        t = gdf[time_col].to_numpy(dtype=float)
        y = gdf[value_col].to_numpy(dtype=float)
        m = np.isfinite(t) & np.isfinite(y)
        t, y = t[m], y[m]
        out: Dict[str, Any] = {"n": int(len(y))}
        if len(y) < min_points or (float(np.nanmax(y)) - float(np.nanmin(y)) <= np.exp(.1)):
            out.update({"success": False, "message": "insufficient or flat data"})
            return out

        # Heuristic initial guesses
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
            lam_guess = float(t[min(k, len(t)-1)])
        else:
            mu_guess = 0.2
            lam_guess = float(np.nanmedian(t))

        p0 = [y0_guess, a_guess, mu_guess, lam_guess]

        # Soft bounds to stabilize fitting; adjust if needed
        lb = [0.0,   1e-5, 1e-5, 0.0]
        ub = [max(2.0, float(np.nanmax(y))*2),
              max(1e-3, float(np.nanmax(y))*2),
              10.0,
              float(np.nanmax(t))*1.5 + 1]

        try:
            result = curve_fit(gompertz_internal, t, y, p0=p0, bounds=(lb, ub), maxfev=20000)
            popt = result[0]
            yhat = gompertz_internal(t, *popt)
            resid = y - yhat
            ss_res = float(np.sum(resid**2))
            ss_tot = float(np.sum((y - np.mean(y))**2))
            r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
            rmse = float(np.sqrt(ss_res / max(1, len(y) - len(popt))))

            out.update({
                "success": True, "message": "ok",
                "y0": popt[0], "A": popt[1], "mu_max": popt[2], "lambda": popt[3],
                "r2": r2, "rmse": rmse
            })
            # predictions aligned to original (including any rows filtered by NaN check)
            out["pred"] = gompertz_internal(gdf[time_col].to_numpy(dtype=float), *popt)
        except (RuntimeError, ValueError) as ex:
            out.update({"success": False, "message": str(ex)[:200]})
        return out

    def gompertz_internal_fixed(t, A, mu_max, lam):
        return gompertz(t, y_0_fixed, A, mu_max, lam, clip_exp)

    def fit_one_fixed(gdf: pd.DataFrame) -> Dict[str, Any]:
        t = gdf[time_col].to_numpy(dtype=float)
        y = gdf[value_col].to_numpy(dtype=float)
        m = np.isfinite(t) & np.isfinite(y)
        t, y = t[m], y[m]
        out: Dict[str, Any] = {"n": int(len(y))}
        if len(y) < min_points or (float(np.nanmax(y)) - float(np.nanmin(y)) <= np.exp(.1)):
            out.update({"success": False, "message": "insufficient or flat data"})
            return out

        # Heuristic initial guesses
        ymax = float(np.nanpercentile(y, 95))
        a_guess = max(ymax - y_0_fixed, 1e-3)

        if len(t) >= 3:
            dt = np.diff(t)
            dy = np.diff(y)
            slopes = dy / np.where(dt == 0, np.nan, dt)
            mu_guess = float(np.nanmax(slopes)) if np.isfinite(slopes).any() else 0.1
            mu_guess = max(mu_guess, 1e-3)
            k = int(np.nanargmax(slopes)) if np.isfinite(slopes).any() else 0
            lam_guess = float(t[min(k, len(t)-1)])
        else:
            mu_guess = 0.2
            lam_guess = float(np.nanmedian(t))

        p0 = [a_guess, mu_guess, lam_guess]

        # Soft bounds to stabilize fitting; adjust if needed
        lb = [1e-5, 1e-5, 0.0]
        ub = [max(1e-3, float(np.nanmax(y))*2),
              10.0,
              float(np.nanmax(t))*1.5 + 1]

        try:
            result = curve_fit(gompertz_internal_fixed, t, y, p0=p0, bounds=(lb, ub), maxfev=20000)
            popt = result[0]
            yhat = gompertz_internal_fixed(t, *popt)
            resid = y - yhat
            ss_res = float(np.sum(resid**2))
            ss_tot = float(np.sum((y - np.mean(y))**2))
            r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
            rmse = float(np.sqrt(ss_res / max(1, len(y) - len(popt))))

            out.update({
                "success": True, "message": "ok",
                "y0": y_0_fixed, "A": popt[0], "mu_max": popt[1], "lambda": popt[2],
                "r2": r2, "rmse": rmse
            })
            # predictions aligned to original (including any rows filtered by NaN check)
            out["pred"] = gompertz_internal_fixed(gdf[time_col].to_numpy(dtype=float), *popt)
        except (RuntimeError, ValueError) as ex:
            out.update({"success": False, "message": str(ex)[:200]})
        return out

    # Fit each series
    param_rows = []
    preds_accum = []
    for keys, g in df.groupby(group_cols, sort=False):
        keys = (keys,) if not isinstance(keys, tuple) else keys
        if y_0_fixed is not None:
            result = fit_one_fixed(g)
        else:
            result = fit_one(g)
        row = {col: val for col, val in zip(group_cols, keys)}
        row.update({k: v for k, v in result.items() if k != "pred"})
        param_rows.append(row)

        h = g.copy()
        h["od600_fit"] = result.get("pred", np.full(len(g), np.nan))
        h["residual"] = h[value_col] - h["od600_fit"]
        preds_accum.append(h)

    params_df = pd.DataFrame(param_rows)
    preds_df = pd.concat(preds_accum, ignore_index=True)

    if save_plot_data:
        plot_and_save(preds_df, params_df, time_col, value_col, group_cols, output_dir)
    return params_df, preds_df


def plot_and_save(preds_df: pd.DataFrame,
                  params_df: pd.DataFrame,
                  time_col: str = "time_h",        # or "time_min" (then your mu_max units change accordingly)
                  value_col: str = "od600",
                  group_cols: List[str] = ["well"],
                  output_dir: Optional[str | Path] = None):
    """Plots and saves growth curves with fitted model for each series."""
    df_combined = preds_df.merge(params_df, on="well", how="left")
    pred_col: str = "od600_fit"
    for keys, g in df_combined.groupby(group_cols, sort=False):
        if g['success'].iloc[0] is False or pd.isna(g['mu_max'].iloc[0]):
            continue
        keys = (keys,) if not isinstance(keys, tuple) else keys
        t = g[time_col].to_numpy(dtype=float)
        y = g[value_col].to_numpy(dtype=float)
        y_hat = g[pred_col].to_numpy(dtype=float)



        plt.figure(figsize=(8, 5))
        plt.scatter(t, y, label=", ".join(f"{col}={val}" for col, val in zip(group_cols, keys)), marker='o')
        plt.plot(t, y_hat, label="Predicted", linestyle="--", color="orange")
        plt.xlabel("Time (h)")
        plt.ylabel(r"$\ln(OD/OD_0)$")
        plt.title(f"Growth Curve for  well {keys[0]} growth_rate:{g['mu_max'].iloc[0]:.4f} 1/h \n"
                  f"y0={g['y0'].iloc[0]:.3f}, A={g['A'].iloc[0]:.3f}, lambda={g['lambda'].iloc[0]:.3f}")
        plt.legend()
        # plt.ylim(bottom=0)
        plt.grid(True)

        # Save figure instead of showing
        filename = f"growth_curve_well_{keys[0]}.png"
        if output_dir is not None:
            output_dir_plots = Path(output_dir) / "plots"
            output_dir_plots.mkdir(parents=True, exist_ok=True)
            filepath = output_dir_plots / filename
        else:
            filepath = filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
        plt.close()  # Close figure to free memory
