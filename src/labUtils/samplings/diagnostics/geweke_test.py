from __future__ import annotations

from os import PathLike
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

try:
    from statsmodels.regression.linear_model import yule_walker
except ImportError:  # pragma: no cover - exercised only when optional dep missing
    yule_walker = None


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int_]


def spectral_density_zero(x: ArrayLike, order: int | None = None) -> float:
    """Estimate the zero-frequency spectral density for a one-dimensional series.

    Fits an AR(p) model with the Yule-Walker equations and returns the
    zero-frequency estimate ``S(0) = sigma^2 / (1 - phi_1 - ... - phi_p)^2``.
    When the series is too short, the autoregressive fit is unstable, or
    ``statsmodels`` is unavailable, this falls back to the sample variance.

    Parameters
    ----------
    x : ArrayLike
        One-dimensional time series, typically a window from an MCMC chain.
    order : int | None, default=None
        Autoregressive order. If omitted, uses ``round(log(n))`` with a minimum
        value of 1.

    Returns
    -------
    float
        Estimated spectral density at zero frequency.
    """
    x = np.asarray(x, dtype=float)
    x = x - x.mean()                              # demean
    n = len(x)
    if n < 4:
        return float(np.var(x, ddof=1))           # too short for AR fit
    if order is None:
        order = max(1, int(np.round(np.log(n))))  # rule-of-thumb
    if yule_walker is None:
        return float(np.var(x, ddof=1))
    try:
        phi, sigma2 = yule_walker(x, order=order, method="mle")
        denominator = 1.0 - phi.sum()
        if abs(denominator) < 1e-8:               # near-unit-root: fall back
            return float(np.var(x, ddof=1))
        return float(sigma2 / denominator ** 2)
    except Exception:
        return float(np.var(x, ddof=1))           # fallback to sample variance


def geweke_zscore(x: ArrayLike, first: float = 0.1, last: float = 0.5) -> float:
    """Compute a Geweke Z-score for a one-dimensional MCMC chain.

    Compares the mean of the first ``first`` fraction of the chain with the mean
    of the last ``last`` fraction. Standard errors are estimated
    from the spectral density at zero to correct for autocorrelation.

    Under stationarity, the two windows should have matching expectations.
    Values with absolute magnitude above roughly 1.96 are commonly treated as
    evidence against convergence at the 5% level.

    Parameters
    ----------
    x : ArrayLike
        One-dimensional Markov chain sample.
    first : float, default=0.1
        Fraction of the chain used for the leading window.
    last : float, default=0.5
        Fraction of the chain used for the trailing window.

    Returns
    -------
    float
        Geweke Z-score comparing the two windows.

    Raises
    ------
    ValueError
        If ``first + last >= 1.0`` because the windows would overlap.
    """
    x = np.asarray(x, dtype=float)
    if first + last >= 1.0:
        raise ValueError(
            f"first ({first}) + last ({last}) must be < 1.0 "
            "to ensure window independence."
        )
    n = len(x)
    n_A = int(first * n)
    n_B = int(last * n)
    window_A = x[:n_A]
    window_B = x[n - n_B:]                        # last `last` fraction

    mean_A, mean_B = window_A.mean(), window_B.mean()
    S0_A = spectral_density_zero(window_A)
    S0_B = spectral_density_zero(window_B)

    se = np.sqrt(S0_A / n_A + S0_B / n_B)
    if se == 0:
        return 0.0
    return float((mean_A - mean_B) / se)


def geweke_sweep(
    x: ArrayLike,
    first: float = 0.1,
    last: float = 0.5,
    intervals: int = 20,
) -> tuple[IntArray, FloatArray]:
    """Run the Geweke diagnostic over multiple starting positions.

    Repeats the Geweke test at ``intervals`` evenly spaced starting positions
    from 0 to ``floor((1 - last) * N)`` and returns arrays suitable for plotting.

    Parameters
    ----------
    x : ArrayLike
        Full one-dimensional MCMC chain.
    first : float, default=0.1
        Fraction for the leading comparison window.
    last : float, default=0.5
        Fraction for the trailing comparison window.
    intervals : int, default=20
        Number of starting positions evaluated across the chain.

    Returns
    -------
    tuple[IntArray, FloatArray]
        Start indices and corresponding Z-scores.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    max_start = int(n * (1.0 - last))
    starts = np.linspace(0, max_start, intervals).astype(int)
    zscores = np.array(
        [geweke_zscore(x[s:], first=first, last=last) for s in starts]
    )
    return starts, zscores


def plot_geweke_sweeps(
    sweeps: Sequence[tuple[ArrayLike, ArrayLike]],
    labels: Sequence[str] | None = None,
    colors: Sequence[Any | None] | None = None,
    title: str = "Geweke Convergence Diagnostic - Multi-Interval Sweep",
    figsize_per_plot: tuple[float, float] = (6, 4),
    save_path: str | PathLike[str] | None = None,
    show: bool = True,
) -> tuple[Any, NDArray[Any]]:
    """Plot one or more Geweke sweep results.

    Parameters
    ----------
    sweeps : Sequence[tuple[ArrayLike, ArrayLike]]
        Sequence of ``(starts, zscores)`` pairs as returned by ``geweke_sweep``.
    labels : Sequence[str] | None, default=None
        Optional subplot titles for each sweep.
    colors : Sequence[Any | None] | None, default=None
        Optional matplotlib-compatible color for each sweep.
    title : str
        Figure-level title.
    figsize_per_plot : tuple[float, float], default=(6, 4)
        Per-subplot width and height in inches.
    save_path : str | PathLike[str] | None, default=None
        If provided, save the figure to this path.
    show : bool, default=True
        Whether to call ``plt.show()``.

    Returns
    -------
    tuple[Any, NDArray[Any]]
        The matplotlib figure and array of axes.

    Raises
    ------
    ValueError
        If no sweeps are provided or if ``labels`` or ``colors`` lengths do not
        match the number of sweeps.
    """
    import matplotlib.pyplot as plt

    if not sweeps:
        raise ValueError("sweeps must contain at least one (starts, zscores) result.")

    n_sweeps = len(sweeps)
    if labels is None:
        labels = [f"Sweep {i + 1}" for i in range(n_sweeps)]
    if colors is None:
        colors = [None] * n_sweeps

    if len(labels) != n_sweeps:
        raise ValueError("labels length must match sweeps length.")
    if len(colors) != n_sweeps:
        raise ValueError("colors length must match sweeps length.")

    fig, axes = plt.subplots(
        1,
        n_sweeps,
        figsize=(figsize_per_plot[0] * n_sweeps, figsize_per_plot[1]),
        sharey=True,
    )
    axes = np.atleast_1d(axes)

    for ax, (starts, zscores), label, color in zip(axes, sweeps, labels, colors):
        ax.axhspan(-1.96, 1.96, alpha=0.15, color="green", label="95% acceptance band")
        ax.axhline(0, color="gray", lw=0.8, ls="-")
        ax.axhline(1.96, color="red", lw=1.2, ls="--", label="+/-1.96 threshold")
        ax.axhline(-1.96, color="red", lw=1.2, ls="--")
        ax.scatter(starts, zscores, color=color, s=60, zorder=4)
        ax.set_xlabel("Start iteration", fontsize=11)
        ax.set_ylabel("Geweke Z-score", fontsize=11)
        ax.set_title(label, fontsize=11)
        ax.legend(fontsize=9)

    fig.suptitle(title, fontsize=13, y=1.02)
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()

    return fig, axes