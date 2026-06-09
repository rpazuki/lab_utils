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



def plot_geweke_diagnostics(
    chains: dict[str, np.ndarray] | list[np.ndarray],
    labels: dict[str, str] | list[str] | None = None,
    colors: dict[str, str] | list[str] | None = None,
    first: float = 0.1,
    last: float = 0.5,
    intervals: int = 20,
    title: str = "Geweke Convergence Diagnostic - Multi-Interval Sweep",
    figsize_per_plot: tuple[float, float] = (6, 4),
    output_file: str | None = None,
    show_plot: bool = True,
) -> tuple[Any, NDArray[Any]]:
    """Plot Geweke convergence diagnostics for multiple chains.

    This is a convenience function that computes Geweke sweeps internally and
    then plots them. It follows the same interface pattern as plot_ess_diagnostics()
    and plot_ipsrf_diagnostics().

    Parameters
    ----------
    chains : dict[str, np.ndarray] | list[np.ndarray]
        Dictionary mapping chain identifiers to chain data, or a list of chains.
    labels : dict[str, str] | list[str] | None, default=None
        Optional labels for each chain. If None, creates labels as "chain1", "chain2", etc.
        Can be a dict mapping chain IDs to labels, or a list of labels in chain order.
    colors : dict[str, str] | list[str] | None, default=None
        Optional matplotlib-compatible colors for each chain. If None, uses default palette.
        Can be a dict mapping chain IDs to colors, or a list of colors in chain order.
    first : float, default=0.1
        Fraction of each chain used for the leading comparison window in Geweke test.
    last : float, default=0.5
        Fraction of each chain used for the trailing comparison window in Geweke test.
    intervals : int, default=20
        Number of starting positions evaluated across each chain in the sweep.
    title : str, optional
        Figure-level title (default: "Geweke Convergence Diagnostic - Multi-Interval Sweep")
    figsize_per_plot : tuple[float, float], optional
        Per-subplot width and height in inches (default: (6, 4))
    output_file : str | None, optional
        If provided, save the figure to this path (default: None)
    show_plot : bool, optional
        Whether to call ``plt.show()`` (default: True)

    Returns
    -------
    tuple[Any, NDArray[Any]]
        The matplotlib figure and array of axes.

    Raises
    ------
    ValueError
        If chains is empty or if ``first + last >= 1.0``.
    """
    import matplotlib.pyplot as plt
    # Normalize chains to dict format
    if isinstance(chains, list):
        chain_ids = [f"chain{i}" for i in range(len(chains))]
        chains_dict = {cid: c for cid, c in zip(chain_ids, chains)}
    else:
        chain_ids = list(chains.keys())
        chains_dict = chains

    if not chains_dict:
        raise ValueError("chains must contain at least one chain.")

    # Normalize labels
    if labels is None:
        labels_dict = {cid: f"chain{i+1}" for i, cid in enumerate(chain_ids)}
    elif isinstance(labels, dict):
        labels_dict = labels
    else:
        labels_dict = {cid: l for cid, l in zip(chain_ids, labels)}

    # Normalize colors
    default_colors = ["steelblue", "seagreen", "tomato", "gold", "purple", "coral",
                      "brown", "pink", "gray", "olive"]
    if colors is None:
        colors_dict = {cid: default_colors[i % len(default_colors)]
                       for i, cid in enumerate(chain_ids)}
    elif isinstance(colors, dict):
        colors_dict = colors
    else:
        colors_dict = {cid: c for cid, c in zip(chain_ids, colors)}

    # Compute sweeps for each chain
    sweeps = []
    labels = []
    colors = []

    for cid in chain_ids:
        chain = chains_dict[cid]
        starts, zscores = geweke_sweep(chain, first=first, last=last, intervals=intervals)
        sweeps.append((starts, zscores))
        labels.append(labels_dict.get(cid, cid))
        colors.append(colors_dict.get(cid))
    
    n_sweeps = len(sweeps)
    
    # Default color palette
    default_colors = ["steelblue", "seagreen", "tomato", "gold", "purple", "coral",
                      "brown", "pink", "gray", "olive"]
    
    # Set up labels
    if labels is None:
        labels = [f"sweep{i+1}" for i in range(n_sweeps)]
    else:
        labels = list(labels)
        if len(labels) != n_sweeps:
            raise ValueError("labels length must match sweeps length.")
    
    # Set up colors
    if colors is None:
        colors = [default_colors[i % len(default_colors)] for i in range(n_sweeps)]
    else:
        colors = list(colors)
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

    if output_file is not None:
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
    if show_plot:
        plt.show()

    return fig, axes