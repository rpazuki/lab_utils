import numpy as np

try:
    from statsmodels.regression.linear_model import yule_walker
except ImportError:  # pragma: no cover - exercised only when optional dep missing
    yule_walker = None


def spectral_density_zero(x, order=None):
    """
    Estimate the spectral density at frequency zero for a 1-D time series.

    Fits an AR(p) model using the Yule-Walker equations and returns:
        S(0) = σ² / (1 − φ₁ − φ₂ − ⋯ − φ_p)²

    Parameters
    ----------
    x     : array-like, 1-D time series (typically a chain window)
    order : int or None; if None uses floor(log(n)) as rule-of-thumb

    Returns
    -------
    float : estimated spectral density at zero frequency
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


def geweke_zscore(x, first=0.1, last=0.5):
    """
    Compute a single Geweke Z-score for a 1-D MCMC chain.

    Compares the mean of the first `first` fraction (window A) with the mean
    of the last `last` fraction (window B). Standard errors are estimated
    from the spectral density at zero to correct for autocorrelation.

    H₀: E[θ]_A = E[θ]_B  (chain has reached stationarity)
    Reject H₀ if |Z| > 1.96  (α = 0.05)

    Parameters
    ----------
    x     : array-like, MCMC chain of length N
    first : float ∈ (0, 1), default 0.1 — fraction used for window A
    last  : float ∈ (0, 1), default 0.5 — fraction used for window B

    Returns
    -------
    z : float
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


def geweke_sweep(x, first=0.1, last=0.5, intervals=20):
    """
    Multi-interval Geweke sweep.

    Repeats the Geweke test at `intervals` evenly-spaced starting positions
    from 0 to floor((1 − last) × N). Returns arrays suitable for plotting.

    Parameters
    ----------
    x         : array-like, full MCMC chain
    first     : float, fraction for window A (default 0.1)
    last      : float, fraction for window B (default 0.5)
    intervals : int, number of starting positions (default 20)

    Returns
    -------
    starts  : ndarray, start indices
    zscores : ndarray, corresponding Z-scores
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
    sweeps,
    labels=None,
    colors=None,
    title="Geweke Convergence Diagnostic - Multi-Interval Sweep",
    figsize_per_plot=(6, 4),
    save_path=None,
    show=True,
):
    """
    Plot one or more Geweke sweep results.

    Parameters
    ----------
    sweeps : sequence of tuple(ndarray, ndarray)
        Each element must be ``(starts, zscores)`` as returned by ``geweke_sweep``.
    labels : sequence of str or None
        Optional subplot titles for each sweep.
    colors : sequence of matplotlib-compatible color or None
        Optional scatter color per sweep.
    title : str
        Figure-level title.
    figsize_per_plot : tuple(float, float)
        Per-subplot width and height in inches.
    save_path : str or None
        If provided, saves the figure to this path.
    show : bool
        Whether to call ``plt.show()``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : ndarray of matplotlib.axes.Axes
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