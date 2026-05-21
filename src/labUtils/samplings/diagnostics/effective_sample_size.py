from typing import TYPE_CHECKING

import numpy as np
# import arviz as az

if TYPE_CHECKING:
    from matplotlib.figure import Figure


# ─────────────────────────────────────────────────────────────────────────────
# §7.1  From-scratch ESS: FFT autocorrelation + Geyer IPS
# ─────────────────────────────────────────────────────────────────────────────

def autocorr_fft(x: np.ndarray) -> np.ndarray:
    """
    Normalised autocorrelation function ρ(k) for k = 0, 1, …, N-1.

    Uses the FFT-based Wiener-Khinchin approach for O(N log N) complexity.
    Zero-pads to length 2N to avoid circular correlation wrap-around.

    Parameters
    ----------
    x : np.ndarray
        1-D array-like (mean-centred internally)

    Returns
    -------
    np.ndarray
        Autocorrelation array of length N, with rho[0] = 1.0 by construction
    """
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    # Pad to 2N for linear (non-circular) correlation
    f = np.fft.rfft(x, n=2 * n)
    acf_full = np.fft.irfft(f * np.conj(f))[:n]  # keep first N lags
    acf_full /= acf_full[0]                        # normalise: rho[0] = 1
    return acf_full


def geyer_ips(rho: np.ndarray) -> float:
    """
    Geyer's Initial Positive Sequence estimator for τ_int (integrated autocorrelation time).

    Forms consecutive pairs Γ_k = ρ(2k-1) + ρ(2k) and sums until
    the first negative Γ_k. The resulting τ_int is consistent and
    always ≥ 1.

    Parameters
    ----------
    rho : np.ndarray
        Normalised ACF array with rho[0] = 1.0

    Returns
    -------
    float
        Integrated autocorrelation time τ_int ≥ 1
    """
    n = len(rho)
    tau = 1.0
    k = 1
    while 2 * k + 1 < n:
        gamma_k = rho[2 * k - 1] + rho[2 * k]
        if gamma_k < 0.0:
            break
        tau += 2.0 * gamma_k
        k += 1
    return max(1.0, tau)


def ess_from_chain(x: np.ndarray) -> float:
    """
    Effective sample size for a 1-D MCMC chain.

    Computed as ESS = N / τ_int, where τ_int is the integrated
    autocorrelation time estimated via FFT autocorrelation + Geyer IPS.

    Parameters
    ----------
    x : np.ndarray
        1-D array representing an MCMC chain of length N

    Returns
    -------
    float
        Effective sample size (ESS)
    """
    x = np.asarray(x, dtype=float)
    rho = autocorr_fft(x)
    tau = geyer_ips(rho)
    return len(x) / tau

def ess_sweep(
    chain: np.ndarray,
    subchains_length: int = 1000,
    fraction: float | None = None,
) -> np.ndarray:
    """
    Compute ESS for fixed-length subchains within a chain.

    The chain is split into consecutive non-overlapping subchains of equal length,
    and ESS is computed independently for each subchain. If ``fraction`` is
    provided, subchain length is set to ``floor(fraction * N)`` where ``N`` is
    chain length. Any remainder samples at the end that do not fill a complete
    subchain are ignored.

    Parameters
    ----------
    chain : np.ndarray
        1-D array representing an MCMC chain of length N.
    subchains_length : int, optional
        Length of each subchain. Must be positive (default: 1000).
    fraction : float | None, optional
        Fraction in ``(0, 1]`` used to derive subchain length as
        ``int(fraction * N)``. If provided, overrides ``subchains_length``
        (default: None).

    Returns
    -------
    np.ndarray
        Array of ESS values, one per full subchain.

    Raises
    ------
    ValueError
        If ``chain`` is empty, if ``subchains_length <= 0``, or if ``fraction``
        is not in ``(0, 1]``.
    """
    chain = np.asarray(chain, dtype=float)
    n = len(chain)

    if n == 0:
        raise ValueError("chain must contain at least one sample.")

    if subchains_length <= 0:
        raise ValueError("subchains_length must be a positive integer.")

    if fraction is not None:
        if not (0.0 < fraction <= 1.0):
            raise ValueError("fraction must be in the interval (0, 1].")
        subchains_length = int(n * fraction)

    subchains_length = max(1, min(subchains_length, n))
    subchains_number = n // subchains_length

    # split the chain as fixed length subchains and compute ESS for each
    ess_values = np.array([
        ess_from_chain(chain[i * subchains_length : (i + 1) * subchains_length])
        for i in range(subchains_number)
    ])
    return ess_values



# ─────────────────────────────────────────────────────────────────────────────
# Visualization function
# ─────────────────────────────────────────────────────────────────────────────

def plot_ess_diagnostics(
    chains: dict[str, np.ndarray] | list[np.ndarray],
    labels: dict[str, str] | list[str] | None = None,
    colors: dict[str, str] | list[str] | None = None,
    max_trace_draws: int = 500,
    max_acf_lags: int = 200,
    figsize: tuple[int, int] = (15, 8),
    output_file: str | None = None,
    show_plot: bool = True,
) -> "Figure":
    """
    Plot trace plots and autocorrelation functions for MCMC chains.

    Creates a 2×N subplot grid where N is the number of chains. Top row shows
    trace plots (first draws), bottom row shows ACF with ESS estimates.

    Parameters
    ----------
    chains : dict[str, np.ndarray] | list[np.ndarray]
        Either a dictionary mapping chain IDs to 1-D arrays, or a list of 1-D arrays.
        Example: {"chain_a": data1, "chain_b": data2} or [data1, data2]
    labels : dict[str, str] | list[str] | None, optional
        Display labels for each chain. If None, creates labels as "chain1", "chain2", etc.
        If dict, maps chain IDs to labels. If list, must match chain order (default: None)
    colors : dict[str, str] | list[str] | None, optional
        Colors for each chain. If None, uses default palette.
        If dict, maps chain IDs to colors. If list, must match chain order (default: None)
    max_trace_draws : int, optional
        Maximum number of draws to show in trace plots (default: 500)
    max_acf_lags : int, optional
        Maximum lag to show in ACF plots (default: 200)
    figsize : tuple[int, int], optional
        Figure size as (width, height) (default: (15, 8))
    output_file : str | None, optional
        Path to save figure. If None, does not save (default: None)
    show_plot : bool, optional
        Whether to display the plot (default: True)

    Returns
    -------
    Figure
        The created figure object
    """
    import matplotlib.pyplot as plt

    # Convert to dict format for consistency
    if isinstance(chains, list):
        chains = {f"chain{i}": c for i, c in enumerate(chains)}
    
    n_chains = len(chains)
    chain_ids = list(chains.keys())
    
    # Default color palette
    default_colors = ["steelblue", "seagreen", "tomato", "gold", "purple", "coral",
                      "brown", "pink", "gray", "olive"]
    
    # Set up labels
    if labels is None:
        labels = {cid: f"chain{i+1}" for i, cid in enumerate(chain_ids)}
    elif isinstance(labels, list):
        labels = {cid: l for cid, l in zip(chain_ids, labels)}
    
    # Set up colors
    if colors is None:
        colors = {cid: default_colors[i % len(default_colors)] for i, cid in enumerate(chain_ids)}
    elif isinstance(colors, list):
        colors = {cid: c for cid, c in zip(chain_ids, colors)}
    
    fig, axes = plt.subplots(2, n_chains, figsize=figsize)
    
    # Handle single chain case (axes is 1D)
    if n_chains == 1:
        axes = axes.reshape(2, 1)
    
    for col, chain_id in enumerate(chain_ids):
        chain = chains[chain_id]
        label = labels.get(chain_id, chain_id)
        color = colors.get(chain_id, "steelblue")
        n = len(chain)
        
        # Top row — trace plot (first max_trace_draws draws)
        ax_tr = axes[0, col]
        trace_end = min(max_trace_draws, n)
        ax_tr.plot(chain[:trace_end], lw=0.7, color=color, alpha=0.85)
        ax_tr.set_title(label, fontsize=10)
        ax_tr.set_xlabel("Iteration", fontsize=9)
        ax_tr.set_ylabel(r"$\theta$", fontsize=9)
        
        # Bottom row — ACF (first max_acf_lags lags)
        ax_ac = axes[1, col]
        rho = autocorr_fft(chain)
        lags = np.arange(min(max_acf_lags, len(rho)))
        ax_ac.bar(lags, rho[: len(lags)], width=1, color=color, alpha=0.7)
        ax_ac.axhline(0, color="black", lw=0.8)
        ci = 1.96 / np.sqrt(n)
        ax_ac.axhline(ci, color="red", lw=1.0, ls="--", label=r"±1.96/$\sqrt{N}$")
        ax_ac.axhline(-ci, color="red", lw=1.0, ls="--")
        ess_val = ess_from_chain(chain)
        ax_ac.set_title(f"ACF  |  ESS = {ess_val:.0f} / {n}", fontsize=10)
        ax_ac.set_xlabel("Lag $k$", fontsize=9)
        ax_ac.set_ylabel(r"$\hat{\rho}(k)$", fontsize=9)
        ax_ac.legend(fontsize=8, loc="upper right")
    
    plt.suptitle(
        "Effective Sample Size — Trace Plots and Autocorrelation Functions",
        fontsize=13,
    )
    plt.tight_layout()
    
    if output_file is not None:
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
    
    if show_plot:
        plt.show()
    
    return fig

