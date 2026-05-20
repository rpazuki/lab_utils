import numpy as np
# import arviz as az
import matplotlib.pyplot as plt


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


# ─────────────────────────────────────────────────────────────────────────────
# Helper samplers for testing and diagnostics
# ─────────────────────────────────────────────────────────────────────────────

def metropolis_normal(
    n: int = 10_000,
    mu: float = 0.0,
    sigma: float = 1.0,
    proposal_sd: float = 0.5,
    seed: int = 42,
) -> np.ndarray:
    """
    Metropolis-Hastings sampler targeting θ ~ N(mu, sigma²).

    Parameters
    ----------
    n : int, optional
        Number of iterations (default: 10000)
    mu : float, optional
        Mean of target distribution (default: 0.0)
    sigma : float, optional
        Standard deviation of target distribution (default: 1.0)
    proposal_sd : float, optional
        Standard deviation of proposal distribution (default: 0.5)
    seed : int, optional
        Random seed for reproducibility (default: 42)

    Returns
    -------
    np.ndarray
        MCMC chain of length n
    """
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = 0.0
    for i in range(1, n):
        prop = chain[i - 1] + rng.normal(0.0, proposal_sd)
        log_alpha = (
            -0.5 * ((prop - mu) ** 2 - (chain[i - 1] - mu) ** 2) / sigma ** 2
        )
        chain[i] = prop if np.log(rng.uniform()) < log_alpha else chain[i - 1]
    return chain


def iid_normal(
    n: int = 10_000,
    mu: float = 0.0,
    sigma: float = 1.0,
    seed: int = 1,
) -> np.ndarray:
    """
    Generate i.i.d. N(mu, sigma²) samples.

    For i.i.d. samples, effective sample size should be approximately equal to n.

    Parameters
    ----------
    n : int, optional
        Number of samples (default: 10000)
    mu : float, optional
        Mean of distribution (default: 0.0)
    sigma : float, optional
        Standard deviation of distribution (default: 1.0)
    seed : int, optional
        Random seed for reproducibility (default: 1)

    Returns
    -------
    np.ndarray
        Array of i.i.d. samples of length n
    """
    return np.random.default_rng(seed).normal(mu, sigma, size=n)


# ─────────────────────────────────────────────────────────────────────────────
# Visualization function
# ─────────────────────────────────────────────────────────────────────────────

def plot_ess_diagnostics(
    chains: dict[str, np.ndarray],
    labels: dict[str, str] | None = None,
    colors: dict[str, str] | None = None,
    max_trace_draws: int = 500,
    max_acf_lags: int = 200,
    figsize: tuple[int, int] = (15, 8),
    output_file: str | None = "ess_comparison.png",
    show_plot: bool = True,
) -> plt.Figure:
    """
    Plot trace plots and autocorrelation functions for MCMC chains.

    Creates a 2×N subplot grid where N is the number of chains. Top row shows
    trace plots (first draws), bottom row shows ACF with ESS estimates.

    Parameters
    ----------
    chains : dict[str, np.ndarray]
        Dictionary mapping chain identifiers to 1-D chain arrays.
        Example: {"iid": chain1, "good": chain2, "slow": chain3}
    labels : dict[str, str] | None
        Optional dictionary mapping chain identifiers to display labels.
        If None, uses chain identifiers as labels.
    colors : dict[str, str] | None
        Optional dictionary mapping chain identifiers to matplotlib colors.
        If None, uses default color palette.
    max_trace_draws : int, optional
        Maximum number of draws to show in trace plots (default: 500)
    max_acf_lags : int, optional
        Maximum lag to show in ACF plots (default: 200)
    figsize : tuple[int, int], optional
        Figure size as (width, height) (default: (15, 8))
    output_file : str | None, optional
        Path to save figure. If None, does not save (default: "ess_comparison.png")
    show_plot : bool, optional
        Whether to display the plot (default: True)

    Returns
    -------
    plt.Figure
        The created figure object
    """
    n_chains = len(chains)
    
    # Set defaults for labels and colors
    if labels is None:
        labels = {key: key for key in chains.keys()}
    if colors is None:
        color_palette = ["steelblue", "seagreen", "tomato", "gold", "purple", "coral"]
        colors = {key: color_palette[i % len(color_palette)] for i, key in enumerate(chains.keys())}
    
    fig, axes = plt.subplots(2, n_chains, figsize=figsize)
    
    # Handle single chain case (axes is 1D)
    if n_chains == 1:
        axes = axes.reshape(2, 1)
    
    for col, (chain_id, chain) in enumerate(chains.items()):
        label = labels.get(chain_id, chain_id)
        color = colors.get(chain_id, "steelblue")
        n = len(chain)
        
        # Top row — trace plot (first max_trace_draws draws)
        ax_tr = axes[0, col]
        trace_end = min(max_trace_draws, n)
        ax_tr.plot(chain[:trace_end], lw=0.7, color=color, alpha=0.85)
        ax_tr.axhline(0, color="gray", lw=0.5, ls="--")
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

