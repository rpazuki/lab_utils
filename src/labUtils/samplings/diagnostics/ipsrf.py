import numpy as np
import matplotlib.pyplot as plt
from typing import Sequence


# ─────────────────────────────────────────────────────────────────────────────
# §9.1  Core IPSRF computation
# ─────────────────────────────────────────────────────────────────────────────

def ipsrf(chains: Sequence[np.ndarray], alpha: float = 0.05) -> float:
    """
    Compute the Interval-based Potential Scale Reduction Factor (IPSRF).

    Implements the Brooks & Gelman (1998) interval-based convergence
    diagnostic. Compares the width of the central 100(1-alpha)%
    credible interval from the *pooled* chains against the *average*
    width of the same interval within individual chains.

    Convergence  →  IPSRF ≈ 1.0
    Non-convergence  →  IPSRF >> 1.0

    Parameters
    ----------
    chains : Sequence[np.ndarray]
        Sequence of 1-D arrays representing MCMC chains for a scalar parameter.
        All chains should have the same length (burn-in discarded).
    alpha : float, optional
        Tail probability; uses the central 100(1-alpha)% interval (default: 0.05).
        Default gives the 95% credible interval (2.5th–97.5th percentile).

    Returns
    -------
    float
        IPSRF value (≥ 1 in practice; target < 1.1)
    """
    chains = [np.asarray(c, dtype=float) for c in chains]
    lo, hi = alpha / 2.0, 1.0 - alpha / 2.0

    # Width for each individual chain
    within_widths = np.array([
        np.quantile(c, hi) - np.quantile(c, lo)
        for c in chains
    ])
    d_bar = within_widths.mean()

    # Width for the pooled chain
    pooled = np.concatenate(chains)
    D = np.quantile(pooled, hi) - np.quantile(pooled, lo)

    if d_bar == 0:
        return float("inf")
    return float(D / d_bar)


def ipsrf_all_params(
    chains_array: np.ndarray, alpha: float = 0.05
) -> np.ndarray:
    """
    Compute IPSRF for every parameter in a multi-chain, multi-parameter run.

    Parameters
    ----------
    chains_array : np.ndarray
        Shape (n_chains, n_draws, n_params)
    alpha : float, optional
        Tail probability (default: 0.05 → 95% interval)

    Returns
    -------
    np.ndarray
        Array of shape (n_params,) with one IPSRF value per parameter
    """
    n_chains, n_draws, n_params = chains_array.shape
    return np.array([
        ipsrf([chains_array[c, :, p] for c in range(n_chains)], alpha=alpha)
        for p in range(n_params)
    ])


def cipsrf(
    chains: Sequence[np.ndarray],
    alpha: float = 0.05,
    min_frac: float = 0.1,
    steps: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Cumulative IPSRF (CIPSRF) — tracks IPSRF as chain length grows.

    For each fraction f in linspace(min_frac, 1.0, steps), uses only
    the first floor(f * n) samples of each chain and computes IPSRF.

    Parameters
    ----------
    chains : Sequence[np.ndarray]
        Sequence of 1-D arrays of the same length n
    alpha : float, optional
        Tail probability (default: 0.05)
    min_frac : float, optional
        Minimum fraction to start at, must be in (0, 1) (default: 0.1)
    steps : int, optional
        Number of evaluation points (default: 50)

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        - fractions: Array of shape (steps,) with fraction of chain used
        - ipsrf_vals: Array of shape (steps,) with IPSRF at each fraction
    """
    chains = [np.asarray(c, dtype=float) for c in chains]
    n = min(len(c) for c in chains)
    fractions = np.linspace(min_frac, 1.0, steps)
    ipsrf_vals = np.array([
        ipsrf([c[: int(f * n)] for c in chains], alpha=alpha)
        for f in fractions
    ])
    return fractions, ipsrf_vals


# ─────────────────────────────────────────────────────────────────────────────
# Visualization function
# ─────────────────────────────────────────────────────────────────────────────

def plot_ipsrf_diagnostics(
    chains: dict[str, np.ndarray] | list[np.ndarray],
    labels: dict[str, str] | list[str] | None = None,
    colors: dict[str, str] | list[str] | None = None,
    max_trace_draws: int = 500,
    figsize: tuple[int, int] = (16, 10),
    output_file: str | None = None,
    show_plot: bool = True,
) -> plt.Figure:
    """
    Plot IPSRF convergence diagnostics for MCMC chains.

    Creates a 2×N subplot grid where N is the number of chains. Top row shows
    trace plots, bottom row shows interval comparison and IPSRF values.

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
    figsize : tuple[int, int], optional
        Figure size as (width, height) (default: (16, 10))
    output_file : str | None, optional
        Path to save figure. If None, does not save (default: None)
    show_plot : bool, optional
        Whether to display the plot (default: True)

    Returns
    -------
    plt.Figure
        The created figure object
    """
    # Convert to dict format for consistency
    if isinstance(chains, list):
        chains = {f"chain{i}": c for i, c in enumerate(chains)}
    
    n_chains = len(chains)
    chain_ids = list(chains.keys())
    chain_list = [chains[cid] for cid in chain_ids]
    
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
    
    # Handle single chain case
    if n_chains == 1:
        axes = axes.reshape(2, 1)
    
    for col, chain_id in enumerate(chain_ids):
        chain = chains[chain_id]
        label = labels.get(chain_id, chain_id)
        color = colors.get(chain_id, "steelblue")
        
        # Top row — trace plot (first max_trace_draws draws)
        ax_tr = axes[0, col]
        trace_end = min(max_trace_draws, len(chain))
        ax_tr.plot(chain[:trace_end], lw=0.8, color=color, alpha=0.85)
        ax_tr.set_title(label, fontsize=10)
        ax_tr.set_xlabel("Iteration", fontsize=9)
        ax_tr.set_ylabel(r"$\theta$", fontsize=9)
        
        # Bottom row — interval plot (95% credible interval)
        ax_int = axes[1, col]
        lo, hi = 0.025, 0.975
        q_lo = np.quantile(chain, lo)
        q_hi = np.quantile(chain, hi)
        
        ipsrf_val = ipsrf([chain], alpha=0.05)
        
        ax_int.barh(0.5, q_hi - q_lo, left=q_lo, height=0.4,
                   color=color, alpha=0.7)
        ax_int.axvline(chain.mean(), color="darkred", lw=2, ls="--", label="Mean")
        ax_int.set_yticks([0.5])
        ax_int.set_yticklabels(["95% CI"], fontsize=9)
        ax_int.set_xlabel(r"$\theta$ value", fontsize=9)
        ax_int.set_title(f"IPSRF = {ipsrf_val:.3f}", fontsize=10)
        ax_int.legend(fontsize=8)
        ax_int.set_ylim(0, 1)
    
    plt.suptitle(
        "Interval-based Potential Scale Reduction Factor (IPSRF) — Convergence Diagnostics",
        fontsize=13,
    )
    plt.tight_layout()
    
    if output_file is not None:
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
    
    if show_plot:
        plt.show()
    
    return fig

