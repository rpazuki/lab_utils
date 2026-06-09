import numpy as np
from statsmodels.regression.linear_model import yule_walker
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# Core functions
# ─────────────────────────────────────────────────────────────────────────────

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


# ─────────────────────────────────────────────────────────────────────────────
# Example 1 — Converged chain (Metropolis-Hastings targeting N(0, 1))
# ─────────────────────────────────────────────────────────────────────────────

def metropolis_normal(n=5000, mu=0.0, sigma=1.0, proposal_sd=0.5, seed=42):
    """Simple Metropolis-Hastings sampler for θ ~ N(mu, sigma²)."""
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = 0.0
    for i in range(1, n):
        proposal = chain[i - 1] + rng.normal(0, proposal_sd)
        log_alpha = -0.5 * ((proposal - mu) ** 2 - (chain[i - 1] - mu) ** 2) / sigma ** 2
        chain[i] = proposal if np.log(rng.uniform()) < log_alpha else chain[i - 1]
    return chain


chain_good = metropolis_normal(n=5000)

# Single Z-score for the whole chain
z_single = geweke_zscore(chain_good)
print(f"Single Geweke Z-score: {z_single:.3f}")
print(f"  → |Z| < 1.96: {'CONVERGED (fail to reject H₀)' if abs(z_single) < 1.96 else 'NON-CONVERGED (reject H₀)'}")

# Multi-interval sweep
starts_good, zscores_good = geweke_sweep(chain_good, intervals=20)


# ─────────────────────────────────────────────────────────────────────────────
# Example 2 — Non-converged chain (slow AR(1) drift from far starting point)
# ─────────────────────────────────────────────────────────────────────────────

def slow_drift_chain(n=5000, start=8.0, target=0.0, phi=0.998, noise=0.1, seed=7):
    """
    AR(1) chain initialised far from stationarity.
    Slowly drifts toward target_mu but has not reached it within n steps.
    """
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = start
    for i in range(1, n):
        chain[i] = phi * chain[i - 1] + (1 - phi) * target + rng.normal(0, noise)
    return chain


chain_bad = slow_drift_chain(n=5000)

starts_bad, zscores_bad = geweke_sweep(chain_bad, intervals=20)
print(f"\nSlow-drift chain — max |Z|: {np.max(np.abs(zscores_bad)):.2f}")
print(f"  → Large |Z| values indicate the chain has NOT converged")


# ─────────────────────────────────────────────────────────────────────────────
# Plot both sweep diagnostics side-by-side
# ─────────────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

for ax, starts, zscores, label, color in zip(
    axes,
    [starts_good, starts_bad],
    [zscores_good, zscores_bad],
    ["Converged MH chain\n(targeting N(0,1))", "Non-converged AR(1) chain\n(slow drift from start=8)"],
    ["steelblue", "tomato"],
):
    ax.axhspan(-1.96, 1.96, alpha=0.15, color="green", label="95% acceptance band")
    ax.axhline(0,     color="gray",  lw=0.8, ls="-")
    ax.axhline( 1.96, color="red",   lw=1.2, ls="--", label="±1.96 threshold")
    ax.axhline(-1.96, color="red",   lw=1.2, ls="--")
    ax.scatter(starts, zscores, color=color, s=60, zorder=4)
    ax.set_xlabel("Start iteration", fontsize=11)
    ax.set_ylabel("Geweke Z-score", fontsize=11)
    ax.set_title(label, fontsize=11)
    ax.legend(fontsize=9)

plt.suptitle("Geweke Convergence Diagnostic — Multi-Interval Sweep", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig("geweke_sweep_comparison.png", dpi=150, bbox_inches="tight")
plt.show()