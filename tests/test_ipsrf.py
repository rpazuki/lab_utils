"""Unit tests for ipsrf module."""

import numpy as np
import pytest
from scipy.stats import t as student_t
from labUtils.samplings.diagnostics.ipsrf import (
    ipsrf,
    ipsrf_all_params,
    cipsrf,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helper samplers (moved from main module)
# ─────────────────────────────────────────────────────────────────────────────

def metropolis_normal(
    n: int = 5_000,
    start: float = 0.0,
    target_mu: float = 0.0,
    target_sigma: float = 1.0,
    proposal_sd: float = 0.5,
    seed: int = 0,
) -> np.ndarray:
    """
    Metropolis-Hastings sampler targeting θ ~ N(target_mu, target_sigma²).

    Parameters
    ----------
    n : int, optional
        Chain length (default: 5000)
    start : float, optional
        Initial θ value (default: 0.0)
    target_mu : float, optional
        Posterior mean (default: 0.0)
    target_sigma : float, optional
        Posterior standard deviation (default: 1.0)
    proposal_sd : float, optional
        Gaussian proposal standard deviation (default: 0.5)
    seed : int, optional
        Random seed (default: 0)

    Returns
    -------
    np.ndarray
        Chain of length n
    """
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = start
    for i in range(1, n):
        proposal = chain[i - 1] + rng.normal(0.0, proposal_sd)
        log_alpha = (
            -0.5
            * ((proposal - target_mu) ** 2 - (chain[i - 1] - target_mu) ** 2)
            / target_sigma ** 2
        )
        chain[i] = (
            proposal if np.log(rng.uniform()) < log_alpha else chain[i - 1]
        )
    return chain


def metropolis_student_t(
    n: int = 5_000,
    start: float = 0.0,
    df: float = 3.0,
    proposal_sd: float = 1.0,
    seed: int = 0,
) -> np.ndarray:
    """
    Metropolis-Hastings sampler targeting a Student-t(df) distribution.

    Parameters
    ----------
    n : int, optional
        Chain length (default: 5000)
    start : float, optional
        Initial θ value (default: 0.0)
    df : float, optional
        Degrees of freedom (default: 3.0)
    proposal_sd : float, optional
        Gaussian proposal standard deviation (default: 1.0)
    seed : int, optional
        Random seed (default: 0)

    Returns
    -------
    np.ndarray
        Chain of length n
    """
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = start
    log_target = student_t(df=df).logpdf

    for i in range(1, n):
        proposal = chain[i - 1] + rng.normal(0.0, proposal_sd)
        log_alpha = log_target(proposal) - log_target(chain[i - 1])
        chain[i] = (
            proposal if np.log(rng.uniform()) < log_alpha else chain[i - 1]
        )
    return chain


# ─────────────────────────────────────────────────────────────────────────────
# Test fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def converged_chains():
    """Scenario 1: Converged normal chains."""
    n_draws = 5_000
    n_chains = 4
    burn_in = 500
    
    rng_init = np.random.default_rng(99)
    start_points = rng_init.uniform(-2, 2, size=n_chains)
    
    chains = [
        metropolis_normal(
            n=n_draws + burn_in,
            start=float(start_points[c]),
            target_mu=0.0,
            target_sigma=1.0,
            proposal_sd=0.5,
            seed=c,
        )[burn_in:]
        for c in range(n_chains)
    ]
    return chains


@pytest.fixture
def non_converged_chains():
    """Scenario 2: Non-converged chains with diverse starting points."""
    n_chains = 4
    start_points = np.array([-8.0, -4.0, +4.0, +8.0])
    
    chains = [
        metropolis_normal(
            n=200,
            start=float(start_points[c]),
            target_mu=0.0,
            target_sigma=1.0,
            proposal_sd=0.5,
            seed=c + 100,
        )
        for c in range(n_chains)
    ]
    return chains


@pytest.fixture
def heavy_tailed_chains():
    """Scenario 3: Converged Student-t chains."""
    n_draws = 5_000
    n_chains = 4
    burn_in = 500
    
    rng_init = np.random.default_rng(99)
    
    chains = [
        metropolis_student_t(
            n=n_draws + burn_in,
            start=float(rng_init.uniform(-1, 1)),
            df=3.0,
            proposal_sd=1.5,
            seed=c + 200,
        )[burn_in:]
        for c in range(n_chains)
    ]
    return chains


@pytest.fixture
def multi_param_chains():
    """Scenario 4: Multi-parameter chains."""
    n_draws = 5_000
    n_chains = 4
    n_params = 3
    burn_in = 500
    
    chains_multi = np.stack([
        np.column_stack([
            metropolis_normal(
                n=n_draws + burn_in,
                start=0.0,
                target_mu=0.0,
                target_sigma=1.0,
                proposal_sd=0.5,
                seed=c,
            )[burn_in:],
            metropolis_normal(
                n=n_draws + burn_in,
                start=1.0,
                target_mu=2.0,
                target_sigma=0.5,
                proposal_sd=0.3,
                seed=c + 50,
            )[burn_in:],
            metropolis_student_t(
                n=n_draws + burn_in,
                start=0.0,
                df=5.0,
                proposal_sd=1.0,
                seed=c + 100,
            )[burn_in:],
        ])
        for c in range(n_chains)
    ])
    return chains_multi


# ─────────────────────────────────────────────────────────────────────────────
# Tests: IPSRF convergence (Scenario 1)
# ─────────────────────────────────────────────────────────────────────────────

class TestIPSRFConverged:
    """Tests for converged chains (section 9.3)."""

    def test_converged_chains_ipsrf_near_one(self, converged_chains):
        """Well-converged chains should have IPSRF ≈ 1."""
        r = ipsrf(converged_chains, alpha=0.05)
        assert r < 1.10, f"IPSRF {r} should be < 1.10 for converged chains"

    def test_converged_chains_ipsrf_positive(self, converged_chains):
        """IPSRF must be positive."""
        r = ipsrf(converged_chains, alpha=0.05)
        assert r > 0, "IPSRF must be positive"

    def test_converged_chains_ipsrf_finite(self, converged_chains):
        """IPSRF must be finite."""
        r = ipsrf(converged_chains, alpha=0.05)
        assert np.isfinite(r), "IPSRF must be finite"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: IPSRF non-convergence (Scenario 2)
# ─────────────────────────────────────────────────────────────────────────────

class TestIPSRFNonConverged:
    """Tests for non-converged chains (section 9.4)."""

    def test_non_converged_chains_ipsrf_large(self, non_converged_chains):
        """Non-converged chains should have IPSRF noticeably above 1."""
        r = ipsrf(non_converged_chains, alpha=0.05)
        assert r > 1.08, f"IPSRF {r} should be > 1.08 for non-converged chains"

    def test_non_converged_chains_worse_than_converged(
        self, converged_chains, non_converged_chains
    ):
        """Non-converged chains should have higher IPSRF than converged."""
        r_conv = ipsrf(converged_chains, alpha=0.05)
        r_non = ipsrf(non_converged_chains, alpha=0.05)
        assert r_non > r_conv, "Non-converged should have higher IPSRF"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: Heavy-tailed target (Scenario 3)
# ─────────────────────────────────────────────────────────────────────────────

class TestIPSRFHeavyTailed:
    """Tests for Student-t chains (section 9.5)."""

    def test_heavy_tailed_chains_ipsrf_reasonable(self, heavy_tailed_chains):
        """Converged heavy-tailed chains should have IPSRF < 1.3."""
        r = ipsrf(heavy_tailed_chains, alpha=0.05)
        assert r < 1.3, f"IPSRF {r} should be < 1.3 for converged heavy-tailed"

    def test_heavy_tailed_chains_ipsrf_still_finite(self, heavy_tailed_chains):
        """IPSRF should handle heavy tails gracefully."""
        r = ipsrf(heavy_tailed_chains, alpha=0.05)
        assert np.isfinite(r), "IPSRF must be finite even for heavy-tailed data"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: Multi-parameter IPSRF (Scenario 4)
# ─────────────────────────────────────────────────────────────────────────────

class TestIPSRFMultiParam:
    """Tests for multi-parameter IPSRF (section 9.7)."""

    def test_multi_param_ipsrf_shape(self, multi_param_chains):
        """ipsrf_all_params should return one value per parameter."""
        ipsrf_vals = ipsrf_all_params(multi_param_chains, alpha=0.05)
        n_params = multi_param_chains.shape[2]
        assert ipsrf_vals.shape == (n_params,), \
            f"Shape {ipsrf_vals.shape} != ({n_params},)"

    def test_multi_param_ipsrf_all_positive(self, multi_param_chains):
        """All IPSRF values must be positive."""
        ipsrf_vals = ipsrf_all_params(multi_param_chains, alpha=0.05)
        assert np.all(ipsrf_vals > 0), "All IPSRF values must be positive"

    def test_multi_param_ipsrf_all_finite(self, multi_param_chains):
        """All IPSRF values must be finite."""
        ipsrf_vals = ipsrf_all_params(multi_param_chains, alpha=0.05)
        assert np.all(np.isfinite(ipsrf_vals)), "All IPSRF values must be finite"

    def test_multi_param_ipsrf_reasonably_converged(self, multi_param_chains):
        """Most multi-param chains should converge reasonably."""
        ipsrf_vals = ipsrf_all_params(multi_param_chains, alpha=0.05)
        n_acceptable = np.sum(ipsrf_vals < 1.3)
        assert n_acceptable >= 2, \
            "At least 2 of 3 params should have IPSRF < 1.3"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: Cumulative IPSRF (section 9.6)
# ─────────────────────────────────────────────────────────────────────────────

class TestCIPSRF:
    """Tests for cumulative IPSRF (section 9.6)."""

    def test_cipsrf_shapes(self, converged_chains):
        """CIPSRF should return matching fractions and values."""
        fracs, vals = cipsrf(converged_chains, alpha=0.05, steps=10)
        assert fracs.shape == vals.shape, "Fractions and values must have same shape"
        assert len(fracs) == 10, "Should have 10 steps"

    def test_cipsrf_fractions_increasing(self, converged_chains):
        """Fractions should be monotonically increasing."""
        fracs, _ = cipsrf(converged_chains, alpha=0.05, steps=10)
        diffs = np.diff(fracs)
        assert np.all(diffs > 0), "Fractions should be monotonically increasing"

    def test_cipsrf_values_generally_decreasing(self, converged_chains):
        """For converged chains, IPSRF should generally decrease."""
        _, vals = cipsrf(converged_chains, alpha=0.05, steps=30)
        # Check that final IPSRF is typically lower than initial
        assert vals[-1] < vals[0], "CIPSRF should decrease as chain grows"

    def test_cipsrf_custom_min_frac(self, converged_chains):
        """Test with custom min_frac."""
        fracs, vals = cipsrf(converged_chains, alpha=0.05, min_frac=0.2, steps=15)
        assert len(fracs) == 15
        assert fracs[0] >= 0.2, "First fraction should be >= min_frac"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: Edge cases and alpha parameter
# ─────────────────────────────────────────────────────────────────────────────

class TestIPSRFAlpha:
    """Tests for different alpha values."""

    def test_ipsrf_alpha_0_01(self, converged_chains):
        """Test with alpha=0.01 (99% interval)."""
        r = ipsrf(converged_chains, alpha=0.01)
        assert np.isfinite(r) and r > 0

    def test_ipsrf_alpha_0_20(self, converged_chains):
        """Test with alpha=0.20 (80% interval)."""
        r = ipsrf(converged_chains, alpha=0.20)
        assert np.isfinite(r) and r > 0

    def test_ipsrf_alpha_effects_on_value(self, converged_chains):
        """Different alphas should give different IPSRF values."""
        r_01 = ipsrf(converged_chains, alpha=0.01)
        r_05 = ipsrf(converged_chains, alpha=0.05)
        assert r_01 != r_05, "Different alphas should yield different values"


class TestIPSRFEdgeCases:
    """Tests for edge cases."""

    def test_ipsrf_empty_quantile_handling(self):
        """Test that quantile edge case doesn't crash."""
        # Small chains to test edge cases
        chains = [
            np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
            np.array([1.5, 2.5, 3.5, 4.5, 5.5]),
        ]
        r = ipsrf(chains, alpha=0.05)
        assert np.isfinite(r)

    def test_ipsrf_single_chain(self):
        """IPSRF should handle single chain."""
        chain = np.random.randn(1000)
        r = ipsrf([chain], alpha=0.05)
        assert np.isfinite(r)
