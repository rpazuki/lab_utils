"""Unit tests for effective_sample_size module."""

import numpy as np
import pytest
# import arviz as az
from labUtils.samplings.diagnostics.effective_sample_size import (
    autocorr_fft,
    geyer_ips,
    ess_from_chain,
    ess_sweep,
)


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

class TestSamplers:
    """Tests for sampler helper functions (section 7.2)."""

    def test_metropolis_normal_shape(self):
        """Test that Metropolis-Hastings sampler returns correct shape."""
        n = 5000
        chain = metropolis_normal(n=n, seed=42)
        assert chain.shape == (n,)
        assert isinstance(chain, np.ndarray)

    def test_metropolis_normal_mean(self):
        """Test that Metropolis-Hastings chain has reasonable mean."""
        chain = metropolis_normal(n=10_000, mu=0.0, sigma=1.0, seed=42)
        assert abs(np.mean(chain)) < 0.1  # rough check

    def test_metropolis_normal_deterministic(self):
        """Test that sampler is deterministic with fixed seed."""
        chain1 = metropolis_normal(n=1000, seed=42)
        chain2 = metropolis_normal(n=1000, seed=42)
        np.testing.assert_array_equal(chain1, chain2)

    def test_metropolis_normal_diff_seeds(self):
        """Test that different seeds produce different chains."""
        chain1 = metropolis_normal(n=1000, seed=42)
        chain2 = metropolis_normal(n=1000, seed=43)
        assert not np.allclose(chain1, chain2)

    def test_iid_normal_shape(self):
        """Test that iid sampler returns correct shape."""
        n = 5000
        samples = iid_normal(n=n, seed=1)
        assert samples.shape == (n,)
        assert isinstance(samples, np.ndarray)

    def test_iid_normal_mean(self):
        """Test that iid samples have correct empirical mean."""
        samples = iid_normal(n=50_000, mu=0.0, sigma=1.0, seed=1)
        assert abs(np.mean(samples)) < 0.05

    def test_iid_normal_deterministic(self):
        """Test that iid sampler is deterministic with fixed seed."""
        samples1 = iid_normal(n=1000, seed=1)
        samples2 = iid_normal(n=1000, seed=1)
        np.testing.assert_array_equal(samples1, samples2)


class TestESSComparison:
    """Tests comparing three ESS scenarios (section 7.3)."""

    @pytest.fixture
    def chains(self):
        """Generate test chains for comparison."""
        n = 10_000
        chain_iid = iid_normal(n=n, seed=1)
        chain_good = metropolis_normal(n=n, proposal_sd=0.5, seed=42)
        chain_slow = metropolis_normal(n=n, proposal_sd=0.05, seed=42)
        return {
            "iid": chain_iid,
            "good": chain_good,
            "slow": chain_slow,
        }

    def test_iid_ess_near_n(self, chains):
        """i.i.d. chain should have ESS ≈ N."""
        ess = ess_from_chain(chains["iid"])
        n = len(chains["iid"])
        # For i.i.d., ESS should be close to N (within ~10%)
        assert ess > 0.9 * n, f"ESS {ess} not close to N={n}"

    def test_well_tuned_ess_reduced(self, chains):
        """Well-tuned MH should have reduced but reasonable ESS."""
        ess = ess_from_chain(chains["good"])
        n = len(chains["good"])
        # For well-tuned MH, ESS should be reduced but not tiny (typically 5-20% of n)
        assert 0.02 * n < ess < 0.5 * n, f"ESS {ess} not in expected range"

    def test_slow_mixing_ess_much_lower(self, chains):
        """Slow-mixing chain should have much lower ESS."""
        ess = ess_from_chain(chains["slow"])
        n = len(chains["slow"])
        # Slow mixing should give very low ESS
        assert ess < 0.2 * n, f"ESS {ess} not low enough for slow mixing"

    def test_iid_vs_good_ess_ratio(self, chains):
        """i.i.d. chain should have higher ESS than well-tuned MH."""
        ess_iid = ess_from_chain(chains["iid"])
        ess_good = ess_from_chain(chains["good"])
        assert ess_iid > ess_good, "i.i.d. should have higher ESS than well-tuned MH"

    def test_good_vs_slow_ess_ratio(self, chains):
        """Well-tuned MH should have higher ESS than slow-mixing."""
        ess_good = ess_from_chain(chains["good"])
        ess_slow = ess_from_chain(chains["slow"])
        assert ess_good > ess_slow, "Well-tuned should have higher ESS than slow"

    # def test_ess_matches_arviz_iid(self, chains):
    #     """Manual ESS should be reasonably close to ArviZ ESS for i.i.d."""
    #     chain = chains["iid"]
    #     ess_manual = ess_from_chain(chain)
    #     ess_arviz = float(az.ess({"theta": chain[np.newaxis, :]})["theta"].values)
    #     # Should be within 10% of each other
    #     ratio = ess_manual / ess_arviz
    #     assert 0.9 < ratio < 1.1, f"ESS ratio {ratio} too far from 1.0"

    # def test_ess_matches_arviz_good_mh(self, chains):
    #     """Manual ESS should be reasonably close to ArviZ for well-tuned MH."""
    #     chain = chains["good"]
    #     ess_manual = ess_from_chain(chain)
    #     ess_arviz = float(az.ess({"theta": chain[np.newaxis, :]})["theta"].values)
    #     # Should be within 20% of each other (more tolerance for MH)
    #     ratio = ess_manual / ess_arviz
    #     assert 0.8 < ratio < 1.2, f"ESS ratio {ratio} too far from 1.0"


class TestAutocorr:
    """Tests for autocorrelation function."""

    def test_autocorr_fft_starts_with_one(self):
        """Autocorr at lag 0 should be 1.0."""
        x = np.random.default_rng(101).normal(size=1000)
        rho = autocorr_fft(x)
        assert rho[0] == pytest.approx(1.0, abs=1e-10)

    def test_autocorr_fft_decreases(self):
        """Autocorr should generally decrease with lag."""
        chain = metropolis_normal(n=5000, seed=42)
        rho = autocorr_fft(chain)
        # For positive correlation, should generally decrease
        assert rho[1] <= rho[0]
        assert rho[10] <= rho[1]

    def test_autocorr_fft_iid_small(self):
        """i.i.d. data should have small autocorr at lag > 0."""
        samples = np.random.default_rng(102).normal(size=10_000)
        rho = autocorr_fft(samples)
        # For i.i.d., rho[k] should be small for k > 0
        assert abs(rho[1]) < 0.1


class TestGeyer:
    """Tests for Geyer's IPS estimator."""

    def test_geyer_ips_ge_one(self):
        """τ_int should always be ≥ 1."""
        rng = np.random.default_rng(103)
        for _ in range(10):
            x = rng.normal(size=1000)
            rho = autocorr_fft(x)
            tau = geyer_ips(rho)
            assert tau >= 1.0

    def test_geyer_ips_iid_near_one(self):
        """For i.i.d., τ_int should be close to 1."""
        samples = np.random.default_rng(104).normal(size=50_000)
        rho = autocorr_fft(samples)
        tau = geyer_ips(rho)
        assert 0.8 < tau < 1.5, f"τ_int={tau} not close to 1 for i.i.d."

    def test_geyer_ips_correlated_large(self):
        """For correlated data, τ_int should be larger."""
        chain = metropolis_normal(n=10_000, proposal_sd=0.05, seed=42)
        rho = autocorr_fft(chain)
        tau = geyer_ips(rho)
        assert tau > 5, f"τ_int={tau} not large for correlated chain"


class TestESSSweep:
    """Tests for ESS over fixed-length subchains."""

    def test_ess_sweep_matches_manual_chunk_computation(self):
        """ess_sweep should match manual ESS computed on each subchain."""
        chain = iid_normal(n=1000, seed=11)
        subchains_length = 200

        ess_values = ess_sweep(chain, subchains_length=subchains_length)
        expected = np.array(
            [
                ess_from_chain(chain[i * 200 : (i + 1) * 200])
                for i in range(len(chain) // 200)
            ]
        )

        np.testing.assert_allclose(ess_values, expected)

    def test_ess_sweep_fraction_overrides_subchains_length(self):
        """fraction should set subchain length and override explicit length."""
        chain = iid_normal(n=1000, seed=22)

        by_fraction = ess_sweep(chain, subchains_length=50, fraction=0.2)
        by_length = ess_sweep(chain, subchains_length=200)

        np.testing.assert_allclose(by_fraction, by_length)

    def test_ess_sweep_uses_single_subchain_when_length_exceeds_chain(self):
        """If requested subchain length exceeds chain, use full chain once."""
        chain = metropolis_normal(n=300, seed=33)
        ess_values = ess_sweep(chain, subchains_length=1000)

        assert ess_values.shape == (1,)
        assert ess_values[0] == pytest.approx(ess_from_chain(chain))

    def test_ess_sweep_invalid_subchains_length_raises(self):
        """subchains_length must be positive."""
        chain = iid_normal(n=100, seed=44)

        with pytest.raises(ValueError, match="subchains_length"):
            ess_sweep(chain, subchains_length=0)

    def test_ess_sweep_invalid_fraction_raises(self):
        """fraction must lie in (0, 1]."""
        chain = iid_normal(n=100, seed=55)

        with pytest.raises(ValueError, match="fraction"):
            ess_sweep(chain, fraction=0.0)

        with pytest.raises(ValueError, match="fraction"):
            ess_sweep(chain, fraction=1.2)

    def test_ess_sweep_empty_chain_raises(self):
        """Empty chain should raise ValueError."""
        chain = np.array([])

        with pytest.raises(ValueError, match="chain"):
            ess_sweep(chain)
