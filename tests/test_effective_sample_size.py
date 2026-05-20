"""Unit tests for effective_sample_size module."""

import numpy as np
import pytest
# import arviz as az
from labUtils.samplings.diagnostics.effective_sample_size import (
    autocorr_fft,
    geyer_ips,
    ess_from_chain,
    metropolis_normal,
    iid_normal,
)


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
        x = np.random.randn(1000)
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
        samples = np.random.randn(10_000)
        rho = autocorr_fft(samples)
        # For i.i.d., rho[k] should be small for k > 0
        assert abs(rho[1]) < 0.1


class TestGeyer:
    """Tests for Geyer's IPS estimator."""

    def test_geyer_ips_ge_one(self):
        """τ_int should always be ≥ 1."""
        for _ in range(10):
            x = np.random.randn(1000)
            rho = autocorr_fft(x)
            tau = geyer_ips(rho)
            assert tau >= 1.0

    def test_geyer_ips_iid_near_one(self):
        """For i.i.d., τ_int should be close to 1."""
        samples = np.random.randn(50_000)
        rho = autocorr_fft(samples)
        tau = geyer_ips(rho)
        assert 0.8 < tau < 1.5, f"τ_int={tau} not close to 1 for i.i.d."

    def test_geyer_ips_correlated_large(self):
        """For correlated data, τ_int should be larger."""
        chain = metropolis_normal(n=10_000, proposal_sd=0.05, seed=42)
        rho = autocorr_fft(chain)
        tau = geyer_ips(rho)
        assert tau > 5, f"τ_int={tau} not large for correlated chain"
