import numpy as np
import pytest

from labUtils.samplings.diagnostics import geweke_test

pytestmark = pytest.mark.skipif(
    geweke_test.yule_walker is None,
    reason="statsmodels is required for AR-based Geweke example expectations.",
)


def _metropolis_normal(n=5000, mu=0.0, sigma=1.0, proposal_sd=0.5, seed=42):
    """Simple Metropolis-Hastings sampler for theta ~ N(mu, sigma^2)."""
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = 0.0
    for i in range(1, n):
        proposal = chain[i - 1] + rng.normal(0, proposal_sd)
        log_alpha = -0.5 * ((proposal - mu) ** 2 - (chain[i - 1] - mu) ** 2) / sigma ** 2
        chain[i] = proposal if np.log(rng.uniform()) < log_alpha else chain[i - 1]
    return chain


def _slow_drift_chain(n=5000, start=8.0, target=0.0, phi=0.998, noise=0.1, seed=7):
    """AR(1) chain initialized far from stationarity with slow drift."""
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = start
    for i in range(1, n):
        chain[i] = phi * chain[i - 1] + (1 - phi) * target + rng.normal(0, noise)
    return chain


def test_geweke_example1_converged_chain():
    """Example 1: a converged MH chain should not reject stationarity."""
    chain = _metropolis_normal(n=5000)

    z_single = geweke_test.geweke_zscore(chain)
    assert abs(z_single) < 1.96

    starts, zscores = geweke_test.geweke_sweep(chain, intervals=20)
    assert len(starts) == 20
    assert len(zscores) == 20


def test_geweke_example2_non_converged_chain():
    """Example 2: a slow-drift chain should show large Geweke Z values."""
    chain = _slow_drift_chain(n=2000)

    starts, zscores = geweke_test.geweke_sweep(chain, intervals=20)
    assert len(starts) == 20
    assert len(zscores) == 20
    assert np.max(np.abs(zscores)) > 1.96
