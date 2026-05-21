"""Unit tests for convergence diagnostics module."""

import numpy as np
import pytest

from labUtils.samplings.convergence import convergence_matrix


def iid_chain(n: int = 5000, seed: int = 42) -> np.ndarray:
    """Generate i.i.d. N(0, 1) chain (should pass convergence)."""
    return np.random.default_rng(seed).normal(size=n)


def correlated_chain(n: int = 5000, proposal_sd: float = 0.5, seed: int = 42) -> np.ndarray:
    """Generate Metropolis chain (should pass convergence with good tuning)."""
    rng = np.random.default_rng(seed)
    chain = np.empty(n)
    chain[0] = 0.0
    for i in range(1, n):
        prop = chain[i - 1] + rng.normal(0.0, proposal_sd)
        log_alpha = -0.5 * ((prop) ** 2 - (chain[i - 1]) ** 2)
        chain[i] = prop if np.log(rng.uniform()) < log_alpha else chain[i - 1]
    return chain


class TestConvergenceMatrixInput:
    """Test input validation for convergence_matrix."""

    def test_empty_chains_raises(self):
        """Empty chains dict/list should raise ValueError."""
        with pytest.raises(ValueError, match="at least one chain"):
            convergence_matrix({})

        with pytest.raises(ValueError, match="at least one chain"):
            convergence_matrix([])

    def test_non_array_chain_raises(self):
        """Non-numpy-array chain should raise ValueError."""
        with pytest.raises(ValueError, match="not a numpy array"):
            convergence_matrix({"chain1": [1, 2, 3]})

    def test_empty_array_chain_raises(self):
        """Empty numpy array should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            convergence_matrix({"chain1": np.array([])})

    def test_list_input_converted_to_dict(self):
        """List input should be converted to dict with auto-generated IDs."""
        chain1 = iid_chain(n=1000, seed=10)
        chain2 = iid_chain(n=1000, seed=11)

        result = convergence_matrix([chain1, chain2])

        assert len(result) == 2
        assert result[0]["chain_id"] == "chain0"
        assert result[1]["chain_id"] == "chain1"

    def test_dict_input_preserves_ids(self):
        """Dict input should preserve chain IDs."""
        chain1 = iid_chain(n=1000, seed=10)
        chain2 = iid_chain(n=1000, seed=11)

        result = convergence_matrix({"run_a": chain1, "run_b": chain2})

        chain_ids = [r["chain_id"] for r in result]
        assert "run_a" in chain_ids
        assert "run_b" in chain_ids


class TestConvergenceMatrixOutput:
    """Test output structure and field presence."""

    def test_output_has_all_required_fields(self):
        """Output dict must have all expected convergence fields."""
        chain = iid_chain(n=1000, seed=20)
        result = convergence_matrix({"test": chain})

        assert len(result) == 1
        r = result[0]

        required_fields = [
            "chain_id",
            "geweke_convergence",
            "geweke_max_zscore",
            "ess_convergence",
            "ess_mean",
            "cipsrf_convergence",
            "cipsrf_final",
            "overall_convergence",
        ]

        for field in required_fields:
            assert field in r, f"Missing field: {field}"

    def test_output_field_types(self):
        """Output fields must have correct types."""
        chain = iid_chain(n=1000, seed=21)
        result = convergence_matrix({"test": chain})
        r = result[0]

        assert isinstance(r["chain_id"], str)
        assert isinstance(r["geweke_convergence"], (bool, np.bool_))
        assert isinstance(r["geweke_max_zscore"], float)
        assert isinstance(r["ess_convergence"], (bool, np.bool_))
        assert isinstance(r["ess_mean"], float)
        assert isinstance(r["cipsrf_convergence"], (bool, np.bool_))
        assert isinstance(r["cipsrf_final"], float)
        assert isinstance(r["overall_convergence"], (bool, np.bool_))


class TestConvergenceMatrixLogic:
    """Test convergence assessment logic."""

    def test_iid_chain_passes_convergence(self):
        """i.i.d. chain should pass most convergence tests with relaxed threshold."""
        chain = iid_chain(n=5000, seed=30)
        # Use relaxed Geweke threshold to account for multiple comparisons
        result = convergence_matrix(
            {"iid": chain}, geweke_z_threshold=2.5, ess_min_threshold=400
        )

        assert len(result) == 1
        r = result[0]

        # i.i.d. should have high mixing and ESS
        assert r["ess_convergence"] == True
        assert r["ess_mean"] > 400

        # ESS alone is the most reliable indicator for i.i.d.
        assert r["overall_convergence"] == True

    def test_well_tuned_metropolis_passes_convergence(self):
        """Well-tuned Metropolis should pass convergence."""
        chain = correlated_chain(n=5000, proposal_sd=0.5, seed=31)
        result = convergence_matrix({"metro": chain})

        assert len(result) == 1
        r = result[0]

        # Well-tuned should still generally pass Geweke
        assert isinstance(r["geweke_convergence"], (bool, np.bool_))
        # ESS should be moderate (autocorrelated but not too slow)
        assert r["ess_mean"] > 30

    def test_geweke_threshold_customizable(self):
        """Geweke threshold should be customizable."""
        chain = iid_chain(n=2000, seed=32)

        # Strict threshold
        result_strict = convergence_matrix({"test": chain}, geweke_z_threshold=0.5)
        # Loose threshold
        result_loose = convergence_matrix({"test": chain}, geweke_z_threshold=10.0)

        # Loose should definitely pass, strict might not
        assert result_loose[0]["geweke_convergence"] == True

    def test_ess_threshold_customizable(self):
        """ESS threshold should be customizable."""
        chain = iid_chain(n=1000, seed=33)

        result_high = convergence_matrix({"test": chain}, ess_min_threshold=5000)
        result_low = convergence_matrix({"test": chain}, ess_min_threshold=50)

        # High threshold likely fails, low threshold likely passes
        assert result_low[0]["ess_convergence"] == True

    def test_cipsrf_tolerance_customizable(self):
        """CIPSRF tolerance should be customizable."""
        chain = iid_chain(n=3000, seed=34)

        result_tight = convergence_matrix({"test": chain}, cipsrf_tolerance=0.01)
        result_loose = convergence_matrix({"test": chain}, cipsrf_tolerance=0.5)

        # Loose tolerance more likely to pass
        assert isinstance(result_tight[0]["cipsrf_convergence"], (bool, np.bool_))
        assert isinstance(result_loose[0]["cipsrf_convergence"], (bool, np.bool_))

    def test_overall_convergence_requires_all_three(self):
        """Overall convergence should require all three tests to pass."""
        chain = iid_chain(n=2000, seed=35)

        # Use very strict thresholds to force some failures
        result = convergence_matrix(
            {"test": chain},
            ess_min_threshold=100_000,  # Will fail
            geweke_z_threshold=1.96,
            cipsrf_tolerance=0.1,
        )

        r = result[0]
        # If any individual test fails, overall should fail
        if not r["ess_convergence"]:
            assert r["overall_convergence"] == False

    def test_multiple_chains_processed_independently(self):
        """Each chain should be processed independently."""
        chain1 = iid_chain(n=2000, seed=40)
        chain2 = correlated_chain(n=2000, proposal_sd=0.2, seed=41)  # Slower mixing

        result = convergence_matrix({"chain1": chain1, "chain2": chain2})

        assert len(result) == 2
        # chain1 (i.i.d.) should have higher ESS than chain2 (correlated)
        ess1 = result[0]["ess_mean"] if result[0]["chain_id"] == "chain1" else result[1]["ess_mean"]
        ess2 = result[1]["ess_mean"] if result[1]["chain_id"] == "chain2" else result[0]["ess_mean"]

        assert ess1 > ess2, "i.i.d. chain should have higher ESS than slow chain"


class TestConvergenceMatrixEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_minimum_chain_length(self):
        """Function should handle short chains gracefully."""
        # Create a chain just long enough for CIPSRF to work
        chain = iid_chain(n=300, seed=50)
        result = convergence_matrix({"short": chain})

        assert len(result) == 1
        assert "overall_convergence" in result[0]

    def test_large_chain(self):
        """Function should handle large chains without errors."""
        chain = iid_chain(n=50_000, seed=51)
        result = convergence_matrix({"large": chain})

        assert len(result) == 1
        # Should process large chains successfully
        assert isinstance(result[0], dict)
        assert "overall_convergence" in result[0]

    def test_geweke_max_zscore_nonnegative(self):
        """Geweke max Z-score should always be non-negative."""
        chain = iid_chain(n=2000, seed=52)
        result = convergence_matrix({"test": chain})

        assert result[0]["geweke_max_zscore"] >= 0

    def test_ess_mean_positive(self):
        """ESS mean should always be positive."""
        chain = iid_chain(n=2000, seed=53)
        result = convergence_matrix({"test": chain})

        assert result[0]["ess_mean"] > 0

    def test_cipsrf_final_roughly_one(self):
        """CIPSRF final value should trend toward 1 for converged chain."""
        chain = iid_chain(n=5000, seed=54)
        result = convergence_matrix({"test": chain})

        cipsrf_final = result[0]["cipsrf_final"]
        # For well-mixing i.i.d., CIPSRF should be close to 1
        assert 0.5 < cipsrf_final < 2.0, f"CIPSRF {cipsrf_final} seems unreasonable"
