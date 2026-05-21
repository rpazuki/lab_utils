
from typing import Any

import numpy as np

from .diagnostics.geweke_test import geweke_sweep
from .diagnostics.effective_sample_size import ess_sweep
from .diagnostics.ipsrf import cipsrf


def convergence_matrix(
    chains: dict[str, np.ndarray] | list[np.ndarray],
    ess_min_threshold: float = 400.0,
    geweke_z_threshold: float = 1.96,
    cipsrf_tolerance: float = 0.1,
    geweke_first: float = 0.1,
    geweke_last: float = 0.5,
    geweke_intervals: int = 20,
    ess_subchain_length: int = 10000,
    cipsrf_steps: int = 20
) -> list[dict[str, Any]]:
    """Assess convergence of MCMC chains using multiple diagnostics.

    Computes three independent convergence measures (Geweke, ESS, CIPSRF) for
    each chain and returns a convergence status dict for each. A chain is
    considered converged if all three diagnostics pass their respective thresholds.

    Parameters
    ----------
    chains : dict[str, np.ndarray] | list[np.ndarray]
        Dictionary mapping chain IDs to 1-D arrays, or a list of 1-D arrays.
        Example: {"chain_a": data1, "chain_b": data2} or [data1, data2]
    ess_min_threshold : float, optional
        Minimum acceptable mean ESS across subchains (default: 400.0).
    geweke_z_threshold : float, optional
        Absolute Z-score threshold for Geweke test; values below this
        indicate convergence (default: 1.96 for 5% significance).
    cipsrf_tolerance : float, optional
        Tolerance band for CIPSRF convergence; values should trend toward
        1.0 and stay within [1 - tolerance, 1 + tolerance] (default: 0.1).
    geweke_first : float, optional
        Fraction of chain to use as "first" segment in Geweke test (default: 0.1).
    geweke_last : float, optional
        Fraction of chain to use as "last" segment in Geweke test (default: 0.5).
    geweke_intervals : int, optional
        Number of intervals to evaluate in Geweke test (default: 20).
    ess_subchain_length : int, optional
        Length of subchains to compute ESS on (default: 10000).
    cipsrf_steps : int, optional
        Number of evaluation points for CIPSRF convergence tracking (default: 20).  

    Returns
    -------
    list[dict[str, Any]]
        List of result dictionaries, one per chain, each containing:
        - chain_id: str, chain identifier
        - geweke_convergence: bool, passed Geweke test
        - geweke_max_zscore: float, maximum absolute Z-score
        - ess_convergence: bool, passed ESS threshold
        - ess_mean: float, mean ESS across subchains
        - cipsrf_convergence: bool, passed CIPSRF bounds
        - cipsrf_final: float, final CIPSRF value
        - overall_convergence: bool, all three tests passed

    Raises
    ------
    ValueError
        If chains is empty, chains have inconsistent lengths, or individual
        chain is not a numpy array.
    """
    # Convert to dict format for consistency
    if isinstance(chains, list):
        chains = {f"chain{i}": c for i, c in enumerate(chains)}
        
    if not chains:
        raise ValueError("chains must contain at least one chain.")

    results = []
    for cid, chain in chains.items():
        if not isinstance(chain, np.ndarray):
            raise ValueError(f"Chain {cid} is not a numpy array.")
        if len(chain) == 0:
            raise ValueError(f"Chain {cid} is empty.")
        
        # Geweke test: check if all Z-scores are within threshold
        sweep_starts, geweke_z_scores = geweke_sweep(
            chain, first=geweke_first, last=geweke_last, intervals=geweke_intervals
        )
        max_zscore = np.max(np.abs(geweke_z_scores))
        geweke_convergence = max_zscore < geweke_z_threshold

        # Effective sample size: check mean ESS across subchains
        ess_sweep_values = ess_sweep(chain, subchains_length=ess_subchain_length)
        ess_mean = np.mean(ess_sweep_values)
        ess_convergence = ess_mean > ess_min_threshold

        # CIPSRF: track convergence over chain growth using first and last thirds
        first_one_third = chain[: len(chain) // 3]
        last_one_third = chain[2 * len(chain) // 3 :]
        cipsrf_ratios, cipsrf_values = cipsrf([first_one_third, last_one_third], steps=cipsrf_steps)
        cipsrf_final = cipsrf_values[-1]
        cipsrf_lower = 1.0 - cipsrf_tolerance
        cipsrf_upper = 1.0 + cipsrf_tolerance
        cipsrf_convergence = cipsrf_lower <= cipsrf_final <= cipsrf_upper

        # Overall convergence: all three tests must pass
        convergence = geweke_convergence and ess_convergence and cipsrf_convergence
        results.append(
            {
                "chain_id": cid,
                "geweke_convergence": geweke_convergence,
                "geweke_max_zscore": float(max_zscore),
                "ess_convergence": ess_convergence,
                "ess_mean": float(ess_mean),
                "cipsrf_convergence": cipsrf_convergence,
                "cipsrf_final": float(cipsrf_final),
                "overall_convergence": convergence,
            }
        )
    return results
        