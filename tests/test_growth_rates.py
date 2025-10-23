"""
The test suite includes comprehensive tests for:

1. Time label parsing (`test_parse_time_label`)
2. BMG export file reading (`test_read_bmg_export`)
3. Metadata reading (`test_read_meta_experiment`)
4. Full data integration (`test_parse_integration`)
5. Data validation reporting (`test_report_function`)
"""
import numpy as np
import pandas as pd

from labUtils.growth_rates import fit_modified_gompertz_per_series, gompertz


def test_gompertz_basic():
    """Test basic functionality of the Gompertz function"""
    # Test with typical values
    t = np.array([0, 2, 4, 6, 8, 10])
    y0, A_0, mu_max, lam = 0.1, 1.0, 0.5, 2.0
    result = gompertz(t, y0, A_0, mu_max, lam, clip_exp=50.0)

    assert isinstance(result, np.ndarray)
    assert len(result) == len(t)
    assert np.all(np.isfinite(result))

    # Test initial value
    assert np.isclose(result[0], y0, rtol=1e-3)

    # Test monotonic increase
    assert np.all(np.diff(result) >= 0)


def test_gompertz_edge_cases():
    """Test edge cases of the Gompertz function"""
    t = np.array([0], dtype=float)

    # Test with zero growth
    result = gompertz(t, y0=0.1, A_0=0, mu_max=0.5, lam=2.0, clip_exp=50.0)
    assert np.isclose(result[0], 0.1)

    # Test with zero initial value
    result = gompertz(t, y0=0, A_0=1.0, mu_max=0.5, lam=2.0, clip_exp=50.0)
    assert result[0] >= 0


def test_fit_modified_gompertz_simple():
    """Test fitting with a simple dataset"""
    # Create synthetic data
    t = np.linspace(0, 10, 20)
    y0, A, mu_max, lam = 0.1, 1.0, 0.5, 2.0
    y = gompertz(t, y0, A, mu_max, lam, clip_exp=50.0)
    # Add some noise
    np.random.seed(42)
    y += np.random.normal(0, 0.02, size=len(t))

    df = pd.DataFrame({
        'time_h': t,
        'od600': y,
        'well': ['A1'] * len(t)
    })

    params_df, preds_df = fit_modified_gompertz_per_series(df)

    # Check basic structure
    assert isinstance(params_df, pd.DataFrame)
    assert isinstance(preds_df, pd.DataFrame)
    assert len(preds_df) == len(df)

    # Check if fitting was successful
    assert params_df['success'].iloc[0]

    # Check if parameters are reasonable
    assert np.isclose(params_df['y0'].iloc[0], y0, rtol=0.2)
    assert np.isclose(params_df['A'].iloc[0], A, rtol=0.2)
    assert np.isclose(params_df['mu_max'].iloc[0], mu_max, rtol=0.2)
    assert np.isclose(params_df['lambda'].iloc[0], lam, rtol=0.2)

    # Check R² value (should be good for this synthetic data)
    assert params_df['r2'].iloc[0] > 0.95


def test_fit_modified_gompertz_multiple_series():
    """Test fitting with multiple series"""
    # Create synthetic data for two wells
    t = np.linspace(0, 10, 20)
    data = []

    for well, params in [('A1', (0.1, 1.0, 0.5, 2.0)),
                        ('A2', (0.2, 0.8, 0.6, 1.5))]:
        y = gompertz(t, *params, clip_exp=50.0)
        y += np.random.normal(0, 0.02, size=len(t))
        for ti, yi in zip(t, y):
            data.append({
                'time_h': ti,
                'od600': yi,
                'well': well
            })

    df = pd.DataFrame(data)
    params_df, preds_df = fit_modified_gompertz_per_series(df)

    # Check if we got results for both series
    assert len(params_df) == 2
    assert params_df['success'].all()
    assert len(preds_df) == len(df)


def test_fit_modified_gompertz_insufficient_data():
    """Test handling of insufficient data"""
    # Create a dataset with too few points
    df = pd.DataFrame({
        'time_h': [0, 1],
        'od600': [0.1, 0.2],
        'well': ['A1', 'A1']
    })

    params_df, _ = fit_modified_gompertz_per_series(df)

    assert len(params_df) == 1
    assert not params_df['success'].iloc[0]
    assert "insufficient" in params_df['message'].iloc[0].lower()


def test_fit_modified_gompertz_flat_data():
    """Test handling of flat data (no growth)"""
    # Create a dataset with no growth
    df = pd.DataFrame({
        'time_h': np.linspace(0, 10, 10),
        'od600': [0.1] * 10,
        'well': ['A1'] * 10
    })

    params_df, _ = fit_modified_gompertz_per_series(df)

    assert len(params_df) == 1
    assert not params_df['success'].iloc[0]
    assert "flat" in params_df['message'].iloc[0].lower()
