"""Test script to verify curve_smoothing works correctly."""

import numpy as np
import pandas as pd

from src.labUtils.growth_rates import fit_max_growth_rate_per_series

# Generate synthetic noisy growth curve data
np.random.seed(42)
time = np.linspace(0, 24, 50)

# Create exponential growth with noise
true_growth_rate = 0.3
y_true = np.exp(true_growth_rate * time)
y_noisy = y_true + np.random.normal(0, 0.1, len(time))

# Create DataFrame
df = pd.DataFrame({"time": time, "OD": y_noisy, "well": ["A1"] * len(time)})

print("Testing fit_max_growth_rate_per_series with smoothing\n")
print("=" * 70)

# Test 1: No smoothing
print("\nTest 1: No smoothing (baseline)")
print("-" * 70)
result1 = fit_max_growth_rate_per_series(
    df, time_col="time", value_col="OD", group_cols=["well"], moving_window_size=5, smoothing_iterations=0
)
print(result1[["well", "mu_max", "r2", "rmse", "success"]].to_string(index=False))
print(f"True growth rate: {true_growth_rate}")

# Test 2: With 1 iteration of smoothing
print("\n\nTest 2: With 1 smoothing iteration")
print("-" * 70)
result2 = fit_max_growth_rate_per_series(
    df,
    time_col="time",
    value_col="OD",
    group_cols=["well"],
    moving_window_size=5,
    smoothing_iterations=1,
    smooth_window_size=3,
)
print(result2[["well", "mu_max", "r2", "rmse", "success"]].to_string(index=False))

# Test 3: With 2 iterations of smoothing
print("\n\nTest 3: With 2 smoothing iterations")
print("-" * 70)
result3 = fit_max_growth_rate_per_series(
    df,
    time_col="time",
    value_col="OD",
    group_cols=["well"],
    moving_window_size=5,
    smoothing_iterations=2,
    smooth_window_size=3,
)
print(result3[["well", "mu_max", "r2", "rmse", "success"]].to_string(index=False))

# Test 4: Edge case - small dataset
print("\n\nTest 4: Edge case - Small dataset with smoothing")
print("-" * 70)
df_small = df.head(10)
result4 = fit_max_growth_rate_per_series(
    df_small,
    time_col="time",
    value_col="OD",
    group_cols=["well"],
    moving_window_size=5,
    smoothing_iterations=1,
    smooth_window_size=3,
)
print(result4[["well", "mu_max", "success", "message"]].to_string(index=False))

print("\n" + "=" * 70)
print("All smoothing tests completed!")
