"""Test script to verify the refactored fit_modified_gompertz_per_series function."""

import numpy as np
import pandas as pd

from src.labUtils.growth_rates import fit_modified_gompertz_per_series, gompertz

# Generate synthetic growth curve data
np.random.seed(42)
time = np.linspace(0, 24, 25)

# True parameters - use larger amplitude for clearer signal
true_y0 = 0.05
true_A = 2.0  # Larger amplitude
true_mu_max = 0.4
true_lambda = 3.0

# Generate data with some noise
y_true = gompertz(time, true_y0, true_A, true_mu_max, true_lambda, clip_exp=100)
y_noisy = y_true + np.random.normal(0, 0.01, len(time))  # Less noise

# Create DataFrame
df = pd.DataFrame({"time": time, "OD": y_noisy, "well": ["A1"] * len(time)})

print("Testing refactored fit_modified_gompertz_per_series function\n")
print("=" * 70)

# Test 1: No fixed parameters (fit all 4)
print("\nTest 1: No fixed parameters (fit all 4 parameters)")
print("-" * 70)
result1 = fit_modified_gompertz_per_series(df, time_col="time", value_col="OD", group_cols=["well"])
print("Columns returned:", result1.columns.tolist())
print(result1.to_string(index=False))
print(f"\nTrue values: y0={true_y0}, A={true_A}, mu_max={true_mu_max}, lambda={true_lambda}")

# Test 2: Fix y0 only (mimicking old y_0_fixed behavior)
print("\n\nTest 2: Fix y0 only (old y_0_fixed=0.05 behavior)")
print("-" * 70)
result2 = fit_modified_gompertz_per_series(
    df, time_col="time", value_col="OD", group_cols=["well"], fixed_params={"y0": 0.05}
)
print("Columns returned:", result2.columns.tolist())
print(result2.to_string(index=False))  # Test 3: Fix multiple parameters (new capability)
print("\n\nTest 3: Fix y0 and lambda (new capability)")
print("-" * 70)
result3 = fit_modified_gompertz_per_series(
    df, time_col="time", value_col="OD", group_cols=["well"], fixed_params={"y0": 0.05, "lambda": 3.0}
)
print("Columns returned:", result3.columns.tolist())
print(result3.to_string(index=False))

# Test 4: Fix all but one parameter
print("\n\nTest 4: Fix y0, A, lambda (fit only mu_max)")
print("-" * 70)
result4 = fit_modified_gompertz_per_series(
    df, time_col="time", value_col="OD", group_cols=["well"], fixed_params={"y0": 0.05, "A": 2.0, "lambda": 3.0}
)
print("Columns returned:", result4.columns.tolist())
print(result4.to_string(index=False))

print("\n" + "=" * 70)
print("All tests completed successfully!")
