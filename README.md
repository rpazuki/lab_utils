# lab_utils

Utilities for experiments in computational biology that operate on laboratory data. This repository collects small, reusable tools and notebooks that help with data cleaning, QC, transformation, and common analysis workflows so collaborators can use and share them across projects.

## Features
- Data ingestion and standardization helpers for common lab file formats (CSV, TSV, Excel).
- Quality-control and normalization routines for experimental measurements.


## Repository layout
- `src/` — library modules (importable utilities).
- `tests/` — automated tests.

## Quick start
1. Clone the repo:
    ```
    git clone <repo-url>
    cd lab_utils
    ```
2. Create a virtual environment and install dependencies:
    ```
    python -m venv .venv
    source .venv/bin/activate   # or .venv\Scripts\activate on Windows
    pip install -r requirements.txt
    ```


## Usage (import example)
```
from labUtils.media_bot import parse
# Step 1: Parse raw data and metadata
df = parse(test_data_path, test_meta_path)
# Step 2: Fit growth curves for each well
params_df, preds_df = fit_modified_gompertz_per_series(
    df,
    time_col="time_h",
    value_col="od600",
    group_cols=["well"],
    min_points=5
)
```


## Contributing
- Open an issue for bugs or feature requests.
- Fork, add tests, and submit a pull request.
- Follow existing code style and add documentation for new utilities.

## License
See the LICENSE file in this repository.

## Support
For questions or contributions, please open an issue or submit a pull request.
