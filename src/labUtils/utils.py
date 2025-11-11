from pathlib import Path

import pandas as pd


def read_csv(path: Path) -> pd.DataFrame:
    """Utility function to read a CSV file into a DataFrame."""
    df = pd.read_csv(path)
    return df
