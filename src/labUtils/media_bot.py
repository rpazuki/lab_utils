##############################################################
#  labUtils: A python package for parsing od600 data with    #
#  provided metadata                                         #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2025                          #
#  License: MIT                                              #
##############################################################

from __future__ import annotations

import re
from io import StringIO
from pathlib import Path
from typing import Union

import pandas as pd

__all__ = [
    'parse_raw_CLARIOstar_export',
    'parse_protocol_metadata',
    'parse',
    'report',
    'parse_time_label',
    'calculate_replicate_statistics',
]


# ---------- helpers ----------
def find_header_row(lines):
    """Locate the measurement table header in a BMG CLARIOstar export."""
    for i, line in enumerate(lines):
        if line.strip().startswith("Well Row,Well Col,Content,Raw Data"):
            return i
    raise ValueError("Could not find the data header row (Well Row,Well Col,Content,...) in the raw CSV.")

_time_pat = re.compile(r"^\s*(?:(?P<h>\d+)\s*h)?\s*(?:(?P<m>\d+)\s*min)?\s*$")

def parse_time_label(lbl: str):
    """Turn '3 h 15 min' -> (3.25 h, 3, 15). Handle '0 h 60 min' etc."""
    if lbl is None:
        return 0.0, 0, 0
    s = str(lbl).strip().replace("\u2009"," ").replace("\xa0"," ")
    m = _time_pat.match(s)
    if not m:
        return None, None, None
    h = int(m.group("h")) if m.group("h") else 0
    mnt = int(m.group("m")) if m.group("m") else 0
    h += mnt // 60
    mnt = mnt % 60
    return h + mnt/60.0, h, mnt

def parse_raw_CLARIOstar_export(path: Path, value_column_name: str = "od") -> pd.DataFrame:
    """Parse CLARIOstar OD600 export and process the header to create a tidy long format dataframe..



    Parameters
    ----------
    path : Path
        Path to the CLARIOstar OD600 export file
    value_column_name : str, optional
        Name for the value column in the output DataFrame (default: "od")

    Returns
    -------
     pd.DataFrame
        Tidy long format dataframe with processed time labels and well information
    """
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    hdr_idx = find_header_row(lines)

    # Read the block starting at the header row
    block = "\n".join(lines[hdr_idx:])
    df = pd.read_csv(StringIO(block), header=None)

    # Row 0 -> header names; Row 1 -> time labels
    time_row = df.iloc[1].tolist()
    # Time labels start from column 3 (after Well Row, Well Col, Content)
    cols_fixed = ["well_row", "well_col", "content"]
    times = []
    for t in time_row[3:]:
        t = "" if (pd.isna(t)) else str(t).strip()
        if t.lower() == "time" or t == "":
            continue
        times.append(t)

    # Fallback: parse the raw second line directly if needed
    n_data_cols = df.shape[1] - 3
    if len(times) < n_data_cols:
        second_line = lines[hdr_idx + 1]
        raw_cells = [c.strip() for c in second_line.split(",")]
        tail = raw_cells[3:]
        if tail and tail[0].lower() == "time":
            tail = tail[1:]
        times = [c for c in tail]

    times = times[:n_data_cols]
    final_cols = cols_fixed + times

    # Drop the two header rows and set final column names
    df = df.iloc[2:].reset_index(drop=True)
    df = df.iloc[:, :len(final_cols)].copy()
    df.columns = final_cols

    # Normalize types
    df["well_row"] = df["well_row"].astype(str).str.strip()
    df["well_col"] = pd.to_numeric(df["well_col"], errors="coerce").astype("Int64")

    # Wide -> long
    value_vars = [c for c in df.columns if c not in ("well_row","well_col","content")]
    long_df = df.melt(id_vars=["well_row","well_col","content"], value_vars=value_vars,
                      var_name="time_label", value_name=value_column_name)
    long_df[value_column_name] = pd.to_numeric(long_df[value_column_name], errors="coerce")
    long_df["well"] = long_df["well_row"].astype(str) + long_df["well_col"].astype("Int64").astype(str)

    # Parse time labels
    parsed = long_df["time_label"].apply(parse_time_label)
    long_df["time_h"] = parsed.apply(lambda x: x[0])
    long_df["time_h_int"] = parsed.apply(lambda x: x[1])
    long_df["time_min_int"] = parsed.apply(lambda x: x[2])
    long_df["time_min"] = (long_df["time_h"] * 60).round(3)

    # Sort by well then time
    long_df = long_df.sort_values(["well_row","well_col","time_min"], kind="stable").reset_index(drop=True)
    return long_df


def split_sections(meta_text: str):
    """Split the custom metadata file into sections keyed by the === Title === line."""
    sections = {}
    lines = [ln.strip() for ln in meta_text.splitlines() if ln.strip() != ""]
    idxs = [i for i, ln in enumerate(lines) if ln.startswith("===") and ln.endswith("===")]
    idxs.append(len(lines))
    for s, e in zip(idxs, idxs[1:]):
        title = lines[s].strip("= ").strip()
        block = [ln for ln in lines[s+1:e] if not (ln.startswith("===") and ln.endswith("==="))]
        sections[title] = "\n".join(block)
    return sections

def parse_protocol_metadata(meta_path: Path) -> pd.DataFrame:
    """Parse the '=== Experiment Data ===' section as a DataFrame and tidy column names."""
    text = meta_path.read_text(encoding="utf-8", errors="ignore")
    sections = split_sections(text)
    key = next((k for k in sections if "Experiment Data" in k), None)
    if key is None:
        raise ValueError('Could not find "Experiment Data" section in the metadata file.')
    df = pd.read_csv(StringIO(sections[key]))
    df.columns = [c.strip() for c in df.columns]
    if "Well" not in df.columns:
        raise ValueError('Experiment Data section is missing a "Well" column.')

    vol_map = {
        "Media Volume (µL)": "media_volume_uL",
        "Water Volume (µL)": "water_volume_uL",
        "Supplement Volume (µL)": "supplement_volume_uL",
        "Volume per Supplement (µL)": "volume_per_supplement_uL",
        "Total Volume (µL)": "total_volume_uL",
    }
    out = df.copy()
    for orig, new in vol_map.items():
        if orig in out.columns:
            out[new] = pd.to_numeric(out[orig], errors="coerce")

    for c in ["Well","Strain","Strain Well","Media Type","Supplements"]:
        if c in out.columns:
            out[c] = out[c].astype(str).str.strip()

    keep = ["Well","Strain","Strain Well","Media Type","Supplements"] + list(vol_map.values())
    keep = [c for c in keep if c in out.columns]
    tidy = out[keep].drop_duplicates(subset=["Well"]).rename(columns={
        "Well": "well",
        "Strain": "strain",
        "Strain Well": "strain_well",
        "Media Type": "media_type",
        "Supplements": "supplements",
    }).reset_index(drop=True)
    tidy["is_blank"] = tidy["strain"].str.lower().eq("blank")
    return tidy


def parse(raw_data: Union[pd.DataFrame, Path],
          meta_data: Union[pd.DataFrame, Path],
          value_column_name: str = "od") -> pd.DataFrame:
    """Parse raw OD600 data and metadata, return merged tidy DataFrame.

    Parameters
    ----------
    raw_data : Union[pd.DataFrame, Path]
        Either a preprocessed plate reader DataFrame from process_bmg_dataframe,
        or a Path to a raw data file that will be processed
    meta_data : Union[pd.DataFrame, Path]
        Either a preprocessed metadata DataFrame from parse_meta_experiment,
        or a Path to a metadata file that will be processed
    value_column_name : str, optional
        Name for the value column in the output DataFrame (default: "od")

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with both raw data and metadata
    """
    # Handle file paths if provided
    if isinstance(raw_data, Path):
        raw_long = parse_raw_CLARIOstar_export(raw_data, value_column_name=value_column_name)
    else:
        raw_long = raw_data

    if isinstance(meta_data, Path):
        meta = parse_protocol_metadata(meta_data)
    else:
        meta = meta_data

    # Merge data
    df = raw_long.merge(meta, on="well", how="left", validate="m:1")

    ordered_cols = [
        "well","well_row","well_col","content","strain","strain_well","is_blank","media_type","supplements",
        "media_volume_uL","water_volume_uL","supplement_volume_uL","volume_per_supplement_uL","total_volume_uL",
        "time_label","time_h","time_min_int","time_h_int","time_min",value_column_name
    ]
    ordered_cols = [c for c in ordered_cols if c in df.columns]
    df = df[ordered_cols].copy()
    return df

def report(raw_data: Union[pd.DataFrame, Path],
           meta_data: Union[pd.DataFrame, Path],
           value_column_name: str = "od") -> pd.DataFrame:
    """Generate a validation report comparing raw data and metadata.

    Parameters
    ----------
    raw_data : Union[pd.DataFrame, Path]
        Either a preprocessed plate reader DataFrame from process_bmg_dataframe,
        or a Path to a raw data file that will be processed
    meta_data : Union[pd.DataFrame, Path]
        Either a preprocessed metadata DataFrame from parse_meta_experiment,
        or a Path to a metadata file that will be processed
    value_column_name : str, optional
        Name of the value column to check for non-numeric values (default: "od")

    Returns
    -------
    pd.DataFrame
        Report of any issues found in the data
    """
    # Handle file paths if provided
    if isinstance(raw_data, Path):
        raw_long = parse_raw_CLARIOstar_export(raw_data, value_column_name=value_column_name)
    else:
        raw_long = raw_data

    if isinstance(meta_data, Path):
        meta = parse_protocol_metadata(meta_data)
    else:
        meta = meta_data
    report_rows = []
    raw_wells = set(raw_long["well"].unique())
    meta_wells = set(meta["well"].unique())
    miss_meta = sorted(raw_wells - meta_wells)
    if miss_meta:
        report_rows.append({"issue": "wells_missing_in_metadata", "count": len(miss_meta),
                            "wells": ";".join(miss_meta[:96])})
    miss_raw = sorted(meta_wells - raw_wells)
    if miss_raw:
        report_rows.append({"issue": "wells_missing_in_raw", "count": len(miss_raw),
                            "wells": ";".join(miss_raw[:96])})
    non_num = raw_long[value_column_name].isna().sum()
    if non_num:
        report_rows.append({"issue": f"non_numeric_{value_column_name}_values", "count": int(non_num)})
    dup_meta = meta["well"].duplicated().sum()
    if dup_meta:
        report_rows.append({"issue": "duplicate_wells_in_metadata", "count": int(dup_meta)})

    return pd.DataFrame(report_rows)


def calculate_replicate_statistics(
    df_parsed: pd.DataFrame,
    direction: str = "alphabetical",
    sample_size: int = 3,
    ddof: int = 1,
    value_column_name: str = "od"
) -> pd.DataFrame:
    """Calculate mean and standard deviation for groups of replicate wells.

    This function groups wells based on their position in the plate and calculates
    statistics across replicates for each measurement time point. Wells are grouped
    either along rows (alphabetical direction) or columns (numerical direction).

    Parameters
    ----------
    df_parsed : pd.DataFrame
        DataFrame returned by the `parse()` function, must contain columns:
        'well', 'well_row', 'well_col', value column, and time-related columns
    direction : str, optional
        Direction for grouping wells:
        - "alphabetical" or "alpha": Group wells along rows (e.g., A1, A2, A3)
        - "numerical" or "num": Group wells along columns (e.g., A1, B1, C1)
        Default is "alphabetical"
    sample_size : int, optional
        Number of wells to group together as replicates. Default is 3
    ddof : int, optional
        Delta Degrees of Freedom for standard deviation calculation:
        - ddof=1 (default): Unbiased estimator (sample standard deviation, N-1)
        - ddof=0: Biased estimator (population standard deviation, N)
    value_column_name : str, optional
        Name of the value column to calculate statistics on (default: "od")
        Must match the column name in the input DataFrame

    Returns
    -------
    pd.DataFrame
        DataFrame with averaged statistics containing:
        - group_id: Identifier for the replicate group (e.g., "A_1-3", "ABC_1")
        - wells: Comma-separated list of wells included in the group (e.g., "A1,A2,A3")
        - well_rows: Comma-separated list of well rows in the group (e.g., "A,A,A")
        - well_cols: Comma-separated list of well columns in the group (e.g., "1,2,3")
        - {value_column_name}_mean: Mean value across replicates
        - {value_column_name}_std: Standard deviation across replicates
        - n_replicates: Number of non-null values used in calculation
        - Time columns and metadata columns (if present in input)

    Raises
    ------
    ValueError
        - If required columns are missing from the DataFrame
        - If direction is not 'alphabetical'/'alpha' or 'numerical'/'num'
        - If the number of rows/columns is not divisible by sample_size
        - If sample_size is less than 1

    Examples
    --------
    >>> # Average 3 wells along rows (A1, A2, A3), (A4, A5, A6), etc.
    >>> stats_df = calculate_replicate_statistics(
    ...     df, direction="alphabetical", sample_size=3
    ... )

    >>> # Average 4 wells along columns (A1, B1, C1, D1), etc.
    >>> stats_df = calculate_replicate_statistics(
    ...     df, direction="numerical", sample_size=4
    ... )

    >>> # Use population standard deviation (biased)
    >>> stats_df = calculate_replicate_statistics(
    ...     df, direction="alphabetical", sample_size=3, ddof=0
    ... )

    Notes
    -----
    - Wells are sorted before grouping to ensure consistent pairing
    - For 'alphabetical' direction: sorts by well_row then well_col
    - For 'numerical' direction: sorts by well_col then well_row
    - NaN values in the value column are excluded from mean and std calculations
    - Metadata from the first well in each group is retained in the output
    """
    # Validate input DataFrame
    required_cols = ["well", "well_row", "well_col", value_column_name]
    missing_cols = [col for col in required_cols if col not in df_parsed.columns]
    if missing_cols:
        raise ValueError(f"DataFrame is missing required columns: {missing_cols}")

    # Validate direction parameter
    direction = direction.lower()
    if direction not in ["alphabetical", "alpha", "numerical", "num"]:
        raise ValueError(
            f"Invalid direction '{direction}'. Must be 'alphabetical'/'alpha' or 'numerical'/'num'"
        )

    # Validate sample_size
    if sample_size < 1:
        raise ValueError(f"sample_size must be at least 1, got {sample_size}")

    # Determine grouping direction
    is_alphabetical = direction in ["alphabetical", "alpha"]

    # Get unique wells and sort them appropriately
    unique_wells = df_parsed[["well", "well_row", "well_col"]].drop_duplicates()

    if is_alphabetical:
        # Sort by row then column (A1, A2, A3, ..., B1, B2, B3, ...)
        unique_wells = unique_wells.sort_values(["well_row", "well_col"]).reset_index(drop=True)
        n_positions = len(unique_wells["well_row"].unique())
        position_name = "rows"
    else:
        # Sort by column then row (A1, B1, C1, ..., A2, B2, C2, ...)
        unique_wells = unique_wells.sort_values(["well_col", "well_row"]).reset_index(drop=True)
        n_positions = len(unique_wells["well_col"].unique())
        position_name = "columns"

    # Check if divisible by sample_size
    if len(unique_wells) % sample_size != 0:
        raise ValueError(
            f"Number of wells ({len(unique_wells)}) is not divisible by sample_size ({sample_size}). "
            f"Total {position_name}: {n_positions}"
        )

    # Create group assignments
    unique_wells["group_num"] = unique_wells.index // sample_size

    # Merge group assignments back to original dataframe
    df_with_groups = df_parsed.merge(unique_wells[["well", "group_num"]], on="well", how="left")

    # Identify time columns for grouping
    time_cols = [col for col in ["time_label", "time_h", "time_min", "time_h_int", "time_min_int"]
                 if col in df_parsed.columns]

    # Group by group_num and time, then calculate statistics
    group_cols = ["group_num"] + time_cols

    # Calculate statistics - only include mean, std, and well information
    agg_dict = {
        value_column_name: [
            (f"{value_column_name}_mean", "mean"),
            (f"{value_column_name}_std", lambda x: x.std(ddof=ddof)),
            ("count", "count")
        ],
        "well": [("wells", lambda x: ",".join(sorted(set(x))))],
        "well_row": [("well_rows", lambda x: ",".join(sorted(set(x))))],
        "well_col": [("well_cols", lambda x: ",".join(sorted(set(map(str, x)))))],
    }

    # Include metadata columns from first well in each group
    metadata_cols = [col for col in df_parsed.columns if col not in
                    required_cols + time_cols + ["content", "group_num"]]
    for col in metadata_cols:
        if col in df_parsed.columns:
            agg_dict[col] = "first"

    stats_df = df_with_groups.groupby(group_cols, as_index=False).agg(agg_dict)

    # Flatten column names - handle both tuple and string column names
    new_cols = []
    for col in stats_df.columns:
        if isinstance(col, tuple):
            if col[0] == value_column_name:
                new_cols.append(col[1])  # Use the renamed name directly (e.g., od_mean, od_std)
            elif col[0] in ["well", "well_row", "well_col"]:
                new_cols.append(col[1])  # Use the renamed names (wells, well_rows, well_cols)
            else:
                new_cols.append(col[0] if not col[1] else f"{col[0]}_{col[1]}")
        else:
            new_cols.append(col)
    stats_df.columns = new_cols

    # Rename count column
    rename_map = {
        "count": "n_replicates",
    }
    stats_df = stats_df.rename(columns=rename_map)

    # Create group_id column
    def create_group_id(row):
        wells_list = row["wells"].split(",")
        if is_alphabetical:
            # For alphabetical: use row letter and column range (e.g., "A_1-3")
            row_letter = wells_list[0][0]  # Get row letter from first well
            col_nums = [int(w[1:]) for w in wells_list]
            if len(col_nums) == 1:
                return f"{row_letter}_{col_nums[0]}"
            return f"{row_letter}_{min(col_nums)}-{max(col_nums)}"
        else:
            # For numerical: use row range and column number (e.g., "ABC_1")
            row_letters = "".join([w[0] for w in wells_list])
            col_num = wells_list[0][1:]  # Get column number from first well
            return f"{row_letters}_{col_num}"

    stats_df["group_id"] = stats_df.apply(create_group_id, axis=1)

    # Reorder columns for better readability
    first_cols = ["group_id",
                  "wells",
                  "well_rows",
                  "well_cols",
                  f"{value_column_name}_mean",
                  f"{value_column_name}_std",
                  "n_replicates"]
    first_cols = [col for col in first_cols if col in stats_df.columns]
    other_cols = [col for col in stats_df.columns if col not in first_cols]
    stats_df = stats_df[first_cols + other_cols]

    # Drop group_num as it's internal
    if "group_num" in stats_df.columns:
        stats_df = stats_df.drop(columns=["group_num"])

    # Sort by group_id and time columns
    sort_cols = ["group_id"]
    # Add time columns in logical order if they exist
    for time_col in ["time_min", "time_h", "time_label"]:
        if time_col in stats_df.columns:
            sort_cols.append(time_col)
            break  # Use first available time column for sorting

    stats_df = stats_df.sort_values(sort_cols).reset_index(drop=True)

    return stats_df
