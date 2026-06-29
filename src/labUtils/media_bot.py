##############################################################
#  labUtils: A python package for parsing od600 data with    #
#  provided metadata                                         #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2025                          #
#  License: MIT                                              #
##############################################################

from __future__ import annotations

import logging
import csv
import re
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Any, TypedDict

import pandas as pd

__all__ = [
    "ColumnSchema",
    "ParseResult",
    "ReplicateStats",
    "CustomReplicateRule",
    "parse_raw_CLARIOstar_export",
    "parse_protocol_metadata",
    "parse",
    "parse_result",
    "report",
    "parse_time_label",
    "calculate_replicate_statistics_by_well",
    "calculate_replicate_statistics_by_custom",
]


# ---------------------------------------------------------------------------
# Output containers (synthetic.py-style)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ColumnSchema:
    """Configurable column names used by parsing / replicate-stats functions.

    Defaults match the current pipeline; pass a custom instance to rename columns
    end-to-end without touching every call site.
    """

    well: str = "well"
    well_row: str = "well_row"
    well_col: str = "well_col"
    strain: str = "strain"
    is_blank: str = "is_blank"
    media_type: str = "media_type"
    value: str = "od"


@dataclass
class ParseResult:
    """Output of `parse_result`. Wraps the merged dataframe plus its inputs."""

    df: pd.DataFrame
    raw_df: pd.DataFrame
    meta_df: pd.DataFrame
    value_column_name: str = "od"


@dataclass
class ReplicateStats:
    """Output of `calculate_replicate_statistics_*` (dataclass-returning sibling)."""

    df: pd.DataFrame
    grouping: str
    sample_size_used: int | dict[str, int] = 0
    direction: str | None = None


class CustomReplicateRule(TypedDict, total=False):
    """Rule definition for custom replicate statistics aggregation."""

    direction: str
    pattern: str
    sample_size: int


# ---------- helpers ----------
def find_header_row(lines):
    """Locate the measurement table header in a BMG CLARIOstar export."""
    logging.info("Function labUtils.media_bot.find_header_row is started")
    for i, line in enumerate(lines):
        norm = re.sub(r"\s*,\s*", ",", line.strip()).lower()
        if norm.startswith("well row,well col,content,raw data"):
            logging.info("Function labUtils.media_bot.find_header_row finished successfully")
            return i
        if norm.startswith("well,content,raw data"):
            logging.info("Function labUtils.media_bot.find_header_row finished successfully")
            return i
    raise ValueError("Could not find the data header row (Well Row,Well Col,Content,...) in the file.")


def _line_is_nonempty(line: str) -> bool:
    """Return True if a CSV line (or an Excel row joined as CSV) has any content."""
    return bool(line.replace(",", "").strip())


_time_pat = re.compile(r"^\s*(?:(?P<h>\d+)\s*h)?\s*(?:(?P<m>\d+)\s*min)?\s*$")


def parse_time_label(lbl: str):
    """Turn '3 h 15 min' -> (3.25 h, 3, 15). Handle '0 h 60 min' etc."""
    if lbl is None:
        return 0.0, 0, 0
    s = str(lbl).strip().replace("\u2009", " ").replace("\xa0", " ")
    m = _time_pat.match(s)
    if not m:
        return None, None, None
    h = int(m.group("h")) if m.group("h") else 0
    mnt = int(m.group("m")) if m.group("m") else 0
    h += mnt // 60
    mnt = mnt % 60
    return h + mnt / 60.0, h, mnt


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
    logging.info("Function labUtils.media_bot.parse_raw_CLARIOstar_export is started")
    if path.suffix.lower() in (".xlsx", ".xls"):
        df_raw = pd.read_excel(path, header=None, dtype=str).fillna("")
        lines = [",".join(str(c) for c in row) for row in df_raw.values.tolist()]
        hdr_idx = find_header_row(lines)
        rows = df_raw.iloc[hdr_idx:].values.tolist()
    else:
        text = path.read_text(encoding="utf-8", errors="ignore")
        lines = text.splitlines()
        hdr_idx = find_header_row(lines)
        block = "\n".join(lines[hdr_idx:])
        rows = list(csv.reader(StringIO(block)))

    rows = [r for r in rows if any(str(cell).strip() != "" for cell in r)]
    if not rows:
        raise ValueError("No rows found after the detected header row.")
    max_cols = max(len(r) for r in rows)
    rows = [r + [""] * (max_cols - len(r)) for r in rows]
    df = pd.DataFrame(rows)

    # Detect raw-data block structure from the header row
    header_row = ["" if pd.isna(x) else str(x).strip() for x in df.iloc[0].tolist()]
    header_norm = [h.lower() for h in header_row]

    raw_data_candidates = [i for i, h in enumerate(header_norm) if h.startswith("raw data")]
    if not raw_data_candidates:
        raise ValueError("Could not identify the 'Raw Data' column in the raw CSV header row.")

    raw_data_idx = raw_data_candidates[0]

    has_split_well = "well row" in header_norm and "well col" in header_norm
    has_combined_well = "well" in header_norm
    content_idx = header_norm.index("content") if "content" in header_norm else None

    if has_split_well:
        well_row_idx = header_norm.index("well row")
        well_col_idx = header_norm.index("well col")
    elif has_combined_well:
        well_idx = header_norm.index("well")
    else:
        raise ValueError("Could not identify well columns. Expected 'Well Row'/'Well Col' or 'Well'.")

    # Row 0 -> header names; Row 1 -> time labels
    time_row = df.iloc[1].tolist()
    times = []
    for t in time_row[raw_data_idx:]:
        t = "" if (pd.isna(t)) else str(t).strip()
        if t.lower() == "time" or t == "":
            continue
        times.append(t)

    # Fallback: parse the raw second line directly if needed
    n_data_cols = df.shape[1] - raw_data_idx
    if len(times) < n_data_cols:
        next_non_empty_idx = next((i for i in range(hdr_idx + 1, len(lines)) if _line_is_nonempty(lines[i])), None)
        if next_non_empty_idx is None:
            raise ValueError("Could not find a non-empty time-label row after the data header row.")
        second_line = lines[next_non_empty_idx]
        raw_cells = [c.strip() for c in second_line.split(",")]
        tail = raw_cells[raw_data_idx:]
        if tail and tail[0].lower() == "time":
            tail = tail[1:]
        times = [c for c in tail]  # noqa: C416

    times = times[:n_data_cols]

    # Drop the two header rows and keep only columns needed for parsing
    n_keep = raw_data_idx + len(times)
    data_rows = df.iloc[2:].reset_index(drop=True)
    data_rows = data_rows.iloc[:, :n_keep].copy()

    out_df = pd.DataFrame(index=data_rows.index)
    if content_idx is not None and content_idx < data_rows.shape[1]:
        out_df["content"] = data_rows.iloc[:, content_idx]
    else:
        out_df["content"] = pd.NA

    if has_split_well:
        out_df["well_row"] = data_rows.iloc[:, well_row_idx].astype(str).str.strip()
        out_df["well_col"] = pd.to_numeric(data_rows.iloc[:, well_col_idx], errors="coerce").astype("Int64")
    else:
        well_text = data_rows.iloc[:, well_idx].astype(str).str.strip()
        extracted = well_text.str.extract(r"^\s*([A-Za-z]+)\s*(\d+)\s*$")
        out_df["well_row"] = extracted[0].str.strip()
        out_df["well_col"] = pd.to_numeric(extracted[1], errors="coerce").astype("Int64")

    if times:
        time_values = data_rows.iloc[:, raw_data_idx : raw_data_idx + len(times)].copy()
        time_values.columns = times
        out_df = pd.concat([out_df.reset_index(drop=True), time_values.reset_index(drop=True)], axis=1)

    # Wide -> long
    value_vars = [c for c in out_df.columns if c not in ("well_row", "well_col", "content")]
    long_df = out_df.melt(
        id_vars=["well_row", "well_col", "content"],
        value_vars=value_vars,
        var_name="time_label",
        value_name=value_column_name,
    )
    long_df[value_column_name] = pd.to_numeric(long_df[value_column_name], errors="coerce")
    long_df["well"] = long_df["well_row"].astype(str) + long_df["well_col"].astype("Int64").astype(str)

    # Parse time labels
    parsed = long_df["time_label"].apply(parse_time_label)
    long_df["time_h"] = parsed.apply(lambda x: x[0])
    long_df["time_h_int"] = parsed.apply(lambda x: x[1])
    long_df["time_min_int"] = parsed.apply(lambda x: x[2])
    long_df["time_min"] = (long_df["time_h"] * 60).round(3)

    # Sort by well then time
    long_df = long_df.sort_values(["well_row", "well_col", "time_min"], kind="stable").reset_index(drop=True)
    logging.info("Function labUtils.media_bot.parse_raw_CLARIOstar_export finished successfully")
    return long_df


def _normalize_section_marker(line: str) -> str:
    """Strip trailing commas/whitespace that some exports append to `=== Title ===`."""
    return re.sub(r"\s*,+\s*$", "", line.strip())
def split_sections(meta_text: str):
    """Split the custom metadata file into sections keyed by the === Title === line."""
    logging.info("Function labUtils.media_bot.split_sections is started")
    sections = {}

    def _normalize_marker_line(line: str) -> str:
        # Some exports add trailing commas after section headers.
        return re.sub(r"\s*,+\s*$", "", line.strip())


def _is_section_marker(line: str) -> bool:
    normalized = _normalize_section_marker(line)
    return normalized.startswith("===") and normalized.endswith("===")


def split_sections(meta_text: str):
    """Split the custom metadata file into sections keyed by the === Title === line."""
    sections = {}
    lines = [ln.strip() for ln in meta_text.splitlines() if ln.strip() != "" and not re.fullmatch(r"[\s,]+", ln)]
    idxs = [i for i, ln in enumerate(lines) if _is_section_marker(ln)]
    idxs.append(len(lines))
    for s, e in zip(idxs, idxs[1:]):
        title = _normalize_section_marker(lines[s]).strip("= ").strip()
        block = [ln for ln in lines[s + 1 : e] if not _is_section_marker(ln)]
        sections[title] = "\n".join(block)
    logging.info("Function labUtils.media_bot.split_sections finished successfully")
    return sections


def parse_protocol_metadata(meta_path: Path) -> pd.DataFrame:
    """Parse the '=== Experiment Data ===' section as a DataFrame and tidy column names."""
    logging.info("Function labUtils.media_bot.parse_protocol_metadata is started")
    text = meta_path.read_text(encoding="utf-8", errors="ignore")
    sections = split_sections(text)
    df = None

    def _try_parse_csv_block(csv_text: str) -> pd.DataFrame | None:
        cleaned_lines = []
        for ln in csv_text.splitlines():
            stripped = ln.strip()
            if not stripped:
                continue
            if re.fullmatch(r"[\s,]+", stripped):
                continue
            marker_candidate = re.sub(r"\s*,+\s*$", "", stripped)
            if marker_candidate.startswith("===") and marker_candidate.endswith("==="):
                continue
            cleaned_lines.append(stripped)
        if not cleaned_lines:
            return None

        try:
            candidate = pd.read_csv(StringIO("\n".join(cleaned_lines)))
        except Exception:
            return None
        candidate.columns = [str(c).strip() for c in candidate.columns]
        return candidate if "Well" in candidate.columns else None
def _make_group_id(row: pd.Series, is_alphabetical: bool) -> str:
    """Format the per-row replicate group label (e.g. `A_1-3` or `ABC_1`).

    Promoted from a nested closure so both `calculate_replicate_statistics_*`
    functions can share the exact same formatting logic.
    """
    wells_list = row["wells"].split(",")
    if is_alphabetical:
        row_letter = wells_list[0][0]
        col_nums = [int(w[1:]) for w in wells_list]
        if len(col_nums) == 1:
            return f"{row_letter}_{col_nums[0]}"
        return f"{row_letter}_{min(col_nums)}-{max(col_nums)}"
    row_letters = "".join([w[0] for w in wells_list])
    col_num = wells_list[0][1:]
    return f"{row_letters}_{col_num}"


def _load_config(config: dict | str | Path) -> dict:
    """Load a media_bot pipeline config (yaml file or dict).

    Mirrors `synthetic._load_config`. Accepts a dict directly or a path to a yaml
    file; unwraps a top-level `media_bot:` key if present. Currently unused by
    the public API — exposed for downstream pipelines that want to declare their
    parsing / replicate-stats settings in yaml.
    """
    if isinstance(config, (str, Path)):
        import yaml  # local import — yaml is optional for the rest of the module

        with open(config, encoding="utf-8") as f:
            data = yaml.safe_load(f)
    else:
        data = config
    if not isinstance(data, dict):
        raise ValueError("media_bot config must be a mapping at top level.")
    if "media_bot" in data:
        data = data["media_bot"]
    return data


def _try_csv_block(csv_text: str) -> pd.DataFrame | None:
    """Best-effort CSV parser used inside `parse_protocol_metadata` fallbacks."""
    cleaned_lines = []
    for ln in csv_text.splitlines():
        stripped = ln.strip()
        if not stripped:
            continue
        if re.fullmatch(r"[\s,]+", stripped):
            continue
        marker_candidate = re.sub(r"\s*,+\s*$", "", stripped)
        if marker_candidate.startswith("===") and marker_candidate.endswith("==="):
            continue
        cleaned_lines.append(stripped)
    if not cleaned_lines:
        return None
    try:
        candidate = pd.read_csv(StringIO("\n".join(cleaned_lines)))
    except Exception:
        return None
    candidate.columns = [str(c).strip() for c in candidate.columns]
    return candidate if "Well" in candidate.columns else None


def parse_protocol_metadata(meta_path: Path) -> pd.DataFrame:
    """Parse the '=== Experiment Data ===' section as a DataFrame and tidy column names."""
    text = meta_path.read_text(encoding="utf-8", errors="ignore")
    sections = split_sections(text)
    df = None

    key = next((k for k in sections if "experiment data" in k.lower()), None)
    if key is None:
        key = next((k for k in sections if "experiment" in k.lower() and "data" in k.lower()), None)
    if key is None:
        for section_title, section_text in sections.items():
            try:
                probe_df = pd.read_csv(StringIO(section_text))
            except Exception:
                continue
            probe_cols = {str(c).strip().lower() for c in probe_df.columns}
            if "well" in probe_cols:
                key = section_title
                break

    # Parse selected section when available
    if key is not None:
        df = pd.read_csv(StringIO(sections[key]))
        df.columns = [c.strip() for c in df.columns]

    # Fallback 1: file might be plain CSV with no section markers
    if df is None:
        df = _try_csv_block(text)

    # Fallback 2: file may have preamble lines before the CSV header
    if df is None:
        raw_lines = text.splitlines()
        start_idx = next(
            (
                i
                for i, ln in enumerate(raw_lines)
                if "," in ln and any(cell.strip().lower() == "well" for cell in ln.split(","))
            ),
            None,
        )
        if start_idx is not None:
            df = _try_csv_block("\n".join(raw_lines[start_idx:]))

    if df is None:
        raise ValueError('Could not find "Experiment Data" section in the metadata file.')

    if "Well" not in df.columns:
        raise ValueError('Experiment Data section is missing a "Well" column.')

    # Guard against section labels or comma-only rows leaking into parsed data.
    df["Well"] = df["Well"].astype(str).str.strip()
    df = df[df["Well"].str.fullmatch(r"[A-Za-z]+\d+", na=False)].copy()

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

    for c in ["Well", "Strain", "Strain Well", "Media Type", "Supplements"]:
        if c in out.columns:
            out[c] = out[c].astype(str).str.strip()

    keep = ["Well", "Strain", "Strain Well", "Media Type", "Supplements"] + list(vol_map.values())
    keep = [c for c in keep if c in out.columns]
    tidy = (
        out[keep]
        .drop_duplicates(subset=["Well"])
        .rename(
            columns={
                "Well": "well",
                "Strain": "strain",
                "Strain Well": "strain_well",
                "Media Type": "media_type",
                "Supplements": "supplements",
            }
        )
        .reset_index(drop=True)
    )
    tidy["is_blank"] = tidy["strain"].str.lower().eq("blank")
    logging.info("Function labUtils.media_bot.parse_protocol_metadata finished successfully")
    return tidy


def parse(raw_data: pd.DataFrame | Path, meta_data: pd.DataFrame | Path, value_column_name: str = "od") -> pd.DataFrame:
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
    logging.info("Function labUtils.media_bot.parse is started")
    # Handle file paths if provided
    if isinstance(raw_data, Path):
        raw_long = parse_raw_CLARIOstar_export(raw_data, value_column_name=value_column_name)
    else:
        raw_long = raw_data

    if isinstance(meta_data, Path):  # noqa: SIM108
        meta = parse_protocol_metadata(meta_data)
    else:
        meta = meta_data

    # Merge data
    df = raw_long.merge(meta, on="well", how="left", validate="m:1")

    ordered_cols = [
        "well",
        "well_row",
        "well_col",
        "content",
        "strain",
        "strain_well",
        "is_blank",
        "media_type",
        "supplements",
        "media_volume_uL",
        "water_volume_uL",
        "supplement_volume_uL",
        "volume_per_supplement_uL",
        "total_volume_uL",
        "time_label",
        "time_h",
        "time_min_int",
        "time_h_int",
        "time_min",
        value_column_name,
    ]
    ordered_cols = [c for c in ordered_cols if c in df.columns]
    df = df[ordered_cols].copy()
    logging.info("Function labUtils.media_bot.parse finished successfully")
    return df


def parse_result(
    raw_data: pd.DataFrame | Path,
    meta_data: pd.DataFrame | Path,
    value_column_name: str = "od",
) -> ParseResult:
    """Dataclass-returning sibling of `parse`.

    Same behavior; returns a `ParseResult` carrying the merged dataframe plus
    references to the resolved raw and metadata frames.
    """
    raw_long = (
        parse_raw_CLARIOstar_export(raw_data, value_column_name=value_column_name)
        if isinstance(raw_data, Path)
        else raw_data
    )
    meta = parse_protocol_metadata(meta_data) if isinstance(meta_data, Path) else meta_data
    df = parse(raw_long, meta, value_column_name=value_column_name)
    return ParseResult(df=df, raw_df=raw_long, meta_df=meta, value_column_name=value_column_name)


def report(
    raw_data: pd.DataFrame | Path, meta_data: pd.DataFrame | Path, value_column_name: str = "od"
) -> pd.DataFrame:
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
    logging.info("Function labUtils.media_bot.report is started")
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
        report_rows.append(
            {"issue": "wells_missing_in_metadata", "count": len(miss_meta), "wells": ";".join(miss_meta[:96])}
        )
    miss_raw = sorted(meta_wells - raw_wells)
    if miss_raw:
        report_rows.append({"issue": "wells_missing_in_raw", "count": len(miss_raw), "wells": ";".join(miss_raw[:96])})
    non_num = raw_long[value_column_name].isna().sum()
    if non_num:
        report_rows.append({"issue": f"non_numeric_{value_column_name}_values", "count": int(non_num)})
    dup_meta = meta["well"].duplicated().sum()
    if dup_meta:
        report_rows.append({"issue": "duplicate_wells_in_metadata", "count": int(dup_meta)})

    logging.info("Function labUtils.media_bot.report finished successfully")
    return pd.DataFrame(report_rows)


def calculate_replicate_statistics_by_well(
    df_parsed: pd.DataFrame,
    direction: str = "alphabetical",
    sample_size: int = 3,
    ddof: int = 1,
    value_column_name: str = "od",
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
        - "alphabetical", "alpha", or "rows": Group wells along rows (e.g., A1, A2, A3)
        - "numerical", "num", or "columns": Group wells along columns (e.g., A1, B1, C1)
        Convention: alphabetical/rows = row direction, numerical/columns = column direction
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
    >>> stats_df = calculate_replicate_statistics_by_well(
    ...     df, direction="alphabetical", sample_size=3
    ... )

    >>> # Average 4 wells along columns (A1, B1, C1, D1), etc.
    >>> stats_df = calculate_replicate_statistics_by_well(
    ...     df, direction="numerical", sample_size=4
    ... )

    >>> # Use population standard deviation (biased)
    >>> stats_df = calculate_replicate_statistics_by_well(
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
    logging.info("Function labUtils.media_bot.calculate_replicate_statistics_by_well is started")
    # Validate input DataFrame
    required_cols = ["well", "well_row", "well_col", value_column_name]
    missing_cols = [col for col in required_cols if col not in df_parsed.columns]
    if missing_cols:
        raise ValueError(f"DataFrame is missing required columns: {missing_cols}")

    # Validate direction parameter
    direction = direction.lower()
    if direction not in ["alphabetical", "alpha", "rows", "numerical", "num", "columns"]:
        raise ValueError(f"Invalid direction '{direction}'. Must be 'alphabetical'/'alpha'/'rows' or 'numerical'/'num'/'columns'")

    # Validate sample_size
    if sample_size < 1:
        raise ValueError(f"sample_size must be at least 1, got {sample_size}")

    # Determine grouping direction
    is_alphabetical = direction in ["alphabetical", "alpha", "rows"]

    # Get unique wells and sort them appropriately
    unique_wells = df_parsed[["well", "well_row", "well_col"]].drop_duplicates()

    if is_alphabetical:
        # Sort by row then column (A1, A2, A3, ..., B1, B2, B3, ...)
        unique_wells = unique_wells.sort_values(["well_row", "well_col"]).reset_index(drop=True)
        n_positions = len(unique_wells["well_row"].unique())
        position_name = "rows"
        logging.info(f"Grouping wells along rows (alphabetical direction). Total unique rows: {n_positions}")
    else:
        # Sort by column then row (A1, B1, C1, ..., A2, B2, C2, ...)
        unique_wells = unique_wells.sort_values(["well_col", "well_row"]).reset_index(drop=True)
        n_positions = len(unique_wells["well_col"].unique())
        position_name = "columns"
        logging.info(f"Grouping wells along columns (numerical direction). Total unique columns: {n_positions}")

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
    time_cols = [
        col for col in ["time_label", "time_h", "time_min", "time_h_int", "time_min_int"] if col in df_parsed.columns
    ]

    # Group by group_num and time, then calculate statistics
    group_cols = ["group_num"] + time_cols

    # Calculate statistics - only include mean, std, and well information
    agg_dict = {
        value_column_name: [
            (f"{value_column_name}_mean", "mean"),
            (f"{value_column_name}_std", lambda x: x.std(ddof=ddof)),
            ("count", "count"),
        ],
        "well": [("wells", lambda x: ",".join(sorted(set(x))))],
        "well_row": [("well_rows", lambda x: ",".join(sorted(set(x))))],
        "well_col": [("well_cols", lambda x: ",".join(sorted(set(map(str, x)))))],
    }

    # Include metadata columns from first well in each group
    metadata_cols = [
        col for col in df_parsed.columns if col not in required_cols + time_cols + ["content", "group_num"]
    ]
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

    stats_df["group_id"] = stats_df.apply(
        lambda row: _make_group_id(row, is_alphabetical=is_alphabetical), axis=1
    )

    # Reorder columns for better readability
    first_cols = [
        "group_id",
        "wells",
        "well_rows",
        "well_cols",
        f"{value_column_name}_mean",
        f"{value_column_name}_std",
        "n_replicates",
    ]
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

    # Romove "_first" suffix from columns names that have it
    first_suffix_cols = [col for col in stats_df.columns if col.endswith("_first")]
    for first_col in first_suffix_cols:
        stats_df = stats_df.rename(columns={first_col: first_col[:-6]})  # Remove '_first' suffix

    logging.info("Function labUtils.media_bot.calculate_replicate_statistics_by_well finished successfully")
    return stats_df


def calculate_replicate_statistics_by_custom(
    df_parsed: pd.DataFrame,
    strain_pattern: str | None = r"([A-Za-z]+)",
    custom_rules: dict[str, CustomReplicateRule] | None = None,
    ddof: int = 1,
    value_column_name: str = "od",
) -> pd.DataFrame:
    """Calculate mean and standard deviation with custom rules per strain group.

    This function combines strain-based grouping with custom direction and sample_size
    rules for each strain. It first extracts strain groups using a regex pattern, then
    applies specific aggregation rules (direction and sample_size) defined for each strain.
    Rules can target specific strains by name or match multiple strains using a regex pattern.

    Parameters
    ----------
    df_parsed : pd.DataFrame
        DataFrame returned by the `parse()` function, must contain columns:
        'well', 'well_row', 'well_col', 'strain', value column, and time-related columns
    strain_pattern : str | None, optional
        Global regex pattern to extract the grouping key from the 'strain' column.
        Default is r"[A-Za-z]+" (extracts alphabetical part only).
        Example: "strainA1", "strainA2" -> grouped as "strainA"
        If None, uses strain values as-is without pattern extraction.
        This global pattern is applied before rule-specific patterns.
    custom_rules : dict[str, dict[str, int | str]] | None, optional
        Dictionary mapping rule names to their aggregation rules.
        Each rule can target specific strains by name or use a pattern to match multiple strains.
        Each rule dictionary can contain:
        - 'direction': "alphabetical"/"alpha"/"rows" or "numerical"/"num"/"columns" (required).
          Convention: alphabetical/rows = row direction, numerical/columns = column direction
        - 'sample_size': int >= 1 (optional for pattern rules)
        - 'pattern': str, optional regex pattern to match strain names

        Format (exact strain match):
        {
            "strainA": {"direction": "alphabetical", "sample_size": 3},
            "strainB": {"direction": "numerical", "sample_size": 4}
        }

        Format (pattern-based matching):
        {
            "control_strains": {"pattern": r"control.*", "direction": "alphabetical", "sample_size": 3},
            "experimental": {"pattern": r"st_\\d+_r\\d+", "direction": "numerical", "sample_size": 4}
        }

        Format (pattern-based matching without fixed sample_size):
        {
            "EGMB": {"pattern": r"EGMB_\\d+", "direction": "alphabetical"}
        }

        When 'pattern' is specified, all strains matching the regex pattern are processed
        as separate groups when 'sample_size' is provided. When 'pattern' is specified
        without 'sample_size', all matching strains are combined under the rule name and
        grouped by the full row/column implied by 'direction', so replicate counts can vary
        between groups. When 'pattern' is not specified, the rule name is treated as an
        exact strain name match.
        If a strain is not matched by any rule, it will be skipped.
        If None, an empty dict is used (no strains processed).
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
        - strain: Identifier for the strain group
        - group_id: Identifier for the replicate group within strain (e.g., "A_1-3", "ABC_1")
        - wells: Comma-separated list of wells included in the group
        - well_rows: Comma-separated list of well rows in the group
        - well_cols: Comma-separated list of well columns in the group
        - {value_column_name}_mean: Mean value across replicates
        - {value_column_name}_std: Standard deviation across replicates
        - n_replicates: Number of non-null values used in calculation
        - Time columns and metadata columns (if present in input)

    Raises
    ------
    ValueError
        - If required columns are missing from the DataFrame
        - If direction in custom_rules is invalid
        - If sample_size in custom_rules is less than 1 when provided
        - If the number of wells for a strain is not divisible by its sample_size
        - If pattern in custom_rules fails to compile as a regex

    Examples
    --------
    >>> # Define custom rules for specific strains by name
    >>> custom_rules = {
    ...     "strainA": {"direction": "alphabetical", "sample_size": 3},
    ...     "strainB": {"direction": "numerical", "sample_size": 4},
    ...     "control": {"direction": "alphabetical", "sample_size": 6}
    ... }
    >>> stats_df = calculate_replicate_statistics_by_custom(
    ...     df, strain_pattern=r"[A-Za-z]+", custom_rules=custom_rules
    ... )

    >>> # Use pattern-based rules to match multiple strains
    >>> custom_rules = {
    ...     "experimental_strains": {
    ...         "pattern": r"st_\\d+_r\\d+",
    ...         "direction": "alphabetical",
    ...         "sample_size": 3
    ...     },
    ...     "control": {"direction": "numerical", "sample_size": 3}
    ... }
    >>> stats_df = calculate_replicate_statistics_by_custom(
    ...     df, strain_pattern=None, custom_rules=custom_rules
    ... )

    Notes
    -----
    - Each strain group is processed independently with its own rules
    - For each strain, wells are grouped according to its direction and sample_size
    - Results from all strains are concatenated into a single DataFrame
    - Only strains matched by custom_rules are processed
    - When using patterns with sample_size, each distinct matching strain is treated as a
      separate group
    - When using patterns without sample_size, all matching strains are combined under the
      rule name and grouped by full rows or columns, depending on direction
    """
    logging.info("Function labUtils.media_bot.calculate_replicate_statistics_by_custom is started")
    # Validate input DataFrame
    required_cols = ["well", "well_row", "well_col", "strain", value_column_name]
    missing_cols = [col for col in required_cols if col not in df_parsed.columns]
    if missing_cols:
        raise ValueError(f"DataFrame is missing required columns: {missing_cols}")

    output_columns = [
        "strain",
        "group_id",
        "wells",
        "well_rows",
        "well_cols",
        f"{value_column_name}_mean",
        f"{value_column_name}_std",
        "n_replicates",
    ]

    # Handle empty custom_rules
    if custom_rules is None:
        custom_rules = {}

    if not custom_rules:
        # Return empty DataFrame with expected columns
        logging.info("Function calculate_replicate_statistics_by_custom finished successfully")
        return pd.DataFrame(columns=output_columns)

    # Create a copy to avoid modifying the original
    df_with_groups = df_parsed.copy()

    logging.info(f"strain_pattern: {strain_pattern}")
    # Apply global pattern extraction if specified
    if strain_pattern is not None:
        extracted = df_with_groups["strain"].astype(str).str.extract(f"{strain_pattern}", expand=True)
        # Use the first captured group if there are multiple
        if isinstance(extracted, pd.DataFrame):
            # join the columns of the extracted DataFrame into a single string for each row
            df_with_groups["strain_group"] = extracted.apply(lambda row: "".join(row.dropna().astype(str)), axis=1)
            # df_with_groups["strain_group"] = extracted.iloc[:, 0]
        else:
            df_with_groups["strain_group"] = extracted
    else:
        df_with_groups["strain_group"] = df_with_groups["strain"]

    # Identify time columns for grouping
    time_cols = [
        col for col in ["time_label", "time_h", "time_min", "time_h_int", "time_min_int"] if col in df_parsed.columns
    ]

    # Build a list of rule applications to process and track conflicts across rules.
    unique_strains = [str(strain) for strain in df_with_groups["strain_group"].dropna().unique()]
    logging.info(f"Unique strains: {unique_strains}")
    assigned_strains: dict[str, str] = {}
    rules_to_process: list[tuple[str, CustomReplicateRule, list[str], bool]] = []
    for rule_key, rules in custom_rules.items():
        pattern_str = rules.get("pattern")
        use_full_direction_groups = pattern_str is not None and rules.get("sample_size") is None

        # Check if this rule has a pattern
        if pattern_str is not None:
            if not isinstance(pattern_str, str):
                raise ValueError(f"Pattern for rule '{rule_key}' must be a string, got {type(pattern_str).__name__}")

            # Pattern-based rule: find all strains matching the pattern
            try:
                pattern = re.compile(pattern_str)
            except re.error as e:
                raise ValueError(f"Invalid regex pattern in rule '{rule_key}': {pattern_str}. Error: {e}")

            # Find all strain_groups matching this pattern
            matching_strains = [strain for strain in unique_strains if pattern.search(strain)]
            for strain in matching_strains:
                if strain in assigned_strains:
                    raise ValueError(
                        f"Strain '{strain}' matches multiple patterns in custom_rules. "
                        f"Please ensure each strain matches at most one pattern."
                    )
                assigned_strains[strain] = rule_key
            logging.info(f"Rule '{rule_key}' matches strains: {matching_strains} pattern_str: {pattern_str} sample_size: {rules.get('sample_size')}")
            if use_full_direction_groups:
                rules_to_process.append((rule_key, rules, matching_strains, True))
            else:
                for strain in matching_strains:
                    rules_to_process.append((strain, rules, [strain], False))
        else:
            # Exact name match: use rule_key as the strain name
            if rule_key in assigned_strains:
                raise ValueError(
                    f"Strain '{rule_key}' matches multiple patterns in custom_rules. "
                    f"Please ensure each strain matches at most one pattern."
                )
            assigned_strains[rule_key] = rule_key
            use_full_direction_groups = rules.get("sample_size") is None
            rules_to_process.append((rule_key, rules, [rule_key], use_full_direction_groups))
            logging.info(f"Rule '{rule_key}' matches exact strain: {rule_key} sample_size: {rules.get('sample_size')}")

    # Process each strain according to its custom rules
    all_results = []
    for output_strain, rules, strain_names, use_full_direction_groups in rules_to_process:
        # Filter data for this rule application
        strain_df = df_with_groups[df_with_groups["strain_group"].isin(strain_names)].copy()

        if len(strain_df) == 0:
            continue  # Skip if no data for this strain

        # Extract rules
        direction_raw = rules.get("direction", "alphabetical")
        if not isinstance(direction_raw, str):
            raise ValueError(
                f"Direction for rule '{output_strain}' must be a string, got {type(direction_raw).__name__}"
            )

        # Validate direction
        direction = direction_raw.lower()
        if direction not in ["alphabetical", "alpha", "rows", "numerical", "num", "columns"]:
            raise ValueError(
                f"Invalid direction '{direction}' for strain '{output_strain}'. "
                f"Must be 'alphabetical'/'alpha'/'rows' or 'numerical'/'num'/'columns'"
            )

        # Determine grouping direction
        is_alphabetical = direction in ["alphabetical", "alpha", "rows"]

        # Get unique wells and sort them appropriately
        unique_wells = strain_df[["well", "well_row", "well_col"]].drop_duplicates()
        logging.info(f"Processing strain '{output_strain}' with {len(unique_wells)} unique wells\n"
                     f"\t\t\t\t\t Direction: {direction}, use_full_direction_groups: {use_full_direction_groups}")
        if is_alphabetical:
            # Sort by row then column (A1, A2, A3, ..., B1, B2, B3, ...)
            unique_wells = unique_wells.sort_values(["well_row", "well_col"]).reset_index(drop=True)
            n_positions = len(unique_wells["well_row"].unique())
            position_name = "rows"
            grouping_col = "well_row"
            logging.info(f"Grouping wells along rows (alphabetical direction). Total unique rows: {n_positions}")
        else:
            # Sort by column then row (A1, B1, C1, ..., A2, B2, C2, ...)
            unique_wells = unique_wells.sort_values(["well_col", "well_row"]).reset_index(drop=True)
            n_positions = len(unique_wells["well_col"].unique())
            position_name = "columns"
            grouping_col = "well_col"
            logging.info(f"Grouping wells along columns (numerical direction). Total unique columns: {n_positions}")

        if use_full_direction_groups:
            # Group full rows/columns so replicate counts can vary across positions.
            unique_wells["group_num"] = pd.factorize(unique_wells[grouping_col], sort=False)[0]
        else:
            sample_size = rules.get("sample_size", 3)
            if not isinstance(sample_size, int):
                raise ValueError(
                    f"sample_size must be an integer for strain '{output_strain}', got {type(sample_size).__name__}"
                )
            if sample_size < 1:
                raise ValueError(f"sample_size must be at least 1 for strain '{output_strain}', got {sample_size}")

            # Check if divisible by sample_size
            if len(unique_wells) % sample_size != 0:
                raise ValueError(
                    f"Strain '{output_strain}': Number of wells ({len(unique_wells)}) is not divisible by "
                    f"sample_size ({sample_size}). Total {position_name}: {n_positions}"
                )

            # Create fixed-size group assignments
            unique_wells["group_num"] = unique_wells.index // sample_size

        # Merge group assignments back to strain dataframe
        strain_df = strain_df.merge(unique_wells[["well", "group_num"]], on="well", how="left")

        # Group by group_num and time, then calculate statistics
        group_cols = ["group_num"] + time_cols

        # Calculate statistics
        agg_dict = {
            value_column_name: [
                (f"{value_column_name}_mean", "mean"),
                (f"{value_column_name}_std", lambda x: x.std(ddof=ddof)),
                ("count", "count"),
            ],
            "well": [("wells", lambda x: ",".join(sorted(set(x))))],
            "well_row": [("well_rows", lambda x: ",".join(sorted(set(x))))],
            "well_col": [("well_cols", lambda x: ",".join(sorted(set(map(str, x)))))],
        }

        # Include metadata columns from first well in each group
        metadata_cols = [
            col
            for col in df_parsed.columns
            if col not in required_cols + time_cols + ["content", "group_num", "strain_group"]
        ]
        for col in metadata_cols:
            if col in strain_df.columns:
                agg_dict[col] = "first"

        stats_df = strain_df.groupby(group_cols, as_index=False).agg(agg_dict)

        # Flatten column names
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
        stats_df = stats_df.rename(columns={"count": "n_replicates"})
        stats_df["strain"] = output_strain

        stats_df["group_id"] = stats_df.apply(
            lambda row, _alpha=is_alphabetical: _make_group_id(row, is_alphabetical=_alpha), axis=1
        )

        # Drop group_num as it's internal
        if "group_num" in stats_df.columns:
            stats_df = stats_df.drop(columns=["group_num"])

        all_results.append(stats_df)
    logging.info(f"Processed {len(all_results)} strain groups with custom rules.")
    # Concatenate all results
    if not all_results:
        # Return empty DataFrame with expected columns
        logging.info("Function labUtils.media_bot.calculate_replicate_statistics_by_custom finished successfully")
        return pd.DataFrame(columns=output_columns)

    final_df = pd.concat(all_results, ignore_index=True)

    # Reorder columns for better readability
    first_cols = [
        "strain",
        "group_id",
        "wells",
        "well_rows",
        "well_cols",
        f"{value_column_name}_mean",
        f"{value_column_name}_std",
        "n_replicates",
    ]
    first_cols = [col for col in first_cols if col in final_df.columns]
    other_cols = [col for col in final_df.columns if col not in first_cols]
    final_df = final_df[first_cols + other_cols]

    # Sort by strain, group_id and time columns
    sort_cols = ["strain", "group_id"]
    for time_col in ["time_min", "time_h", "time_label"]:
        if time_col in final_df.columns:
            sort_cols.append(time_col)
            break  # Use first available time column for sorting

    final_df = final_df.sort_values(sort_cols).reset_index(drop=True)

    # Remove "_first" suffix from columns names that have it
    first_suffix_cols = [col for col in final_df.columns if col.endswith("_first")]
    for first_col in first_suffix_cols:
        final_df = final_df.rename(columns={first_col: first_col[:-6]})  # Remove '_first' suffix



    logging.info("Function labUtils.media_bot.calculate_replicate_statistics_by_custom finished successfully")
    return final_df
