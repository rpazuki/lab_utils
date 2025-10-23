# parse_od600_with_metadata.py
from __future__ import annotations

import re
from io import StringIO
from pathlib import Path

import pandas as pd


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

def read_bmg_export(path: Path) -> pd.DataFrame:
    """Read CLARIOstar OD600 export in wide format, return tidy long dataframe."""
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
                      var_name="time_label", value_name="od600")
    long_df["od600"] = pd.to_numeric(long_df["od600"], errors="coerce")
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

def read_meta_experiment(meta_path: Path) -> pd.DataFrame:
    """Read the '=== Experiment Data ===' section as a DataFrame and tidy column names."""
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


def parse(data_path: Path, meta_path: Path) -> pd.DataFrame:
    """Parse raw OD600 data and metadata, return merged tidy DataFrame."""
    raw_long = read_bmg_export(data_path)
    meta = read_meta_experiment(meta_path)

    df = raw_long.merge(meta, on="well", how="left", validate="m:1")

    ordered_cols = [
        "well","well_row","well_col","content","strain","is_blank","media_type","supplements",
        "media_volume_uL","water_volume_uL","supplement_volume_uL","volume_per_supplement_uL","total_volume_uL",
        "time_label","time_h","time_min_int","time_h_int","time_min","od600"
    ]
    ordered_cols = [c for c in ordered_cols if c in df.columns]
    df = df[ordered_cols].copy()
    return df

def report(data_path: Path, meta_path: Path) -> pd.DataFrame:
    """Generate a validation report for the parsed data."""
    raw_long = read_bmg_export(data_path)
    meta = read_meta_experiment(meta_path)

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
    non_num = raw_long["od600"].isna().sum()
    if non_num:
        report_rows.append({"issue": "non_numeric_od600_values", "count": int(non_num)})
    dup_meta = meta["well"].duplicated().sum()
    if dup_meta:
        report_rows.append({"issue": "duplicate_wells_in_metadata", "count": int(dup_meta)})

    return pd.DataFrame(report_rows)
