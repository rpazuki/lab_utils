import ast
from pathlib import Path
from typing import Literal

import pandas as pd
import yaml

from labUtils.pipelines import DFPipeline, build_pipeline_from_yaml_string


def list_folders(path: Path) -> list[Path]:
    """List all folders in the given directory path."""
    path = Path(path)
    if not path.exists() or not path.is_dir():
        raise ValueError(f"Provided path is not a valid directory: {path}")

    folders = [item for item in path.iterdir() if item.is_dir()]
    return folders


def read_csv(path: Path) -> pd.DataFrame:
    """Utility function to read a CSV file into a DataFrame."""
    df = pd.read_csv(path)
    return df


def load_file_mapping(file_path):
    """
    Load file mapping from various formats (CSV, Python dict, YAML).

    Args:
        file_path (str): Path to the mapping file

    Returns:
        dict: Dictionary mapping metadata files to raw data files
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Mapping file not found: {file_path}")

    # Determine file format based on extension
    if file_path.suffix.lower() == ".csv":
        # Load CSV with two columns
        df = pd.read_csv(file_path)
        if len(df.columns) != 2:
            raise ValueError("CSV file must have exactly 2 columns (metadata_file, raw_data_file)")

        # Convert to dictionary using first column as key, second as value
        return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

    elif file_path.suffix.lower() in [".yaml", ".yml"]:
        # Load YAML file
        with open(file_path, encoding="utf-8") as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError("YAML file must contain a dictionary")

        return data

    elif file_path.suffix.lower() == ".py":
        # Load Python dictionary from .py file
        with open(file_path, encoding="utf-8") as f:
            content = f.read()

        # Try to extract dictionary from the file
        # Look for variable assignment like: file_pars = {...}
        try:
            # Parse the file content as Python code
            tree = ast.parse(content)
            for node in tree.body:
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id in ["file_pars", "file_mapping", "mapping"]:
                            return ast.literal_eval(node.value)

            # If no variable found, try to evaluate the entire content as a dict
            return ast.literal_eval(content.strip())

        except (ValueError, SyntaxError) as e:
            raise ValueError(f"Could not parse Python dictionary from file: {e}") from e

    else:
        # Try to parse as JSON/dict literal
        try:
            with open(file_path, encoding="utf-8") as f:
                content = f.read().strip()
            return ast.literal_eval(content)
        except (ValueError, SyntaxError) as e:
            raise ValueError(f"Unsupported file format or invalid content: {e}") from e


def create_file_mapping_from_patterns(data_dir, raw_pattern, meta_pattern):
    """
    Create file mapping by finding files matching patterns and pairing them.

    Args:
        data_dir (Path): Directory to search for files
        raw_pattern (str): Pattern to match raw data files
        meta_pattern (str): Pattern to match metadata files

    Returns:
        dict: Dictionary mapping metadata files to raw data files
    """
    data_dir = Path(data_dir)

    # Find all files matching patterns
    raw_data_files = list(data_dir.glob(raw_pattern))
    meta_data_files = list(data_dir.glob(meta_pattern))

    if not raw_data_files:
        raise ValueError(f"No raw data files found matching pattern: {raw_pattern}")

    if not meta_data_files:
        raise ValueError(f"No metadata files found matching pattern: {meta_pattern}")

    # Sort files to ensure consistent pairing
    raw_data_files.sort()
    meta_data_files.sort()

    if len(raw_data_files) != len(meta_data_files):
        print(f"Warning: Found {len(raw_data_files)} raw data files but {len(meta_data_files)} metadata files")
        print("Files will be paired in order, remaining files will be skipped")

    # Create mapping by pairing files (metadata -> raw data)
    file_mapping = {}
    min_length = min(len(raw_data_files), len(meta_data_files))

    for i in range(min_length):
        meta_file = meta_data_files[i].name
        raw_file = raw_data_files[i].name
        file_mapping[raw_file] = meta_file

    return file_mapping


def build_pipeline_from_yaml(
    yaml_path: str | Path,
    pipeline_name: str,
    output_dir: str | Path | None = None,
    input_sources: dict[str, str] | None = None,
    process_arg_mapping: dict[str, dict[str, str]] | None = None,
) -> tuple[DFPipeline, dict]:
    """Build a DFPipeline from a YAML configuration file.

    Parameters
    ----------
    yaml_path : str | Path
        Path to the YAML configuration file
    pipeline_name : str
        Name of the pipeline to build from the YAML file
    output_dir : str | Path, optional
        Output directory to prepend to all output file paths from YAML.
        If None, uses paths as specified in YAML.
    input_sources : dict[str, str], optional
        Dictionary mapping input names to source paths. This overrides the 'src'
        field in the YAML for specified inputs. Allows reusing the same pipeline
        configuration with different input files.
        Example: {'raw_data': 'file1.csv', 'meta_data': 'metadata1.csv'}
    process_arg_mapping : dict[str, dict[str, str]], optional
        Dictionary mapping process names to argument mappings. This allows
        overriding or remapping arguments for specific processes.

    Returns
    -------
    DFPipeline
        A configured DFPipeline ready to execute

    Raises
    ------
    ValueError
        If the YAML structure is invalid or pipeline_name is not found
    FileNotFoundError
        If the YAML file does not exist

    Examples
    --------
    # Use default sources from YAML
    pipeline = build_pipeline_from_yaml('config.yaml', 'pipeline_1')

    # Override input sources
    pipeline = build_pipeline_from_yaml(
        'config.yaml',
        'pipeline_1',
        input_sources={'raw_data': 'data/file1.csv'}
    )

    # Process multiple files with same pipeline
    for file in data_files:
        pipeline = build_pipeline_from_yaml(
            'config.yaml',
            'pipeline_1',
            input_sources={'raw_data': file},
            output_dir=f'results/{file.stem}'
        )
        result = pipeline()
    """
    # Load YAML file
    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")

    with open(yaml_path, encoding="utf-8") as f:
        yaml_string = f.read()

    # Delegate to build_pipeline_from_yaml_string
    return build_pipeline_from_yaml_string(
        yaml_string=yaml_string,
        pipeline_name=pipeline_name,
        output_dir=output_dir,
        input_sources=input_sources,
        process_arg_mapping=process_arg_mapping,
    )


def build_pipeline_from_lib_yaml(
    lib_pipelines_yaml_name: str,
    pipeline_name: str,
    output_dir: str | Path | None = None,
    input_sources: dict[str, str] | None = None,
    process_arg_mapping: dict[str, dict[str, str]] | None = None,
) -> tuple[DFPipeline, dict]:
    """Build a data processing pipeline from a library YAML file."""
    yaml_path = Path(__file__).parent / "yamls" / lib_pipelines_yaml_name
    return build_pipeline_from_yaml(
        yaml_path,
        pipeline_name,
        output_dir=output_dir,
        input_sources=input_sources,
        process_arg_mapping=process_arg_mapping,
    )


def smart_join(
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    on_cols: list[str] = ["well"],  # noqa: B006
    how: Literal["inner"] = "inner",
) -> pd.DataFrame:
    """Join two DataFrames with intelligent handling of duplicate columns.

    For columns with the same name in both DataFrames (excluding join keys):
    - If values are identical, keeps only one copy
    - If values differ, keeps both with right column renamed as 'column_name_right'

    Parameters
    ----------
    left_df : pd.DataFrame
        Left DataFrame to join
    right_df : pd.DataFrame
        Right DataFrame to join
    on_cols : list[str] | str
        Column name(s) to join on. Can be a single column name or list of column names.
    how : str, optional
        Type of join to perform ('inner', 'left', 'right', 'outer'), default 'inner'

    Returns
    -------
    pd.DataFrame
        Joined DataFrame with intelligent duplicate column handling

    Examples
    --------
    >>> df1 = pd.DataFrame({'id': [1, 2], 'name': ['A', 'B'], 'value': [10, 20]})
    >>> df2 = pd.DataFrame({'id': [1, 2], 'name': ['A', 'B'], 'score': [100, 200]})
    >>> result = smart_join(df1, df2, on='id')
    # 'name' column appears once since values are identical

    >>> df1 = pd.DataFrame({'id': [1, 2], 'status': ['active', 'inactive']})
    >>> df2 = pd.DataFrame({'id': [1, 2], 'status': ['pending', 'complete']})
    >>> result = smart_join(df1, df2, on='id')
    # Result has 'status' and 'status_right' since values differ
    """
    # Ensure 'on' is a list
    if isinstance(on_cols, str):
        on_cols = [on_cols]

    # Perform the join with suffixes to identify duplicates
    merged = left_df.merge(right_df, on=on_cols, how=how, suffixes=("", "_right"))

    # Find columns that got the '_right' suffix
    right_suffix_cols = [col for col in merged.columns if col.endswith("_right")]

    # Process each column with '_right' suffix
    cols_to_drop = []
    for right_col in right_suffix_cols:
        # Get the original column name (without '_right')
        original_col = right_col[:-6]  # Remove '_right' suffix

        # Check if both columns exist (they should)
        if original_col in merged.columns:
            # Compare values (handling NaN equality)
            left_vals = merged[original_col]
            right_vals = merged[right_col]

            # Check if columns are identical (considering NaN == NaN as True)
            are_identical = ((left_vals == right_vals) | (left_vals.isna() & right_vals.isna())).all()

            if are_identical:
                # Values are identical, drop the '_right' version
                cols_to_drop.append(right_col)

    # Drop columns where values were identical
    if cols_to_drop:
        merged = merged.drop(columns=cols_to_drop)

    return merged


def collate_by_strain(
    folders_list: list[Path],
    csv_file_name: str = "growth_rates.csv",
    strain_col: str = "strain",
    create_output_folder: bool = True,
    output_dir: str | Path | None = None,
    groupby_pattern: str | None = r"[A-Za-z]+",
) -> pd.DataFrame:
    """Collate data by strain, averaging values across wells for each time point.

    Parameters
    ----------
    folders_list : list[Path]
        List of folder paths containing parsed data files
    csv_file_name : str, optional
        Name of the CSV file to read from each folder, default 'growth_rates.csv'
    strain_col : str, optional
        Column name for strain identifiers, default 'strain'
    create_output_folder : bool, optional
        Whether to create an output folder for collated files, default True
    output_dir : str | Path, optional
        Directory to save collated files. If None, uses parent of first folder in folders_list
    groupby_pattern : str | None, optional
        Regex pattern to extract grouping key from strain_col values.
        Default is "[A-Za-z]+" (one or more letters - extracts alphabetical part only).
        Example: "strainA1", "strainA2" -> grouped as "strainA"
        If None, uses strain_col values as-is without pattern extraction.

    Returns
    -------
    pd.DataFrame
        Collated DataFrame with mean values per strain and time point
    """
    assert len(folders_list) > 0, "folders_list must contain at least one folder path"
    if output_dir is not None:
        parent_folder = Path(output_dir)
        parent_folder.mkdir(exist_ok=True)
    else:
        parent_folder: Path = folders_list[0].parent
    df = pd.read_csv(folders_list[0] / csv_file_name)
    df["experiment"] = folders_list[0].name
    for folder in folders_list[1:]:
        if folder == output_dir:
            continue
        df_new = pd.read_csv(folder / csv_file_name)
        df_new["experiment"] = folder.name
        df = pd.concat([df, df_new], ignore_index=True)
    # Apply pattern extraction if specified
    if groupby_pattern is not None:
        df["_groupby_key"] = df[strain_col].astype(str).str.extract(f"({groupby_pattern})", expand=False)
        groupby_col = "_groupby_key"
    else:
        groupby_col = strain_col

    report = []
    grouped = df.groupby([groupby_col])
    for g in grouped.groups:
        group_df = grouped.get_group((g,)).copy()
        # Drop the temporary groupby key if it was created
        if groupby_pattern is not None:
            group_df = group_df.drop(columns=["_groupby_key"])
        output_file_name = f"{g}_growth_rates_collated.csv"
        if create_output_folder:
            output_folder = parent_folder / str(g)
            output_folder.mkdir(exist_ok=True)
            output_file_name = output_folder / output_file_name
        else:
            output_file_name = parent_folder / output_file_name
        group_df.to_csv(output_file_name, index=False)
        report.append(str(output_file_name))

    return pd.DataFrame({"output_file": report})
