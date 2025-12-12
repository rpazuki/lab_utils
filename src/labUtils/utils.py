import ast
import re
from pathlib import Path
from typing import Literal

import numpy as np
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


def smart_join_drop_right(
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    on_cols: list[str] = ["well"],  # noqa: B006
    how: Literal["inner"] = "inner",
) -> pd.DataFrame:
    """Join two DataFrames and always drop right-side duplicate columns.

    This function merges two DataFrames and automatically drops all columns
    from the right DataFrame that have the same name as columns in the left
    DataFrame (excluding join keys). This is useful when you want to prioritize
    the left DataFrame's values for overlapping columns.

    Parameters
    ----------
    left_df : pd.DataFrame
        Left DataFrame to join (values from this DataFrame are kept for duplicate columns)
    right_df : pd.DataFrame
        Right DataFrame to join (duplicate columns from this DataFrame are dropped)
    on_cols : list[str] | str
        Column name(s) to join on. Can be a single column name or list of column names.
    how : str, optional
        Type of join to perform ('inner', 'left', 'right', 'outer'), default 'inner'

    Returns
    -------
    pd.DataFrame
        Joined DataFrame with all '_right' suffixed columns removed

    Examples
    --------
    >>> df1 = pd.DataFrame({'id': [1, 2], 'name': ['A', 'B'], 'value': [10, 20]})
    >>> df2 = pd.DataFrame({'id': [1, 2], 'name': ['C', 'D'], 'score': [100, 200]})
    >>> result = smart_join_drop_right(df1, df2, on='id')
    # Result has 'name' from df1 (['A', 'B']), 'name' from df2 is dropped
    # Result also has 'value' and 'score' columns

    >>> df1 = pd.DataFrame({'id': [1, 2], 'status': ['active', 'inactive']})
    >>> df2 = pd.DataFrame({'id': [1, 2], 'status': ['pending', 'complete']})
    >>> result = smart_join_drop_right(df1, df2, on='id')
    # Result has only 'status' with values from df1 (['active', 'inactive'])
    """
    # Ensure 'on' is a list
    if isinstance(on_cols, str):
        on_cols = [on_cols]

    # Perform the join with suffixes to identify duplicates
    merged = left_df.merge(right_df, on=on_cols, how=how, suffixes=("", "_right"))

    # Find all columns that got the '_right' suffix and drop them
    right_suffix_cols = [col for col in merged.columns if col.endswith("_right")]

    # Drop all '_right' columns
    if right_suffix_cols:
        merged = merged.drop(columns=right_suffix_cols)

    return merged


def collate_by_strain(
    folders_list: list[Path],
    csv_input_file_name: str = "growth_rates.csv",
    strain_col: str = "strain",
    create_output_folder: bool = True,
    output_dir: str | Path | None = None,
    groupby_pattern: str | None = r"[A-Za-z]+",
    csv_output_file_name: str = "growth_rates.csv",
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
    csv_output_file_name: str = "growth_rates.csv"
        Name of the output CSV file for collated data, default 'growth_rates.csv'

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
    df = pd.read_csv(folders_list[0] / csv_input_file_name)
    df["experiment"] = folders_list[0].name
    for folder in folders_list[1:]:
        if folder == output_dir:
            continue
        df_new = pd.read_csv(folder / csv_input_file_name)
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
        if create_output_folder:
            output_folder = parent_folder / str(g)
            output_folder.mkdir(exist_ok=True)
            output_file_path: Path = output_folder / csv_output_file_name
        else:
            output_file_path: Path = parent_folder / f"{g}_{csv_output_file_name}"
        group_df.to_csv(output_file_path, index=False)
        report.append(str(output_file_path))

    return pd.DataFrame({"output_file": report})


# =====================================================================
#
#      Moleculear conversions
#
# =====================================================================

COMPOUND_CACHE_BY_NAME = {}
COMPOUND_CACHE_BY_CID = {}
COMPOUND_CACHE_BY_SMILES = {}


def get_compound_by_name(compound_name: str, chache: bool = True):
    try:
        from pubchempy import get_compounds

        if chache and compound_name in COMPOUND_CACHE_BY_NAME:
            return COMPOUND_CACHE_BY_NAME[compound_name]

        compound = get_compounds(compound_name, "name")
        if chache and len(compound) > 0:
            COMPOUND_CACHE_BY_NAME[compound_name] = compound[0]

        if len(compound) > 0:
            return compound[0]
        else:
            print(
                f"Warning: pubchempy did not find the compound name:'{compound_name}', cannot fetch molecular weight."
            )
            return None
    except ImportError as e:
        print("Warning: pubchempy not available.")
        print(f"ImportError: {e}")
        return None
    except Exception as e:
        print(f"Warning: An error occurred while fetching molecular weight for name:'{compound_name}': {e}")
        return None


def get_compound_by_id(compound_id: str, chache: bool = True):
    try:
        from pubchempy import get_compounds

        if chache and compound_id in COMPOUND_CACHE_BY_CID:
            return COMPOUND_CACHE_BY_CID[compound_id]

        compound = get_compounds(compound_id, "cid")
        if chache and len(compound) > 0:
            COMPOUND_CACHE_BY_CID[compound_id] = compound[0]
        if len(compound) > 0:
            return compound[0]
        else:
            print(f"Warning: pubchempy did not find the compound cid:'{compound_id}', cannot fetch molecular weight.")
            return None
    except ImportError as e:
        print("Warning: pubchempy not available.")
        print(f"ImportError: {e}")
        return None
    except Exception as e:
        print(f"Warning: An error occurred while fetching molecular weight for cid:'{compound_id}': {e}")
        return None


def get_compound_by_smiles(smiles: str, chache: bool = True):
    try:
        from pubchempy import get_compounds

        if chache and smiles in COMPOUND_CACHE_BY_SMILES:
            return COMPOUND_CACHE_BY_SMILES[smiles]

        compound = get_compounds(smiles, "smiles")
        if chache and len(compound) > 0:
            COMPOUND_CACHE_BY_SMILES[smiles] = compound[0]
        if len(compound) > 0:
            return compound[0]
        else:
            print(f"Warning: pubchempy did not find the compound smiles:'{smiles}', cannot fetch molecular weight.")
            return None
    except ImportError as e:
        print("Warning: pubchempy not available.")
        print(f"ImportError: {e}")
        return None
    except Exception as e:
        print(f"Warning: An error occurred while fetching molecular weight for smiles:'{smiles}': {e}")
        return None


# Minimal periodic table (standard atomic weights)
ATOMIC_WEIGHTS = {
    "H": 1.00794,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3050,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973761,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955912,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938045,
    "Fe": 55.845,
    "Co": 58.933195,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.92160,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.96,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.9054519,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.242,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.500,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.054,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98040,
    "Po": 209.0,
    "At": 210.0,
    "Rn": 222.0,
}

EL_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
ION_NAME = re.compile(r"\[([^\]]+?)\]")


def parse_ion_formula(component: str) -> tuple[dict[str, int], int]:
    """
    Parse a SMILES component like [Na+], [Cl-], [NH4+], [Ca++], [Fe(CN)6]4- into:
      - composition dict {element: count}
      - net charge (int)
    Heuristic: strips brackets, charges, parentheses; counts element tokens.
    """
    s = component.strip()
    # Extract trailing charge like 4- or 3+
    m = re.search(r"(\d+)([+-])$", s)
    if m:
        charge = int(m.group(1)) * (1 if m.group(2) == "+" else -1)
        s = re.sub(r"(\d+)([+-])$", "", s)
    else:
        # Count +/- signs inside component
        charge = s.count("+") - s.count("-")

    # Remove brackets and common non-element symbols
    s = s.replace("[", "").replace("]", "")
    s = re.sub(r"[+\-@.]", "", s)
    s = s.replace("(", "").replace(")", "")

    comp: dict[str, int] = {}
    for elem, count in EL_TOKEN.findall(s):
        n = int(count) if count else 1
        comp[elem] = comp.get(elem, 0) + n

    return comp, charge


def molar_mass(composition: dict[str, int]) -> float:
    mass = 0.0
    for elem, n in composition.items():
        if elem not in ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown element '{elem}' in composition")
        mass += ATOMIC_WEIGHTS[elem] * n
    return mass


def ions_from_smiles(smiles: str) -> list[dict]:
    """Split a salt SMILES at '.' and parse components into ions with masses."""

    parts = [p for p in smiles.split(".") if p]
    ions = []
    for p in parts:
        comp, ch = parse_ion_formula(p)
        mass = molar_mass(comp)

        # Try to get compound info using the original SMILES component
        p_cid = get_compound_by_smiles(p)

        # Determine the simplified component for secondary lookup
        # Only simplify if it's a simple ion in brackets like [Na+], [Cl-], [NH4+]
        # For complex structures like OP(=O)([O-])[O-], keep the original
        m = ION_NAME.search(p)

        # Check if the entire component is just a simple bracketed ion
        # (starts with '[' and ends with ']' or charge notation)
        is_simple_ion = p.startswith("[") and re.search(r"\][0-9]*[+-]*$", p)

        if m and is_simple_ion:
            # Extract the content inside brackets and remove charge symbols
            inner = m.group(1)
            # Remove trailing + or - signs
            inner_clean = re.sub(r"[+\-]+$", "", inner)
            smiles_component = inner_clean
            p_cid_secondary = get_compound_by_name(smiles_component)
        else:
            # For complex structures, keep the original SMILES
            smiles_component = p
            p_cid_secondary = None

        ions.append(
            {
                "smiles_component": smiles_component,
                "iupac_name": p_cid.iupac_name if p_cid is not None else "Unknown",
                "iupac_name_secondary": p_cid_secondary.iupac_name if p_cid_secondary is not None else None,
                "pubchem_cid": p_cid.cid if p_cid is not None else None,
                "composition": comp,
                "charge": ch,
                "molar_mass": mass,
            }
        )
    return ions


def mass_fractions(ions: list[dict]) -> list[dict]:
    total = sum(i["molar_mass"] for i in ions)
    return [{**i, "mass_fraction": (i["molar_mass"] / total) if total > 0 else None} for i in ions]


# --- Examples ---
# NaCl
# print(mass_fractions(ions_from_smiles("[Na+].[Cl-]")))
# NH4Cl
# print(mass_fractions(ions_from_smiles("[NH4+].[Cl-]")))
# CaCl2 (one common ionic representation)
# print(mass_fractions(ions_from_smiles("[Ca++].[Cl-].[Cl-]")))
# import pubchempy as pcp
# pcp.get_compounds("[Na+]", "smiles")
# get_compound_by_smiles("[Na+].[Cl-]")
# ions_from_smiles("OP(=O)([O-])[O-].[Na+].[Na+]")
# comp = get_compound_by_smiles("OP(=O)([O-])[O-]", False)
# print(comp)


def od600_to_gCDW(od600: float | list[float] | np.ndarray, conversion_rate: float = 0.4) -> float | np.ndarray:  # noqa: N802
    """
    Convert OD600 to grams of cell dry weight (gCDW) for E. coli.

    Parameters
    ----------
    od600 : float | list[float] | np.ndarray
        Optical density at 600 nm (OD600) measurement

    conversion_rate : float
        Conversion factor from OD600 to gCDW per litre (default: 0.4 gCDW/L per OD600)
    """
    if isinstance(od600, (list, np.ndarray)):
        od600 = np.array(od600)  # Convert to numpy array for vectorized operations

    return od600 * conversion_rate


def mg_to_mmol(
    mass_mg: float | list[float] | np.ndarray,
    molar_mass_g_per_mol: float,
) -> float | np.ndarray:
    """
    Convert mass in miligrams to millimoles.

    Parameters
    ----------
    mass_mg : float | list[float] | np.ndarray
        Mass in miligrams

    molar_mass_g_per_mol : float
        Molar mass of the substance in grams per mole

    Returns
    -------
    float | np.ndarray
        Amount in millimoles
    """
    if isinstance(mass_mg, (list, np.ndarray)):
        mass_mg = np.array(mass_mg)  # Convert to numpy array for vectorized operations

    return mass_mg / molar_mass_g_per_mol  # mg to mmol


def find_molecular_weight(
    compound_name: str,
) -> float | None:
    """
    Find the molecular weight of a compound using PubChem.

    Parameters
    ----------
    compound_name : str
        Common name of the compound

    Returns
    -------
    float | None
        Molecular weight in g/mol, or None if not found

    Examples
    --------
    >>> mw = find_molecular_weight("glucose" )
    >>> print(mw)  # 180.16
    """
    compound = get_compound_by_name(compound_name)
    if compound is not None:
        return compound.molecular_weight
    else:
        return None


def find_molecular_weight_by_id(
    compound_id: str,
) -> float | None:
    """
    Find the molecular weight of a compound using PubChem.

    Parameters
    ----------
    compound_id : str
        PubChem ID of the compound
    Returns
    -------
    float | None
        Molecular weight in g/mol, or None if not found

    Examples
    --------
    >>> mw = find_molecular_weight_by_id("5793" )
    >>> print(mw)  # 180.16
    """
    compound = get_compound_by_id(compound_id)
    if compound is not None:
        return compound.molecular_weight
    else:
        return None
