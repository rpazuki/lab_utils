import ast
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml

from labUtils.pipelines import DFPipeline, build_pipeline_from_yaml_string


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
    if file_path.suffix.lower() == '.csv':
        # Load CSV with two columns
        df = pd.read_csv(file_path)
        if len(df.columns) != 2:
            raise ValueError("CSV file must have exactly 2 columns (metadata_file, raw_data_file)")

        # Convert to dictionary using first column as key, second as value
        return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

    elif file_path.suffix.lower() in ['.yaml', '.yml']:
        # Load YAML file
        with open(file_path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError("YAML file must contain a dictionary")

        return data

    elif file_path.suffix.lower() == '.py':
        # Load Python dictionary from .py file
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Try to extract dictionary from the file
        # Look for variable assignment like: file_pars = {...}
        try:
            # Parse the file content as Python code
            tree = ast.parse(content)
            for node in tree.body:
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id in ['file_pars', 'file_mapping', 'mapping']:
                            return ast.literal_eval(node.value)

            # If no variable found, try to evaluate the entire content as a dict
            return ast.literal_eval(content.strip())

        except (ValueError, SyntaxError) as e:
            raise ValueError(f"Could not parse Python dictionary from file: {e}") from e

    else:
        # Try to parse as JSON/dict literal
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
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
    output_dir: Optional[str | Path] = None,
    input_sources: Optional[dict[str, str]] = None
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

    with open(yaml_path, 'r', encoding='utf-8') as f:
        yaml_string = f.read()

    # Delegate to build_pipeline_from_yaml_string
    return build_pipeline_from_yaml_string(
        yaml_string=yaml_string,
        pipeline_name=pipeline_name,
        output_dir=output_dir,
        input_sources=input_sources
    )

def build_pipeline_from_lib_yaml(
    lib_pipelines_yaml_name: str,
    pipeline_name: str,
    output_dir: Optional[str | Path] = None,
    input_sources: Optional[dict[str, str]] = None
) -> tuple[DFPipeline, dict]:
    """Build a data processing pipeline from a library YAML file."""
    yaml_path = Path(__file__).parent / "yamls" / lib_pipelines_yaml_name
    return build_pipeline_from_yaml(
        yaml_path,
        pipeline_name,
        output_dir=output_dir,
        input_sources=input_sources
    )
