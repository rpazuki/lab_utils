##############################################################
#  labUtils: A python module to run pipelines for laboratory #
#  data analysis and automation. The pipelines structure are #
#  defined in a separate yaml file.                          #
#  It will parse based on the provided file names            #
#                                                            #
#  Author: Roozbeh H. Pazuki":"2025                          #
#  License: MIT                                              #
#                                                            #
#  Examples:                                                 #
# Using mapping file (higher precedence):                    #
# growth_rate_pipeline(                                      #
#     mapping_file="file_mapping.csv"                        #
# )                                                          #
#                                                            #
# growth_rate_pipeline(                                      #
#     mapping_file="mapping.yaml",                           #
#     data_dir="Alfie/data/2025_10_27_plates_reader"         #
# )                                                          #
#                                                            #
# Using pattern-based discovery:                             #
# growth_rate_pipeline(                                      #
#     data_dir="Alfie/data",                                 #
#     raw_data_pattern="mediabot*.csv",                      #
#     meta_data_pattern="protocol_*.csv"                     #
# )                                                          #
#                                                            #
# growth_rate_pipeline(                                      #
#     data_dir="C:/path/to/data",                            #
#     pipeline_name="pipeline_2"                             #
# )                                                          #
#                                                            #
# Behavior:                                                  #
#                                                            #
# With mapping_file: Function loads exact file pairs from    #
# the mapping file                                           #
#                                                            #
# Without mapping_file: Function discovers files using       #
# patterns and pairs them alphabetically                     #
#                                                            #
# Pattern pairing: Files are sorted and paired in order      #
# (1st metadata with 1st raw data, etc.)                     #
# Mismatch handling: Warns if different numbers of files are #
# found, pairs what it can.                                  #
##############################################################

import argparse
import ast
import logging
import os
from pathlib import Path
from typing import Any, Dict, Optional, Union

import pandas as pd
import yaml

from labUtils.pipelines import build_pipeline_from_yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

log = logging.getLogger(__name__)


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


def run_media_bot_growth_rate_pipeline(
    data_dir: str = "",
    pipeline_name: str = "pipeline_1",
    mapping_file: Optional[str] = None,
    raw_data_pattern: str = "mediabot*.csv",
    meta_data_pattern: str = "protocol_metadata_*.csv",
    pipeline_yaml_path: str = str(Path(__file__).parent / "growth_rates_pipeline.yaml"),
    output_base_dir: Optional[Union[str, Path]] = None,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Process laboratory data files through a pipeline.

    Parameters
    ----------
    data_dir : str, optional
        Absolute path or relative path from DATA_DIRECTORY or current directory
        to directory containing input data files (default: "")
    pipeline_name : str, optional
        Name of the pipeline to execute (default: "pipeline_1")
    mapping_file : str, optional
        Path to file containing metadata-to-raw-data mapping (CSV, YAML, or Python dict).
        If not provided, will use pattern-based discovery (default: None)
    raw_data_pattern : str, optional
        Pattern to match raw data files (default: "mediabot*.csv")
    meta_data_pattern : str, optional
        Pattern to match metadata files (default: "protocol_metadata_*.csv")
    pipeline_yaml_path : str, optional
        Path to the pipeline YAML configuration file (default: "growth_rates_pipeline.yaml")
    output_base_dir : str or Path, optional
        Base directory for output files. If None, uses "./processed" (default: None)
    verbose : bool, optional
        Whether to print progress messages (default: True)

    Returns
    -------
    dict
        Dictionary containing processing results and statistics

    Raises
    ------
    FileNotFoundError
        If mapping file or pipeline YAML file is not found
    ValueError
        If file patterns don't match any files or invalid mapping format
    """
    # Directory containing input files
    if os.path.isabs(data_dir):
        data_dir_path = Path(data_dir)
    elif 'DATA_DIRECTORY' in globals():
        data_dir_path = Path(globals()['DATA_DIRECTORY']) / data_dir
    elif 'DATA_DIRECTORY' in locals():
        data_dir_path = Path(locals()['DATA_DIRECTORY']) / data_dir
    else:
        data_dir_path = Path(os.path.realpath(os.path.curdir)) / data_dir

    # Set output base directory
    if output_base_dir is None:
        output_base_dir = Path(os.path.realpath(os.path.curdir)) / "processed"
    else:
        output_base_dir = Path(output_base_dir)

    # Determine file mapping approach (mapping file has precedence)
    if mapping_file:
        # Load file mapping from the specified file
        try:
            file_pars = load_file_mapping(mapping_file)
            if verbose:
                print(f"Loaded {len(file_pars)} file mappings from {mapping_file}")
        except (FileNotFoundError, ValueError, IOError) as e:
            if verbose:
                print(f"Error loading mapping file: {e}")
            raise
    else:
        # Create file mapping from patterns
        try:
            file_pars = create_file_mapping_from_patterns(
                data_dir_path,
                raw_data_pattern,
                meta_data_pattern
            )
            if verbose:
                print(f"Found {len(file_pars)} file pairs using patterns:")
                print(f"  Raw data pattern: {raw_data_pattern}")
                print(f"  Metadata pattern: {meta_data_pattern}")
        except (ValueError, FileNotFoundError, OSError) as e:
            if verbose:
                print(f"Error creating file mapping from patterns: {e}")
            raise

    # Track processing results
    results = {
        'successful_files': [],
        'failed_files': [],
        'total_files': len(file_pars),
        'pipeline_outputs': {}
    }

    # Process each file through the same pipeline
    for raw_file, meta_data_file in file_pars.items():
        if verbose:
            print(f"\nProcessing: {raw_file}")

        # Create unique output directory for this file
        output_dir = Path(output_base_dir) / Path(raw_file).stem
        output_dir.mkdir(parents=True, exist_ok=True)

        # Override the raw_data input source
        # Metadata file could also be overridden if needed
        input_sources = {
            'raw_data': str(data_dir_path / raw_file),
            'meta_data': str(data_dir_path / meta_data_file)
        }

        try:
            # Build pipeline with overridden inputs
            pipeline = build_pipeline_from_yaml(
                pipeline_yaml_path,
                pipeline_name,
                output_dir=output_dir,
                input_sources=input_sources
            )

            # Execute the pipeline
            result = pipeline()

            # Track successful processing
            results['successful_files'].append(raw_file)
            results['pipeline_outputs'][raw_file] = {
                'output_dir': str(output_dir),
                'outputs': list(result.keys())
            }

            if verbose:
                print(f"  ✓ Success! Results saved to {output_dir}")
                print(f"  Generated outputs: {list(result.keys())}")

        except (ValueError, FileNotFoundError, KeyError, RuntimeError) as e:
            # Track failed processing
            results['failed_files'].append({
                'file': raw_file,
                'error': str(e)
            })

            if verbose:
                print(f"  ✗ Error: {e}")
            continue

    if verbose:
        print("\n" + "="*60)
        print("Batch processing completed!")
        print(f"Successfully processed: {len(results['successful_files'])}/{results['total_files']} files")
        if results['failed_files']:
            print(f"Failed files: {len(results['failed_files'])}")
        print("="*60)

    return results


# Example usage
if __name__ == "__main__":
    # Example 1: Using ma # results = growth_rate_pipeline(
    #     mapping_file="file_mapping.csv"
    # )

    # Example 2: Using pattern-based discovery
    # results = growth_rate_pipeline(
    #     data_dir="data",
    #     raw_data_pattern="mediabot*.csv",
    #     meta_data_pattern="protocol_*.csv"
    # )
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process laboratory data files through a pipeline')
    parser.add_argument('--data-dir',
                        type=str,
                        default="",
                        help='Absolute, relative path from DATA_DIRECTORY, or current directory path to directory containing input data files')
    parser.add_argument('--pipeline-name',
                        type=str,
                        default="pipeline_1",
                        help='Name of the pipeline to execute (default: %(default)s)')
    parser.add_argument('--mapping-file',
                        type=str,
                        default=None,
                        help='Path to file containing metadata-to-raw-data mapping (CSV, YAML, or Python dict). If not provided, will use pattern-based discovery.')

    parser.add_argument('--raw-data-pattern',
                        type=str,
                        default="mediabot*.csv",
                        help='Pattern to match raw data files (default: %(default)s)')

    parser.add_argument('--meta-data-pattern',
                        type=str,
                        default="protocol_metadata_*.csv",
                        help='Pattern to match metadata files (default: %(default)s)')

    parser.add_argument('--pipeline-yaml-path',
                        type=str,
                        default="./preprocessing_pipeline.yaml",
                        help='Path to the pipeline YAML configuration file (default: %(default)s)')

    parser.add_argument('--output-base-dir',
                        type=str,
                        default="./processed",
                        help='Base directory for output files (default: %(default)s)')

    parser.add_argument('--verbose',
                        action='store_true',
                        default=True,
                        help='Print progress messages (default: %(default)s)')

    parser.add_argument('--quiet',
                        action='store_true',
                        help='Disable verbose output (overrides --verbose)')

    args = parser.parse_args()


    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Run the growth rate pipeline
        _ = run_media_bot_growth_rate_pipeline(
            data_dir=args.data_dir,
            pipeline_name=args.pipeline_name,
            mapping_file=args.mapping_file,
            raw_data_pattern=args.raw_data_pattern,
            meta_data_pattern=args.meta_data_pattern,
            pipeline_yaml_path=args.pipeline_yaml_path,
            output_base_dir=args.output_base_dir,
            verbose=args.verbose and not args.quiet
        )

    except (ValueError, FileNotFoundError, KeyError, RuntimeError) as e:
        print(f"  ✗ Error: {e}")
