#!/usr/bin/env python
"""Run a pipeline from a YAML configuration file."""
import argparse
import logging
import sys
from pathlib import Path

from labUtils.pipelines import build_pipeline_from_yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
log = logging.getLogger(__name__)


def main():
    """Main function to run a pipeline from YAML configuration.

    Accepts command line arguments for:
    - YAML configuration file path
    - Pipeline name to execute
    - Output directory for results
    """
    parser = argparse.ArgumentParser(
        description='Run a data processing pipeline from a YAML configuration file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python run_pipeline.py config.yaml pipeline_1 --output-dir ./results
  python run_pipeline.py config.yaml pipeline_1 -o ./output -v
        """
    )

    parser.add_argument(
        'yaml_file',
        type=str,
        help='Path to the YAML configuration file'
    )

    parser.add_argument(
        'pipeline_name',
        type=str,
        help='Name of the pipeline to execute from the YAML file'
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        default='./output',
        help='Output directory for results (default: ./output)'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate YAML file exists
    yaml_path = Path(args.yaml_file)
    if not yaml_path.exists():
        log.error(f"YAML file not found: {yaml_path}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Output directory: {output_dir.absolute()}")

    try:
        # Build the pipeline from YAML
        log.info(f"Loading pipeline '{args.pipeline_name}' from {yaml_path}")
        pipeline = build_pipeline_from_yaml(yaml_path, args.pipeline_name, output_dir)
        log.info(f"Pipeline created successfully with {len(pipeline.processes)} processes")

        # Execute the pipeline
        log.info("Executing pipeline...")
        result = pipeline()

        log.info("Pipeline execution completed successfully")
        log.debug(f"Result keys: {list(result.keys())}")

        # Print summary
        print("\n" + "="*60)
        print("Pipeline Execution Summary")
        print("="*60)
        print(f"Pipeline: {args.pipeline_name}")
        print(f"Output directory: {output_dir.absolute()}")
        print(f"Result payload contains {len(result)} items:")
        for key in result.keys():
            value = result[key]
            value_type = type(value).__name__
            print(f"  - {key}: {value_type}")
        print("="*60)

        return 0

    except FileNotFoundError as e:
        log.error(f"File not found: {e}")
        sys.exit(1)
    except ValueError as e:
        log.error(f"Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        log.error(f"Error executing pipeline: {e}", exc_info=args.verbose)
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
