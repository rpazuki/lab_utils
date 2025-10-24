#!/usr/bin/env python
"""Example: Process multiple files through the same pipeline."""
from pathlib import Path

from labUtils.pipelines import build_pipeline_from_yaml

# Directory containing input files
data_dir = Path("data")
output_base_dir = Path("results")

# Find all CSV files that match a pattern
raw_data_files = list(data_dir.glob("experiment_*.csv"))

# Configuration
yaml_config = Path("src/labUtils/pipeline_temp.yaml")
pipeline_name = "pipeline_1"

print(f"Found {len(raw_data_files)} files to process")

# Process each file through the same pipeline
for raw_file in raw_data_files:
    print(f"\nProcessing: {raw_file.name}")

    # Create unique output directory for this file
    output_dir = output_base_dir / raw_file.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    # Override the raw_data input source
    # Metadata file could also be overridden if needed
    input_sources = {
        'raw_data': str(raw_file),
        # 'meta_data': str(data_dir / f"{raw_file.stem}_metadata.csv")
    }

    try:
        # Build pipeline with overridden inputs
        pipeline = build_pipeline_from_yaml(
            yaml_config,
            pipeline_name,
            output_dir=output_dir,
            input_sources=input_sources
        )

        # Execute the pipeline
        result = pipeline()

        print(f"  ✓ Success! Results saved to {output_dir}")
        print(f"  Generated outputs: {list(result.keys())}")

    except Exception as e:
        print(f"  ✗ Error: {e}")
        continue

print("\n" + "="*60)
print("Batch processing completed!")
print("="*60)
