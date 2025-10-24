#!/usr/bin/env python
"""Example: Process a single file with input override."""
from pathlib import Path

from labUtils.pipelines import build_pipeline_from_yaml

# Example 1: Use default sources from YAML
print("Example 1: Using default sources from YAML")
pipeline = build_pipeline_from_yaml(
    'src/labUtils/pipeline_temp.yaml',
    'pipeline_1'
)
print(f"✓ Pipeline created with default inputs\n")

# Example 2: Override one input source
print("Example 2: Override raw_data input")
pipeline = build_pipeline_from_yaml(
    'src/labUtils/pipeline_temp.yaml',
    'pipeline_1',
    input_sources={'raw_data': 'tests/alternative_data.csv'}
)
print(f"✓ Pipeline created with raw_data override\n")

# Example 3: Override multiple inputs and set output directory
print("Example 3: Override multiple inputs + output directory")
pipeline = build_pipeline_from_yaml(
    'src/labUtils/pipeline_temp.yaml',
    'pipeline_1',
    output_dir='results/experiment_001',
    input_sources={
        'raw_data': 'data/experiment_001_raw.csv',
        'meta_data': 'data/experiment_001_meta.csv'
    }
)
print(f"✓ Pipeline created with multiple overrides\n")

print("="*60)
print("Command-line usage examples:")
print("="*60)
print()
print("# Use default inputs from YAML:")
print("python run_pipeline.py config.yaml pipeline_1")
print()
print("# Override one input:")
print("python run_pipeline.py config.yaml pipeline_1 \\")
print("  -i raw_data=data/file1.csv")
print()
print("# Override multiple inputs with custom output directory:")
print("python run_pipeline.py config.yaml pipeline_1 \\")
print("  -i raw_data=data/exp1.csv \\")
print("  -i meta_data=data/exp1_meta.csv \\")
print("  -o results/exp1")
print()
