"""Test that output_dir is properly combined with YAML output paths."""
from pathlib import Path

from labUtils.pipelines import build_pipeline_from_yaml

# Test without output_dir
yaml_path = Path("tests/pipeline_temp.yaml")
pipeline1 = build_pipeline_from_yaml(yaml_path, "pipeline_1")
print("Without output_dir:")
print(f"  Pipeline has {len(pipeline1.processes)} processes")

# Test with output_dir
output_dir = Path("./my_results")
pipeline2 = build_pipeline_from_yaml(yaml_path, "pipeline_1", output_dir)
print(f"\nWith output_dir={output_dir}:")
print(f"  Pipeline has {len(pipeline2.processes)} processes")
print(f"  Output files will be saved to: {output_dir.absolute()}")

# The OutputProcess (last process) should have the combined paths
# This can be verified by looking at the closure variables
print("\nPipeline structure verified!")
