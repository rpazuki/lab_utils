"""Example script to test YAML-based pipeline construction."""
from pathlib import Path

from labUtils.pipelines import build_pipeline_from_yaml

# Build the pipeline from YAML
yaml_path = Path("src/labUtils/pipeline_temp.yaml")
pipeline = build_pipeline_from_yaml(yaml_path, "pipeline_1")

# Execute the pipeline
result = pipeline()

print(f"Pipeline created successfully with {len(pipeline.processes)} processes")
print("Process types:")
for i, proc in enumerate(pipeline.processes):
    print(f"  {i+1}. {type(proc).__name__}")
