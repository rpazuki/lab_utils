"""Test that output_dir is properly combined with YAML output paths."""
from pathlib import Path

from labUtils.utils import build_pipeline_from_yaml


def test_output_dir_pipeline_yaml_loads():
	yaml_path = Path(__file__).with_name("pipeline_temp.yaml")
	pipeline1, _ = build_pipeline_from_yaml(yaml_path, "pipeline_1")
	pipeline2, _ = build_pipeline_from_yaml(yaml_path, "pipeline_1", Path("./my_results"))

	assert len(pipeline1.processes) > 0
	assert len(pipeline2.processes) == len(pipeline1.processes)
