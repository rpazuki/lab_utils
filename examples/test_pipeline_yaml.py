"""Example script to test YAML-based pipeline construction."""
from pathlib import Path

from labUtils.utils import build_pipeline_from_yaml


def test_pipeline_yaml_builds_pipeline():
    yaml_path = Path(__file__).with_name("pipeline_temp.yaml")
    pipeline, _ = build_pipeline_from_yaml(yaml_path, "pipeline_1")

    assert len(pipeline.processes) > 0
