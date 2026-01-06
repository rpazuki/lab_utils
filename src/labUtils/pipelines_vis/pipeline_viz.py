"""
Pipeline visualization utilities for labUtils.

This module provides functions to visualize pipeline YAML files as diagrams.
"""

from pathlib import Path


def visualize_pipeline(
    yaml_path: str | Path, format: str = "html", output: str | Path | None = None, pipeline_name: str | None = None
) -> str:
    """
    Visualize a pipeline YAML file as a diagram.

    Parameters
    ----------
    yaml_path : str | Path
        Path to the pipeline YAML file
    format : str, optional
        Output format: 'mermaid', 'graphviz', or 'html' (default: 'html')
    output : str | Path | None, optional
        Output file path. If None, returns the diagram as a string.
    pipeline_name : str | None, optional
        Specific pipeline name to visualize. If None, uses the first pipeline.

    Returns
    -------
    str
        The generated diagram content (empty string if saved to file)

    Examples
    --------
    >>> # Generate HTML visualization
    >>> visualize_pipeline('pipeline.yaml', format='html', output='pipeline.html')

    >>> # Get Mermaid diagram as string
    >>> diagram = visualize_pipeline('pipeline.yaml', format='mermaid')
    >>> print(diagram)
    """
    from labUtils.yamls.visualize_pipeline import (  # type: ignore
        extract_pipeline_structure,
        generate_graphviz_dot,
        generate_html_interactive,
        generate_mermaid_diagram,
        parse_pipeline_yaml,
    )

    yaml_path = Path(yaml_path)
    pipeline_data = parse_pipeline_yaml(yaml_path)
    pipelines = extract_pipeline_structure(pipeline_data)

    if not pipelines:
        raise ValueError("No pipelines found in YAML file")

    # Select pipeline
    if pipeline_name:
        selected_pipeline = next((p for p in pipelines if p["name"] == pipeline_name), None)
        if not selected_pipeline:
            available = ", ".join(p["name"] for p in pipelines)
            raise ValueError(f"Pipeline '{pipeline_name}' not found. Available: {available}")
    else:
        selected_pipeline = pipelines[0]

    # Generate visualization
    if format == "mermaid":
        content = generate_mermaid_diagram(selected_pipeline)
    elif format == "graphviz":
        content = generate_graphviz_dot(selected_pipeline)
    elif format == "html":
        content = generate_html_interactive(selected_pipeline)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'mermaid', 'graphviz', or 'html'")

    # Save or return
    if output:
        output_path = Path(output)
        output_path.write_text(content, encoding="utf-8")
        return ""
    else:
        return content


def list_pipelines(yaml_path: str | Path) -> list[str]:
    """
    List all pipeline names in a YAML file.

    Parameters
    ----------
    yaml_path : str | Path
        Path to the pipeline YAML file

    Returns
    -------
    list[str]
        List of pipeline names

    Examples
    --------
    >>> pipelines = list_pipelines('growth_rates_pipeline.yaml')
    >>> print(pipelines)
    ['growth_rate_fit_pipeline', 'growth_rate_replicates_fit_pipeline']
    """
    from labUtils.yamls.visualize_pipeline import extract_pipeline_structure, parse_pipeline_yaml  # type: ignore

    yaml_path = Path(yaml_path)
    pipeline_data = parse_pipeline_yaml(yaml_path)
    pipelines = extract_pipeline_structure(pipeline_data)
    return [p["name"] for p in pipelines]
