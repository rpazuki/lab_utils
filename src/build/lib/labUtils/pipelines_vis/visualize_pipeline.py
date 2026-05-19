"""
Pipeline Visualization Script for labUtils

This script reads YAML pipeline configurations and generates visual diagrams
in multiple formats (Mermaid, Graphviz, HTML).

Usage:
    python visualize_pipeline.py <yaml_file> [--format mermaid|graphviz|html] [--output <file>]

Example:
    python visualize_pipeline.py ../yamls/growth_rates_pipeline.yaml --format mermaid --output pipeline.md
"""

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml


def parse_pipeline_yaml(yaml_path: Path) -> dict[str, Any]:
    """Parse a pipeline YAML file and return its structure."""
    with open(yaml_path, encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data


def extract_pipeline_structure(pipeline_data: dict) -> list[dict]:
    """Extract structured pipeline information from YAML data."""
    pipelines = []

    if "pipelines" not in pipeline_data:
        return pipelines

    for pipeline_def in pipeline_data["pipelines"]:
        for pipeline_name, pipeline_config in pipeline_def.items():
            pipeline_info = {
                "name": pipeline_name,
                "inputs": [],
                "processes": [],
                "outputs": [],
            }

            # Extract inputs
            if "Inputs" in pipeline_config:
                for input_def in pipeline_config["Inputs"]:
                    for input_name, input_config in input_def.items():
                        input_info = {"name": input_name, "details": {}}
                        if isinstance(input_config, list):
                            for item in input_config:
                                if isinstance(item, dict):
                                    input_info["details"].update(item)
                        pipeline_info["inputs"].append(input_info)

            # Extract processes
            if "Processes" in pipeline_config:
                for process_def in pipeline_config["Processes"]:
                    for process_name, process_config in process_def.items():
                        process_info = {
                            "name": process_name,
                            "package": process_config.get("package", ""),
                            "method": process_config.get("method", ""),
                            "parameters": process_config.get("parameters", {}),
                        }
                        pipeline_info["processes"].append(process_info)

            # Extract outputs
            if "Outputs" in pipeline_config:
                for output_def in pipeline_config["Outputs"]:
                    for output_name, output_file in output_def.items():
                        pipeline_info["outputs"].append({"name": output_name, "file": output_file})

            pipelines.append(pipeline_info)

    return pipelines


def generate_mermaid_diagram(pipeline: dict) -> str:
    """Generate a Mermaid flowchart from pipeline structure."""
    lines = ["```mermaid", "graph TD"]

    # Add inputs
    for inp in pipeline["inputs"]:
        input_id = inp["name"]
        method = inp["details"].get("method", "")
        src = inp["details"].get("src", "")
        label = f"{input_id}"
        if method:
            label += f"<br/>{method}"
        lines.append(f'    {input_id}[("{label}")]:::input')

    # Add processes
    for i, proc in enumerate(pipeline["processes"]):
        proc_id = proc["name"]
        method = proc["method"]
        label = f"{proc_id}<br/>{method}"
        lines.append(f'    {proc_id}["{label}"]:::process')

        # Connect to dependencies
        params = proc["parameters"]
        for param_name, param_value in params.items():
            if isinstance(param_value, str):
                # Check if this references another node
                if any(inp["name"] == param_value for inp in pipeline["inputs"]) or any(
                    p["name"] == param_value for p in pipeline["processes"][:i]
                ):
                    lines.append(f"    {param_value} --> {proc_id}")

    # Add outputs
    for out in pipeline["outputs"]:
        output_id = f"{out['name']}_out"
        label = f"{out['file']}"
        lines.append(f'    {output_id}[("{label}")]:::output')
        lines.append(f"    {out['name']} --> {output_id}")

    # Add styling
    lines.extend(
        [
            "",
            "    classDef input fill:#e1f5ff,stroke:#01579b,stroke-width:2px",
            "    classDef process fill:#fff9c4,stroke:#f57f17,stroke-width:2px",
            "    classDef output fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px",
            "```",
        ]
    )

    return "\n".join(lines)


def generate_graphviz_dot(pipeline: dict) -> str:
    """Generate a Graphviz DOT file from pipeline structure."""
    lines = ["digraph Pipeline {", "    rankdir=LR;", "    node [shape=box, style=rounded];", ""]

    # Define node styles
    lines.extend(
        [
            '    node [fillcolor="#e1f5ff", style="filled,rounded"];',
            "    // Inputs",
        ]
    )

    # Add inputs
    for inp in pipeline["inputs"]:
        input_id = inp["name"]
        method = inp["details"].get("method", "")
        label = f"{input_id}"
        if method:
            label += f"\\n{method}"
        lines.append(f'    "{input_id}" [label="{label}"];')

    lines.extend(
        [
            "",
            '    node [fillcolor="#fff9c4", style="filled,rounded"];',
            "    // Processes",
        ]
    )

    # Add processes
    for i, proc in enumerate(pipeline["processes"]):
        proc_id = proc["name"]
        method = proc["method"]
        label = f"{proc_id}\\n{method}"
        lines.append(f'    "{proc_id}" [label="{label}"];')

    lines.extend(
        [
            "",
            '    node [fillcolor="#c8e6c9", style="filled,rounded"];',
            "    // Outputs",
        ]
    )

    # Add outputs
    for out in pipeline["outputs"]:
        output_id = f"{out['name']}_out"
        label = out["file"]
        lines.append(f'    "{output_id}" [label="{label}"];')

    lines.append("\n    // Connections")

    # Add edges from processes
    for i, proc in enumerate(pipeline["processes"]):
        proc_id = proc["name"]
        params = proc["parameters"]
        for param_name, param_value in params.items():
            if isinstance(param_value, str):
                # Check if this references another node
                if any(inp["name"] == param_value for inp in pipeline["inputs"]) or any(
                    p["name"] == param_value for p in pipeline["processes"][:i]
                ):
                    lines.append(f'    "{param_value}" -> "{proc_id}";')

    # Add edges to outputs
    for out in pipeline["outputs"]:
        output_id = f"{out['name']}_out"
        lines.append(f'    "{out["name"]}" -> "{output_id}";')

    lines.append("}")
    return "\n".join(lines)


def generate_html_interactive(pipeline: dict) -> str:
    """Generate an interactive HTML visualization using Mermaid.js."""
    mermaid_code = generate_mermaid_diagram(pipeline)
    # Remove the backticks from mermaid code
    mermaid_code = mermaid_code.replace("```mermaid\n", "").replace("```", "")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pipeline Visualization: {pipeline["name"]}</title>
    <script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }}
        .info {{
            background: #e3f2fd;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .mermaid {{
            background: white;
            padding: 20px;
            border-radius: 5px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Pipeline: {pipeline["name"]}</h1>
        <div class="info">
            <strong>Inputs:</strong> {len(pipeline["inputs"])} |
            <strong>Processes:</strong> {len(pipeline["processes"])} |
            <strong>Outputs:</strong> {len(pipeline["outputs"])}
        </div>
        <div class="mermaid">
{mermaid_code}
        </div>
    </div>
    <script>
        mermaid.initialize({{ startOnLoad: true, theme: 'default' }});
    </script>
</body>
</html>"""
    return html


def main():
    parser = argparse.ArgumentParser(description="Visualize labUtils pipeline YAML files")
    parser.add_argument("yaml_file", type=str, help="Path to the pipeline YAML file")
    parser.add_argument(
        "--format",
        type=str,
        choices=["mermaid", "graphviz", "html"],
        default="mermaid",
        help="Output format (default: mermaid)",
    )
    parser.add_argument("--output", type=str, help="Output file path (default: stdout or pipeline.html)")
    parser.add_argument("--pipeline", type=str, help="Specific pipeline name to visualize (default: first pipeline)")

    args = parser.parse_args()

    yaml_path = Path(args.yaml_file)
    if not yaml_path.exists():
        print(f"Error: File not found: {yaml_path}", file=sys.stderr)
        sys.exit(1)

    # Parse the YAML file
    try:
        pipeline_data = parse_pipeline_yaml(yaml_path)
        pipelines = extract_pipeline_structure(pipeline_data)

        if not pipelines:
            print("Error: No pipelines found in YAML file", file=sys.stderr)
            sys.exit(1)

        # Select pipeline
        if args.pipeline:
            selected_pipeline = next((p for p in pipelines if p["name"] == args.pipeline), None)
            if not selected_pipeline:
                print(f"Error: Pipeline '{args.pipeline}' not found", file=sys.stderr)
                print(f"Available pipelines: {', '.join(p['name'] for p in pipelines)}", file=sys.stderr)
                sys.exit(1)
        else:
            selected_pipeline = pipelines[0]

        # Generate visualization
        if args.format == "mermaid":
            output = generate_mermaid_diagram(selected_pipeline)
            default_ext = ".md"
        elif args.format == "graphviz":
            output = generate_graphviz_dot(selected_pipeline)
            default_ext = ".dot"
        else:  # html
            output = generate_html_interactive(selected_pipeline)
            default_ext = ".html"

        # Write output
        if args.output:
            output_path = Path(args.output)
        else:
            if args.format == "html":
                output_path = yaml_path.with_suffix(default_ext).with_name(f"{yaml_path.stem}_pipeline{default_ext}")
            else:
                output_path = None

        if output_path:
            output_path.write_text(output, encoding="utf-8")
            print(f"Pipeline visualization saved to: {output_path}")
        else:
            print(output)

    except Exception as e:
        print(f"Error processing pipeline: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
