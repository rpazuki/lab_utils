# Pipeline Visualization

Visualize labUtils YAML pipeline configurations as interactive diagrams.

## Installation

```bash
pip install pyyaml
```

## Usage

### Command Line

Generate an interactive HTML visualization:
```bash
python src/labUtils/scripts/visualize_pipeline.py src/labUtils/yamls/growth_rates_pipeline.yaml --format html --output pipeline.html
```

Generate Mermaid markdown:
```bash
python src/labUtils/scripts/visualize_pipeline.py src/labUtils/yamls/growth_rates_pipeline.yaml --format mermaid --output pipeline.md
```

Generate Graphviz DOT file:
```bash
python src/labUtils/scripts/visualize_pipeline.py src/labUtils/yamls/growth_rates_pipeline.yaml --format graphviz --output pipeline.dot
```

### Python API

```python
from labUtils.pipeline_viz import visualize_pipeline, list_pipelines

# Generate HTML visualization
visualize_pipeline('src/labUtils/yamls/growth_rates_pipeline.yaml',
                  format='html',
                  output='my_pipeline.html')

# Get Mermaid diagram as string
mermaid_diagram = visualize_pipeline('src/labUtils/yamls/growth_rates_pipeline.yaml',
                                     format='mermaid')
print(mermaid_diagram)

# List all pipelines in a file
pipelines = list_pipelines('src/labUtils/yamls/growth_rates_pipeline.yaml')
print(pipelines)

# Visualize a specific pipeline
visualize_pipeline('src/labUtils/yamls/growth_rates_pipeline.yaml',
                  pipeline_name='growth_rate_replicates_fit_pipeline',
                  format='html',
                  output='replicates_pipeline.html')
```

## Output Formats

### HTML (Recommended)
- Interactive diagram using Mermaid.js
- Opens directly in web browser
- No additional tools needed
- Best for sharing and presentations

### Mermaid
- Markdown with Mermaid syntax
- Can be rendered in VS Code with Mermaid extensions
- GitHub/GitLab will render it automatically
- Portable and version-control friendly

### Graphviz
- DOT format for Graphviz
- Requires Graphviz installation to render
- Most customizable
- Professional quality output

```bash
# Install Graphviz to render DOT files
# Windows: choco install graphviz
# Mac: brew install graphviz
# Linux: sudo apt install graphviz

# Render DOT file to PNG
dot -Tpng pipeline.dot -o pipeline.png
```

## Features

- **Automatic dependency detection**: Traces data flow through inputs → processes → outputs
- **Color-coded nodes**:
  - Blue: Input data sources
  - Yellow: Processing steps
  - Green: Output files
- **Method information**: Shows which package/method each step uses
- **Interactive HTML**: Zoom, pan, and explore large pipelines

## VS Code Integration

To view Mermaid diagrams in VS Code:
1. Install the "Markdown Preview Mermaid Support" extension
2. Generate a `.md` file with mermaid format
3. Open the file and use Ctrl+Shift+V to preview

## Examples

See the generated visualizations:
- `growth_rates_pipeline.html` - Complete growth rate fitting pipeline
- `amn_pipeline.html` - Metabolic mapping pipeline


## Ready to run
- `python visualize_pipeline.py ../yamls/growth_rates_pipeline.yaml --format html --output growth_rates_pipeline.html`
- `python visualize_pipeline.py ../yamls/collateing_pipeline.yaml --format html --output collateing_pipeline.html`
- `python visualize_pipeline.py ../yamls/amn_pipeline.yaml --format html --output amn_pipeline.html`
