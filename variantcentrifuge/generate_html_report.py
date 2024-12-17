# File: variantcentrifuge/generate_html_report.py
# Location: variantcentrifuge/generate_html_report.py

import json
from pathlib import Path
from typing import Optional
from jinja2 import Environment, FileSystemLoader

def generate_html_report(
    variants_json: str,
    summary_json: str,
    output_dir: str
) -> None:
    """
    Generate an interactive HTML report of variants and summary statistics.

    Parameters
    ----------
    variants_json : str
        Path to the variants JSON file containing detailed variant records.
    summary_json : str
        Path to the summary JSON file containing summary metrics.
    output_dir : str
        Directory to write the generated HTML report and assets.

    Returns
    -------
    None
    """
    # Load data
    with open(variants_json, "r") as vf:
        variants = json.load(vf)
    with open(summary_json, "r") as sf:
        summary = json.load(sf)

    # Determine columns based on the first variant's key order
    if variants:
        # Python 3.7+ guarantees that dict keys preserve insertion order.
        # The order should match the order in the original JSON.
        columns = list(variants[0].keys())
    else:
        columns = []

    templates_dir = Path(__file__).parent / "templates"
    if not templates_dir.exists():
        raise FileNotFoundError(f"Templates directory not found at: {templates_dir}")

    env = Environment(loader=FileSystemLoader(str(templates_dir)))
    template = env.get_template("index.html")

    # Render template with columns and variants
    html_content = template.render(
        variants=variants,
        summary=summary,
        columns=columns
    )

    output_path = Path(output_dir) / "index.html"
    with open(output_path, "w") as out_f:
        out_f.write(html_content)
