# File: variantcentrifuge/generate_html_report.py
# Location: variantcentrifuge/generate_html_report.py

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from jinja2 import Environment, FileSystemLoader


def generate_html_report(
    variants_json: str, summary_json: str, output_dir: str, cfg: Dict[str, Any]
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
    cfg : Dict[str, Any]
        The main configuration dictionary containing report settings.

    Returns
    -------
    None
    """
    # Load data
    with open(variants_json, "r", encoding="utf-8") as vf:
        variants_data = json.load(vf)
    with open(summary_json, "r", encoding="utf-8") as sf:
        summary = json.load(sf)

    default_hidden_columns = cfg.get("html_report_default_hidden_columns", [])
    link_configs = cfg.get("links", {})  # For identifying link columns

    column_data_for_template = []
    if variants_data:
        # Get ordered columns from the first variant
        original_columns = list(variants_data[0].keys())

        for col_name in original_columns:
            # Skip the igv_links column as it will be handled specially
            if col_name == "igv_links":
                continue

            display_name = col_name.replace("_", " ")  # Replace underscores for better readability
            is_standard_link = col_name in link_configs

            column_data_for_template.append(
                {
                    "original_name": col_name,
                    "display_name": display_name,
                    "is_standard_link_column": is_standard_link,
                    "is_igv_link_column": False,
                    # Use original column name as link text for standard links
                    "link_display_text": col_name if is_standard_link else None,
                }
            )

        # Add special IGV Links column
        column_data_for_template.append(
            {
                "original_name": "igv_links",
                "display_name": "IGV Links",
                "is_standard_link_column": False,
                "is_igv_link_column": True,
                "link_display_text": None,
            }
        )
    else:
        # If variants_data is empty, we'll pass an empty column_data list
        pass

    templates_dir = Path(__file__).parent / "templates"
    if not templates_dir.exists():
        raise FileNotFoundError(f"Templates directory not found at: {templates_dir}")

    env = Environment(loader=FileSystemLoader(str(templates_dir)))
    template = env.get_template("index.html")

    # Render template with the new column structure
    html_content = template.render(
        variants=variants_data,
        summary=summary,
        column_data=column_data_for_template,
        default_hidden_columns=default_hidden_columns,
    )

    output_path = Path(output_dir) / "index.html"
    with open(output_path, "w", encoding="utf-8") as out_f:
        out_f.write(html_content)
