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

    # Get hover-to-expand settings from config
    truncate_settings = cfg.get("html_report_truncate_settings", {})
    # Default from config or a hardcoded fallback
    default_max_width = truncate_settings.get("default_max_width_px", 150)
    columns_for_hover_expand = truncate_settings.get("columns_for_hover_expand", [])
    column_specific_max_widths = truncate_settings.get("column_specific_max_widths_px", {})

    column_data_for_template = []
    if variants_data:
        # Get ordered columns from the first variant
        original_columns = list(variants_data[0].keys())

        for col_name in original_columns:
            # Skip the igv_links column as it will be handled specially
            if col_name == "igv_links":
                continue

            # Replace underscores for better readability
            display_name = col_name.replace("_", " ")
            is_standard_link = col_name in link_configs

            # Determine if this column should have hover-expand behavior
            apply_hover_behavior = col_name in columns_for_hover_expand
            hover_max_width = None
            if apply_hover_behavior:
                hover_max_width = column_specific_max_widths.get(col_name, default_max_width)

            column_data_for_template.append(
                {
                    "original_name": col_name,
                    "display_name": display_name,
                    "is_standard_link_column": is_standard_link,
                    "is_igv_link_column": False,
                    # Use original column name as link text for standard links
                    "link_display_text": col_name if is_standard_link else None,
                    # Hover-expand settings
                    "apply_hover_expand": apply_hover_behavior,
                    "max_width_px": hover_max_width,
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
                # IGV links don't need hover-expand as they're handled specially
                "apply_hover_expand": False,
                "max_width_px": None,
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
