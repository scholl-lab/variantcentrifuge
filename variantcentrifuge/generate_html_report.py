"""HTML report generation for VariantCentrifuge."""

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup

from .checkpoint_models import AtomicFileOperation

DANGEROUS_TAG_NAMES = {
    "script",
    "iframe",
    "object",
    "embed",
    "svg",
    "img",
    "body",
    "style",
    "link",
    "meta",
    "base",
    "form",
    "input",
    "button",
    "video",
    "audio",
    "source",
    "track",
    "frame",
    "frameset",
    "applet",
    "marquee",
    "details",
    "template",
    "math",
}


def _encode_symbolic_allele_brackets(url: Any) -> str:
    """Percent-encode symbolic VCF alleles like <DEL>, <DUP>, <INS:ME:ALU> in URL path."""
    if not url:
        return ""
    u = str(url)

    def _repl(match: re.Match[str]) -> str:
        tag = match.group(1)
        if tag.lower() in DANGEROUS_TAG_NAMES:
            return match.group(0)
        return f"%3C{tag}%3E"

    u = re.sub(r"<([A-Za-z0-9_:.-]+)>", _repl, u)
    if "<" not in u and not any(c in u for c in " \t\r\n\"'\\`"):
        u = u.replace(">", "%3E")
    return u


def is_safe_url(url: Any) -> bool:
    """Check if a URL uses a safe scheme (http, https, or relative path)."""
    if not url:
        return False
    u = _encode_symbolic_allele_brackets(str(url).strip())
    if any(ord(c) < 32 or ord(c) == 127 or c in " \t\r\n\"'<>\\" for c in u):
        return False
    if u.startswith("//"):
        return False
    u_lower = u.lower()
    if u_lower.startswith(("http://", "https://")):
        return True
    colon_idx = u_lower.find(":")
    slash_idx = u_lower.find("/")
    return colon_idx == -1 or (slash_idx != -1 and slash_idx < colon_idx)


def _coerce_positive_int(val: Any, default: int | None = None) -> int | None:
    """Coerce value to a positive finite integer (> 0), returning default on failure."""
    if isinstance(val, bool):
        return default
    try:
        n = int(val)
        return n if n > 0 else default
    except (ValueError, TypeError, OverflowError):
        return default


def _normalize_sentinel(val: Any) -> str:
    """Normalize missing or sentinel genomic values to empty string."""
    if val is None:
        return ""
    s = str(val).strip()
    if s in {"", ".", "N/A", "NA", "null", "None"}:
        return ""
    return s


def compute_variant_key(v: dict[str, Any]) -> str:
    """Generate a stable allele-aware key for a variant record.

    Incorporates CHROM, POS, REF, ALT, and transcript FEATUREID when available.
    Duplicate annotation rows share the key, while distinct alternate alleles
    or missing/repeated VAR_ID values produce distinct keys.
    """
    chrom = _normalize_sentinel(v.get("CHROM"))
    pos = _normalize_sentinel(v.get("POS"))
    ref = _normalize_sentinel(v.get("REF"))
    alt = _normalize_sentinel(v.get("ALT"))
    feature = _normalize_sentinel(
        v.get("FEATUREID") or v.get("Feature_ID") or v.get("TRANSCRIPT") or v.get("FEATURE")
    )

    if chrom and pos and ref and alt:
        if feature:
            return f"{chrom}:{pos}:{ref}:{alt}:{feature}"
        return f"{chrom}:{pos}:{ref}:{alt}"

    row_bytes = json.dumps(v, sort_keys=True, default=str).encode("utf-8")
    row_hash = hashlib.sha256(row_bytes).hexdigest()[:12]

    var_id = _normalize_sentinel(v.get("VAR_ID") or v.get("ID"))
    if var_id:
        return f"id:{var_id}:{row_hash}"

    return f"row:{row_hash}"


def compute_report_id(
    variants_data: list[dict[str, Any]],
    vcf_source: Any = "Unknown",
    reference_genome: Any = "Unknown",
    filter_expression: Any = "None",
    cfg_report_id: Any = None,
) -> str:
    """Generate a deterministic report dataset ID / fingerprint.

    If an explicit report_id is provided in configuration, it is used directly.
    Otherwise, a deterministic SHA256 hex digest is derived from provenance and
    the variants data.
    """
    if cfg_report_id:
        return str(cfg_report_id)

    hasher = hashlib.sha256()
    hasher.update(str(vcf_source).encode("utf-8"))
    hasher.update(str(reference_genome).encode("utf-8"))
    hasher.update(str(filter_expression).encode("utf-8"))
    variants_json = json.dumps(variants_data, sort_keys=True, default=str)
    hasher.update(variants_json.encode("utf-8"))
    return hasher.hexdigest()[:16]


def _load_assets() -> dict[str, str]:
    """Load vendored JS/CSS assets for inline embedding.

    Returns a dict with keys like 'js/datatables.min' and 'css/datatables.min'
    to avoid collisions between same-stem files in different subdirectories.
    """
    assets_dir = Path(__file__).parent / "assets"
    assets: dict[str, str] = {}
    for subdir in ("js", "css"):
        asset_path = assets_dir / subdir
        if asset_path.exists():
            for f in sorted(asset_path.iterdir()):
                if f.is_file():
                    # Namespace by subdirectory to prevent key collisions
                    assets[f"{subdir}/{f.stem}"] = Markup(f.read_text(encoding="utf-8"))
    return assets


def generate_html_report(
    variants_json: str, summary_json: str, output_dir: str, cfg: dict[str, Any]
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
    with open(variants_json, encoding="utf-8") as vf:
        variants_data = json.load(vf)
    with open(summary_json, encoding="utf-8") as sf:
        summary = json.load(sf)

    for v in variants_data:
        if "igv_links" in v and isinstance(v["igv_links"], list):
            for link in v["igv_links"]:
                if isinstance(link, dict) and "report_path" in link and link["report_path"]:
                    link["report_path"] = _encode_symbolic_allele_brackets(link["report_path"])

    if "reference" in cfg:
        from .links import resolve_links

        resolve_links(cfg)

    default_hidden_columns = cfg.get("html_report_default_hidden_columns", [])
    link_configs = cfg.get("links", {})  # For identifying link columns

    for v in variants_data:
        for col_name in link_configs:
            if col_name in v and isinstance(v[col_name], str) and v[col_name]:
                v[col_name] = _encode_symbolic_allele_brackets(v[col_name])

    # Get hover-to-expand settings from config
    truncate_settings = cfg.get("html_report_truncate_settings", {})
    # Default from config or a hardcoded fallback
    default_max_width = (
        _coerce_positive_int(truncate_settings.get("default_max_width_px"), default=150) or 150
    )
    columns_for_hover_expand = truncate_settings.get("columns_for_hover_expand", [])
    column_specific_max_widths = truncate_settings.get("column_specific_max_widths_px", {})

    column_data_for_template = []
    if variants_data:
        # Get ordered columns from the first variant
        original_columns = list(variants_data[0].keys())

        for col_name in original_columns:
            # Skip the igv_links column and internal metadata keys
            if col_name == "igv_links" or col_name.startswith("_"):
                continue

            # Replace underscores for better readability
            display_name = col_name.replace("_", " ")
            is_standard_link = col_name in link_configs

            # Determine if this column should have hover-expand behavior
            apply_hover_behavior = col_name in columns_for_hover_expand
            hover_max_width = None
            if apply_hover_behavior:
                raw_w = column_specific_max_widths.get(col_name)
                hover_max_width = (
                    _coerce_positive_int(raw_w, default=default_max_width)
                    if raw_w is not None
                    else default_max_width
                )

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

    env = Environment(
        loader=FileSystemLoader(str(templates_dir)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    env.filters["is_safe_url"] = is_safe_url
    env.filters["sanitize_url"] = _encode_symbolic_allele_brackets
    template = env.get_template("index.html")

    # Import version for template
    from variantcentrifuge.version import __version__

    # Load vendored assets for inline embedding
    assets = _load_assets()

    # Extract pipeline metadata from config
    filter_expression = cfg.get("filters") or cfg.get("filter") or "None"
    # Try to get VCF file path from various possible keys
    vcf_path = cfg.get("vcf_file") or cfg.get("input_vcf")
    if vcf_path:
        import os

        vcf_source = os.path.basename(vcf_path)
    else:
        vcf_source = "Unknown"
    reference_genome = cfg.get("reference", "Unknown")

    if cfg.get("pseudonymize"):
        from .stages.metadata_helpers import _sanitize_text_identifiers

        raw_map = cfg.get("pseudonymization_mapping")
        mapping = {str(k): str(v) for k, v in raw_map.items()} if isinstance(raw_map, dict) else {}
        redact = not mapping
        if vcf_source != "Unknown":
            vcf_source = _sanitize_text_identifiers(
                str(vcf_source), mapping, redact_unmapped=redact
            )
        if filter_expression != "None":
            filter_expression = _sanitize_text_identifiers(
                str(filter_expression), mapping, redact_unmapped=redact
            )
        if isinstance(summary.get("metadata"), dict):
            summary["metadata"] = {
                k: _sanitize_text_identifiers(str(v), mapping, redact_unmapped=redact)
                if isinstance(v, str)
                else v
                for k, v in summary["metadata"].items()
            }
            with (
                AtomicFileOperation(summary_json) as temp_sjf,
                open(temp_sjf, "w", encoding="utf-8") as sjf_out,
            ):
                json.dump(summary, sjf_out, indent=2)

    # Compute deterministic report fingerprint and assign stable variant keys
    report_id = compute_report_id(
        variants_data,
        vcf_source=vcf_source,
        reference_genome=reference_genome,
        filter_expression=filter_expression,
        cfg_report_id=cfg.get("report_id"),
    )
    for v in variants_data:
        v["_variant_key"] = compute_variant_key(v)

    # Render template with the new column structure
    html_content = template.render(
        variants=variants_data,
        summary=summary,
        column_data=column_data_for_template,
        default_hidden_columns=default_hidden_columns,
        generation_date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        version=__version__,
        assets=assets,
        filter_expression=filter_expression,
        vcf_source=vcf_source,
        reference_genome=reference_genome,
        report_id=report_id,
    )

    output_path = Path(output_dir) / "index.html"
    with (
        AtomicFileOperation(output_path) as temp_out,
        open(temp_out, "w", encoding="utf-8") as out_f,
    ):
        out_f.write(html_content)
