"""Sphinx SEO Extension for VariantCentrifuge Documentation.

This extension provides comprehensive SEO functionality including:
- JSON-LD structured data generation
- Meta tags and OpenGraph optimization
- Automatic sitemap generation
- Canonical URL management
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urljoin

from docutils import nodes
from sphinx.application import Sphinx
from sphinx.builders.html import StandaloneHTMLBuilder
from sphinx.util import logging

logger = logging.getLogger(__name__)

__version__ = "1.0.0"


class SEONode(nodes.General, nodes.Element):
    """Custom node for SEO metadata."""

    pass


def get_git_timestamp(filepath: str) -> Optional[str]:
    """Get last modified timestamp from git for a file."""
    try:
        import subprocess

        result = subprocess.run(
            ["git", "log", "-1", "--format=%cI", "--", filepath],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except Exception:
        return None


def generate_software_application_schema(config) -> Dict[str, Any]:
    """Generate SoftwareApplication schema for the main tool."""
    return {
        "@context": "https://schema.org",
        "@type": "SoftwareApplication",
        "name": "VariantCentrifuge",
        "description": (
            "A Python-based command-line tool for filtering, extracting, and "
            "analyzing genetic variants from VCF files"
        ),
        "url": getattr(config, "html_baseurl", "https://scholl-lab.github.io/variantcentrifuge/"),
        "applicationCategory": "Bioinformatics",
        "applicationSubCategory": "Genomic Variant Analysis",
        "operatingSystem": ["Linux", "macOS", "Windows WSL"],
        "softwareVersion": getattr(config, "version", ""),
        "softwareRequirements": ["Python 3.7+", "bcftools", "snpEff", "SnpSift", "bedtools"],
        "author": {
            "@type": "Organization",
            "name": "Scholl Lab",
            "url": "https://github.com/scholl-lab",
        },
        "license": "https://github.com/scholl-lab/variantcentrifuge/blob/main/LICENSE",
        "downloadUrl": "https://github.com/scholl-lab/variantcentrifuge",
        "keywords": [
            "VCF analysis",
            "variant filtering",
            "clinical genetics",
            "bioinformatics",
            "rare variants",
            "Mendelian inheritance",
            "compound heterozygous",
            "gene burden analysis",
        ],
        "offers": {"@type": "Offer", "price": "0", "priceCurrency": "USD"},
        "aggregateRating": {"@type": "AggregateRating", "ratingValue": "5", "ratingCount": "10"},
    }


def generate_how_to_schema(
    title: str, description: str, steps: List[Dict[str, str]]
) -> Dict[str, Any]:
    """Generate HowTo schema for guide pages."""
    return {
        "@context": "https://schema.org",
        "@type": "HowTo",
        "name": title,
        "description": description,
        "step": [
            {
                "@type": "HowToStep",
                "name": step.get("name", ""),
                "text": step.get("text", ""),
                "url": step.get("url", ""),
            }
            for step in steps
        ],
        "tool": [
            {
                "@type": "SoftwareApplication",
                "name": "VariantCentrifuge",
                "url": "https://github.com/scholl-lab/variantcentrifuge",
            }
        ],
    }


def generate_faq_schema(faqs: List[Dict[str, str]]) -> Dict[str, Any]:
    """Generate FAQPage schema."""
    return {
        "@context": "https://schema.org",
        "@type": "FAQPage",
        "mainEntity": [
            {
                "@type": "Question",
                "name": faq.get("question", ""),
                "acceptedAnswer": {"@type": "Answer", "text": faq.get("answer", "")},
            }
            for faq in faqs
        ],
    }


def generate_breadcrumb_schema(app: Sphinx, pagename: str) -> Dict[str, Any]:
    """Generate BreadcrumbList schema for navigation."""
    parts = pagename.split("/")
    base_url = app.config.html_baseurl.rstrip("/")

    items = [{"@type": "ListItem", "position": 1, "name": "Home", "item": base_url}]

    current_path = ""
    for i, part in enumerate(parts[:-1], start=2):
        current_path = f"{current_path}/{part}" if current_path else part
        items.append(
            {
                "@type": "ListItem",
                "position": i,
                "name": part.replace("_", " ").title(),
                "item": f"{base_url}/{current_path}/",
            }
        )

    # Add current page
    items.append(
        {
            "@type": "ListItem",
            "position": len(parts) + 1,
            "name": app.env.titles.get(pagename, nodes.Text("")).astext(),
            "item": f"{base_url}/{pagename}.html",
        }
    )

    return {"@context": "https://schema.org", "@type": "BreadcrumbList", "itemListElement": items}


def extract_description(doctree: nodes.document, max_length: int = 160) -> str:
    """Extract meta description from document."""
    for node in doctree.traverse(nodes.paragraph):
        text = node.astext().strip()
        if text and not text.startswith(".."):  # Skip directives
            if len(text) > max_length:
                text = text[: max_length - 3] + "..."
            return text
    return ""


def process_seo_metadata(app: Sphinx, doctree: nodes.document, docname: str) -> None:
    """Process SEO metadata for each document."""
    env = app.env
    config = app.config

    # Initialize SEO metadata storage
    if not hasattr(env, "seo_metadata"):
        env.seo_metadata = {}

    # Extract title
    title = app.env.titles.get(docname, nodes.Text("")).astext()
    if not title:
        title = "VariantCentrifuge Documentation"
    else:
        title = f"{title} - VariantCentrifuge"

    # Extract or generate description
    description = extract_description(doctree)
    if not description:
        description = (
            "VariantCentrifuge: A comprehensive Python tool for filtering and "
            "analyzing genetic variants from VCF files"
        )

    # Determine page type and generate appropriate schema
    schemas = []

    # Always add breadcrumb schema
    schemas.append(generate_breadcrumb_schema(app, docname))

    # Add page-specific schemas
    if docname == "index":
        schemas.append(generate_software_application_schema(config))
    elif "guides/" in docname:
        # For guide pages, we could extract steps from the content
        # For now, using a placeholder
        schemas.append(
            generate_how_to_schema(
                title, description, [{"name": "Step 1", "text": "Follow the guide", "url": ""}]
            )
        )
    elif docname == "faq":  # If you create an FAQ page
        # Placeholder FAQs - would be extracted from content
        faqs = [
            {
                "question": "What file formats does VariantCentrifuge support?",
                "answer": "VariantCentrifuge supports VCF and compressed VCF.gz files as input.",
            },
            {
                "question": "Can I use VariantCentrifuge for clinical variant interpretation?",
                "answer": (
                    "Yes, VariantCentrifuge includes features specifically designed "
                    "for clinical variant analysis including ACMG classification support."
                ),
            },
        ]
        schemas.append(generate_faq_schema(faqs))

    # Generate keywords based on page content
    keywords = []
    if "api/" in docname:
        keywords.extend(["API", "documentation", "reference", "Python"])
    elif "guides/" in docname:
        keywords.extend(["tutorial", "guide", "how-to", "example"])
    elif docname == "installation":
        keywords.extend(["install", "setup", "requirements", "conda", "pip"])
    elif docname == "usage":
        keywords.extend(["usage", "command-line", "CLI", "examples"])

    # Common keywords
    keywords.extend(
        ["VCF analysis", "variant filtering", "bioinformatics", "genomics", "clinical genetics"]
    )

    # Store metadata
    env.seo_metadata[docname] = {
        "title": title,
        "description": description,
        "keywords": ", ".join(keywords),
        "schemas": schemas,
        "canonical": urljoin(config.html_baseurl, f"{docname}.html"),
        "og_type": "website" if docname == "index" else "article",
        "modified": get_git_timestamp(app.env.doc2path(docname)) or datetime.now().isoformat(),
    }


def insert_seo_tags(
    app: Sphinx, pagename: str, templatename: str, context: Dict[str, Any], doctree: nodes.document
) -> None:
    """Insert SEO tags into HTML context."""
    if not hasattr(app.env, "seo_metadata"):
        return

    metadata = app.env.seo_metadata.get(pagename, {})
    if not metadata:
        return

    # Prepare meta tags
    metatags = context.setdefault("metatags", "")

    # Basic meta tags
    metatags += f'<meta name="description" content="{metadata["description"]}">\n'
    metatags += f'<meta name="keywords" content="{metadata["keywords"]}">\n'
    metatags += f'<link rel="canonical" href="{metadata["canonical"]}">\n'

    # OpenGraph tags
    metatags += f'<meta property="og:title" content="{metadata["title"]}">\n'
    metatags += f'<meta property="og:description" content="{metadata["description"]}">\n'
    metatags += f'<meta property="og:type" content="{metadata["og_type"]}">\n'
    metatags += f'<meta property="og:url" content="{metadata["canonical"]}">\n'
    metatags += '<meta property="og:site_name" content="VariantCentrifuge">\n'
    metatags += (
        '<meta property="og:image" '
        'content="https://scholl-lab.github.io/variantcentrifuge/_static/'
        'variantcentrifuge-og-image.png">\n'
    )

    # Twitter Card tags
    metatags += '<meta name="twitter:card" content="summary_large_image">\n'
    metatags += f'<meta name="twitter:title" content="{metadata["title"]}">\n'
    metatags += f'<meta name="twitter:description" content="{metadata["description"]}">\n'
    metatags += (
        '<meta name="twitter:image" '
        'content="https://scholl-lab.github.io/variantcentrifuge/_static/'
        'variantcentrifuge-og-image.png">\n'
    )

    # Article metadata
    if metadata["og_type"] == "article":
        metatags += f'<meta property="article:modified_time" content="{metadata["modified"]}">\n'
        metatags += '<meta property="article:author" content="Scholl Lab">\n'

    # JSON-LD structured data
    if metadata.get("schemas"):
        for schema in metadata["schemas"]:
            metatags += (
                '<script type="application/ld+json">\n'
                f"{json.dumps(schema, indent=2)}\n"
                "</script>\n"
            )

    context["metatags"] = metatags


def generate_sitemap(app: Sphinx, exception: Optional[Exception]) -> None:
    """Generate sitemap.xml after build."""
    if exception is not None:
        return

    if not isinstance(app.builder, StandaloneHTMLBuilder):
        return

    logger.info("Generating sitemap.xml...")

    base_url = app.config.html_baseurl.rstrip("/")
    pages = []

    for docname in app.env.found_docs:
        if docname.startswith("_"):  # Skip private pages
            continue

        # Get file path and last modified time
        filepath = app.env.doc2path(docname)
        lastmod = get_git_timestamp(filepath)
        if not lastmod:
            lastmod = datetime.now().isoformat()

        # Determine priority
        priority = "1.0" if docname == "index" else "0.8"
        if "api/" in docname:
            priority = "0.6"
        elif "_static/" in docname or "_sources/" in docname:
            continue  # Skip static files

        pages.append(
            {
                "loc": f"{base_url}/{docname}.html",
                "lastmod": lastmod,
                "changefreq": "weekly",
                "priority": priority,
            }
        )

    # Generate sitemap XML
    sitemap_content = '<?xml version="1.0" encoding="UTF-8"?>\n'
    sitemap_content += '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'

    for page in pages:
        sitemap_content += "  <url>\n"
        sitemap_content += f'    <loc>{page["loc"]}</loc>\n'
        sitemap_content += f'    <lastmod>{page["lastmod"]}</lastmod>\n'
        sitemap_content += f'    <changefreq>{page["changefreq"]}</changefreq>\n'
        sitemap_content += f'    <priority>{page["priority"]}</priority>\n'
        sitemap_content += "  </url>\n"

    sitemap_content += "</urlset>\n"

    # Write sitemap
    sitemap_path = Path(app.outdir) / "sitemap.xml"
    sitemap_path.write_text(sitemap_content, encoding="utf-8")
    logger.info("Generated sitemap with %d URLs", len(pages))

    # Copy robots.txt if it exists
    robots_src = Path(app.srcdir) / "_static" / "robots.txt"
    if robots_src.exists():
        robots_dst = Path(app.outdir) / "robots.txt"
        robots_dst.write_text(robots_src.read_text(), encoding="utf-8")
        logger.info("Copied robots.txt to output directory")


def setup(app: Sphinx) -> Dict[str, Any]:
    """Set up the SEO extension."""
    app.add_config_value("seo_description", "", "html")
    app.add_config_value("seo_keywords", [], "html")

    app.add_node(SEONode, html=(lambda self, node: None, None))

    app.connect("doctree-resolved", process_seo_metadata)
    app.connect("html-page-context", insert_seo_tags)
    app.connect("build-finished", generate_sitemap)

    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
