"""Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

import os
import sys

from variantcentrifuge.version import __version__

# Add the project root to the path so autodoc can find the variantcentrifuge module
sys.path.insert(0, os.path.abspath("../.."))

# Add custom extensions directory
sys.path.insert(0, os.path.abspath("_ext"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "VariantCentrifuge"
copyright = "2024-2026, Bernt Popp"
author = "Bernt Popp"

# The full version, including alpha/beta/rc tags
# This will automatically get the version from your package
release = __version__
version = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # Core library for autodoc
    "sphinx.ext.autosummary",  # Generate autodoc summaries
    "sphinx.ext.napoleon",  # Support for Google-style docstrings
    "sphinx.ext.viewcode",  # Add links to highlighted source code
    "sphinx.ext.intersphinx",  # Link to other project documentations
    "myst_parser",  # Enable writing docs in Markdown
    "sphinx_autodoc_typehints",  # Automatically document typehints
    "sphinxcontrib.mermaid",  # Mermaid diagram support
    "sphinx_seo",  # Custom SEO extension
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["guides/index.md"]

# The suffix(es) of source filenames.
source_suffix = {
    ".rst": None,
    ".md": "markdown",
}

# The master toctree document.
master_doc = "index"

language = "en"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# SEO Configuration
html_baseurl = "https://scholl-lab.github.io/variantcentrifuge/"
html_title = "VariantCentrifuge - Clinical Variant Analysis Tool"

# Theme options for Furo
html_theme_options = {
    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#2980b9",
        "color-brand-content": "#2980b9",
    },
    "dark_css_variables": {
        "color-brand-primary": "#3498db",
        "color-brand-content": "#3498db",
    },
}

# -- Extension configuration -------------------------------------------------

# Napoleon settings for docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
}

# Autosummary settings
autosummary_generate = True

# MyST-Parser configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
}

# Type hints configuration
typehints_fully_qualified = False
always_document_param_types = True
typehints_document_rtype = True

# Build settings - continue on warnings
suppress_warnings = [
    "autodoc.import_object",
    "docutils.parsers.rst",
    "docutils.frontend",
    "docutils.nodes",
    "docutils.utils",
    "myst",
    "toc.not_included",
]
nitpicky = False

# Configure to not treat warnings as errors
keep_warnings = True

# Linkcheck settings
linkcheck_ignore = [
    # GitHub login redirects are expected
    r"https://github\.com/login.*",
    # Ignore local file links that may not exist during CI
    r"file://.*",
]
linkcheck_retries = 2
linkcheck_timeout = 10
