# File: variantcentrifuge/links.py
# Location: variantcentrifuge/variantcentrifuge/links.py

"""
Module for adding configurable links to output tables.

This module provides functionality to add new columns containing URL links
to each variant line in the final output. The URLs and their substitution
logic are defined in the configuration file (config.json) under "links" or
"links_by_build".

For example, if you have a link defined as::

    "links": {
        "SpliceAI": "https://spliceailookup.broadinstitute.org/#variant={CHROM}-{POS}-{REF}-{ALT}&hg=38&bc=basic&distance=500&mask=0&ra=0"
    }

This will create a new column named "SpliceAI" with a URL where {CHROM}, {POS}, {REF}, {ALT}
are replaced by the corresponding values from each row.
"""

from typing import Any

DEFAULT_LINKS_BY_BUILD: dict[str, dict[str, str]] = {
    "GRCh37": {
        "ClinVar": "https://www.ncbi.nlm.nih.gov/clinvar/?term={CHROM}-{POS}-{REF}-{ALT}",
        "gnomAD_2": "https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}?dataset=gnomad_r2_1",
        "Varsome": "https://varsome.com/variant/hg19/{CHROM}-{POS}-{REF}-{ALT}",
        "Franklin": "https://franklin.genoox.com/clinical-db/variant/snp/{CHROM}-{POS}-{REF}-{ALT}-hg19",
        "SpliceAI": "https://spliceailookup.broadinstitute.org/#variant={CHROM}-{POS}-{REF}-{ALT}&hg=19&bc=basic&distance=500&mask=0&ra=0",
        "autopvs1": "https://autopvs1.bgi.com/variant/hg19/{CHROM}-{POS}-{REF}-{ALT}",
    },
    "GRCh38": {
        "ClinVar": "https://www.ncbi.nlm.nih.gov/clinvar/?term={CHROM}-{POS}-{REF}-{ALT}",
        "gnomAD_4": "https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}?dataset=gnomad_r4",
        "Varsome": "https://varsome.com/variant/hg38/{CHROM}-{POS}-{REF}-{ALT}",
        "Franklin": "https://franklin.genoox.com/clinical-db/variant/snp/{CHROM}-{POS}-{REF}-{ALT}-hg38",
        "SpliceAI": "https://spliceailookup.broadinstitute.org/#variant={CHROM}-{POS}-{REF}-{ALT}&hg=38&bc=basic&distance=500&mask=0&ra=0",
        "autopvs1": "https://autopvs1.bgi.com/variant/hg38/{CHROM}-{POS}-{REF}-{ALT}",
    },
}


def normalize_reference(reference: str | None) -> str | None:
    """
    Normalize a reference genome identifier to canonical assembly name.

    Returns 'GRCh37' for hg19/GRCh37 references and 'GRCh38' for hg38/GRCh38
    references, ignoring version suffixes (e.g., GRCh37.75, GRCh38.99).
    Returns None for unsupported or non-human assemblies.
    """
    if not reference:
        return None
    ref_clean = reference.strip().lower()
    if ref_clean.startswith("grch37") or ref_clean.startswith("hg19"):
        return "GRCh37"
    if ref_clean.startswith("grch38") or ref_clean.startswith("hg38"):
        return "GRCh38"
    return None


def resolve_links(config: dict[str, Any]) -> dict[str, str]:
    """
    Resolve external link configurations based on genome build and explicit configuration.

    If 'links' is explicitly configured with custom templates (i.e. not matching the bundled
    defaults), those custom templates are preserved.
    Otherwise, links are resolved from 'links_by_build' based on 'reference'.
    If the reference build is unsupported or unknown, bundled links return an empty dict to prevent
    silently generating hg19 links from GRCh38 or non-human coordinates.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing 'reference', optional 'links', and optional 'links_by_build'.

    Returns
    -------
    dict of str to str
        The resolved link templates mapping link names to URL templates.
    """
    links_by_build = config.get("links_by_build") or DEFAULT_LINKS_BY_BUILD
    canonical_build = normalize_reference(config.get("reference"))
    build_defaults = dict(links_by_build.get(canonical_build, {})) if canonical_build else {}

    raw_links = config.get("links")

    if raw_links is None:
        resolved = build_defaults
    elif isinstance(raw_links, dict):
        all_defaults: dict[str, set[str]] = {}
        for b_links in links_by_build.values():
            for k, v in b_links.items():
                all_defaults.setdefault(k, set()).add(v)
        for b_links in DEFAULT_LINKS_BY_BUILD.values():
            for k, v in b_links.items():
                all_defaults.setdefault(k, set()).add(v)

        resolved = {}
        gnomad_defaults = (
            all_defaults.get("gnomAD", set())
            | all_defaults.get("gnomAD_2", set())
            | all_defaults.get("gnomAD_4", set())
        )
        for name, url_template in raw_links.items():
            if name in ("gnomAD", "gnomAD_2", "gnomAD_4") and url_template in gnomad_defaults:
                if canonical_build:
                    gnomad_key = next((k for k in build_defaults if k.startswith("gnomAD")), None)
                    if gnomad_key and gnomad_key in build_defaults:
                        resolved[gnomad_key] = build_defaults[gnomad_key]
            elif name in all_defaults and url_template in all_defaults[name]:
                if canonical_build and name in build_defaults:
                    resolved[name] = build_defaults[name]
            else:
                resolved[name] = url_template
    else:
        resolved = {}

    config["links"] = resolved
    return resolved


def add_links_to_table(lines: list[str], link_configs: dict[str, str]) -> list[str]:
    """
    Add link columns to the table lines based on provided link configurations.

    Parameters
    ----------
    lines : list of str
        The lines of the final table (including header).
    link_configs : dict of str to str
        A dictionary where keys are link column names and values are URL templates.
        URL templates contain placeholders like {CHROM}, {POS}, {REF}, {ALT} that
        are replaced with values from the corresponding columns.

    Returns
    -------
    list of str
        The modified lines with new link columns appended.

    Notes
    -----
    - Assumes the first line is a header line.
    - Each link column is appended to the end of the table.
    - Placeholders in the link template must match header column names exactly (case-sensitive).
    - If a placeholder column doesn't exist in the header, the link will be empty.
    """
    if not link_configs:
        return lines

    if not lines:
        return lines

    header = lines[0].rstrip("\n").split("\t")
    # Map column names to their indices for substitution
    col_index = {name: i for i, name in enumerate(header)}

    # Extend the header with new link columns
    new_header = header[:]
    link_names = sorted(link_configs.keys())
    for link_name in link_names:
        new_header.append(link_name)

    modified_lines = ["\t".join(new_header)]

    for line in lines[1:]:
        line = line.rstrip("\n")
        if not line.strip():
            # Empty or blank line
            modified_lines.append(line)
            continue
        fields = line.split("\t")

        # For each link, perform substitution
        for link_name in link_names:
            template = link_configs[link_name]
            # Replace placeholders with corresponding column values
            # If a column is missing, produce empty string for that link
            link_url = template
            for placeholder in col_index:
                placeholder_token = "{" + placeholder + "}"
                if placeholder_token in template:
                    if col_index[placeholder] < len(fields):
                        link_url = link_url.replace(
                            placeholder_token, fields[col_index[placeholder]]
                        )
                    else:
                        # Column index out of range, set link_url empty
                        link_url = ""
                        break
            fields.append(link_url)

        modified_lines.append("\t".join(fields))

    return modified_lines
