# File: variantcentrifuge/replacer.py
# Location: variantcentrifuge/variantcentrifuge/replacer.py

"""
Genotype replacement module.

This module dynamically searches for a "GT" column in the header line of the
input TSV and applies configurable genotype replacement logic.

2. For each variant row:
   - Parse the genotype subfields (split by `cfg["extract_fields_separator"]`, usually ":")
   - Skip non-variant genotypes ("0/0" or "./.").
   - If a `genotype_replacement_map` is provided (e.g. {r"[2-9]": "1"}), apply regex-based replacements
     to transform the allele strings.

3. If `cfg["append_extra_sample_fields"]` is True, the user-supplied fields in `cfg["extra_sample_fields"]`
   are appended to each genotype. For instance, if extra fields are "GEN[*].DP" and "GEN[*].AD",
   and the delimiter is ":", then an output entry might look like:
       sampleName(0/1:100:52,48)

4. The sample names come from `cfg["sample_list"]` (comma-separated). We pair each subfield in GT
   with the corresponding sample. We do not list or output samples if the genotype is "0/0" or "./.".

5. The final joined string for each variant is placed back into the "GT" column, separated from other
   variants by `cfg["separator"]` (commonly ";").

Header Normalization:
---------------------
SnpSift or other tools may remove "GEN[*]." or "ANN[*]." prefixes in headers. To handle this,
we normalize both the final TSV header columns and the user-supplied extra fields so they match.

Example:
--------
If the user passed --append-extra-sample-fields GEN[*].DP GEN[*].AD and the TSV header has columns
"DP", "AD" (prefix removed by SnpSift), the code will still find them by normalizing them to "DP", "AD"
internally.
"""

import logging
import re
from typing import Any, Dict, Iterator

logger = logging.getLogger("variantcentrifuge")


def _normalize_snpeff_column_name(col_name: str) -> str:
    """
    Normalize known SnpEff or SnpSift prefixes from a single column name.

    For example, "GEN[*].DP" or "ANN[*].Effect" might become "DP" or "Effect".
    Modify or extend this if your pipeline uses additional/different prefixes.
    """
    return (
        col_name.replace("ANN[*].", "")
        .replace("ANN[0].", "")
        .replace("GEN[*].", "")
        .replace("GEN[0].", "")
        .replace("NMD[*].", "NMD_")
        .replace("NMD[0].", "NMD_")
        .replace("LOF[*].", "LOF_")
        .replace("LOF[0].", "LOF_")
    )


def replace_genotypes(lines: Iterator[str], cfg: Dict[str, Any]) -> Iterator[str]:
    """
    Replace and optionally append extra sample fields to genotypes in a TSV file.

    Requirements:
    -------------
    - A "GT" column must exist in the header.
    - `cfg["sample_list"]` is a comma-separated list of sample names matching the subfields in "GT".
    - If `cfg["append_extra_sample_fields"]` is True and `cfg["extra_sample_fields"]` is a list of
      columns, we will look them up by matching normalized column names in the header.

    Steps:
    ------
    1. Identify the index of the "GT" column in the TSV header.
    2. Build a map from normalized column names to the actual index in the header.
    3. For each variant line:
       - Split the "GT" column subfields (one subfield per sample).
       - Optionally apply genotype replacements from `cfg["genotype_replacement_map"]`.
       - Filter out genotypes that are "0/0" or "./.".
       - If extra fields exist, append them in parentheses next to each genotype, e.g. "DP:100:52,48".
       - Rejoin them with `cfg["separator"]` (often ";") for final placement in the "GT" column.

    Returns
    -------
    Iterator[str]
        Iterator of updated lines (strings).
    """
    # Basic config
    separator = cfg.get("separator", ";")
    genotype_replacement_map = cfg.get("genotype_replacement_map", {r"[2-9]": "1"})

    # Extra fields config
    append_extra_fields = cfg.get("append_extra_sample_fields", False)
    extra_sample_fields = cfg.get("extra_sample_fields", [])
    extra_field_delimiter = cfg.get("extra_sample_field_delimiter", ":")

    # The real multi-sample subfield separator used by SnpSift for columns in the final TSV
    snpsift_sep = cfg.get("extract_fields_separator", ":")

    logger.debug("replace_genotypes called with config: %s", cfg)
    logger.debug("Extra fields => %s", extra_sample_fields)

    # Sample list
    sample_list_str = cfg.get("sample_list", "")
    samples = [s.strip() for s in sample_list_str.split(",") if s.strip()]
    logger.debug("Parsed samples => %s", samples)
    if not samples:
        logger.warning("No samples found in cfg[sample_list]. Returning lines unchanged.")
        yield from lines
        return

    # Compile genotype replacement regex patterns
    compiled_patterns = []
    for pat, repl in genotype_replacement_map.items():
        compiled = re.compile(pat)
        compiled_patterns.append((compiled, repl))

    # We will store the index of the GT column and the indexes of each extra field
    gt_idx = None
    extra_field_indices = {}

    first_line = True
    for line_idx, line in enumerate(lines, start=1):
        line = line.rstrip("\n")

        # Header line => identify GT and any extra fields
        if first_line:
            header_cols = line.split("\t")
            logger.debug("Header columns at line_idx=%d => %s", line_idx, header_cols)

            # Build a map from normalized name => real index
            normalized_to_index = {}
            for i, c in enumerate(header_cols):
                norm_c = _normalize_snpeff_column_name(c)
                normalized_to_index[norm_c] = i

            # Find the GT column
            if "GT" in header_cols:
                gt_idx = header_cols.index("GT")
                logger.debug("Found 'GT' column at index %d (exact match).", gt_idx)
            else:
                # Fallback: maybe the column was normalized by something else
                # e.g. "GEN[*].GT" => "GT"
                # But if the user had a different naming, we can try the normalized map
                # Usually it's just "GT", though
                if "GT" in normalized_to_index:
                    gt_idx = normalized_to_index["GT"]
                    logger.debug("Found 'GT' column at index %d (via normalized map).", gt_idx)

            if gt_idx is None:
                logger.warning("No GT column found in header. Returning lines unmodified.")
                yield line
                # yield the rest of the lines
                yield from lines
                return

            # Identify requested extra fields
            if append_extra_fields and extra_sample_fields:
                # For each raw field the user provided, find its normalized match in the header
                for raw_field in extra_sample_fields:
                    nf = _normalize_snpeff_column_name(raw_field)
                    if nf in normalized_to_index:
                        ef_col_idx = normalized_to_index[nf]
                        extra_field_indices[raw_field] = ef_col_idx
                        logger.debug(
                            "Mapped raw extra field '%s' (normalized='%s') to col index %d",
                            raw_field,
                            nf,
                            ef_col_idx,
                        )
                    else:
                        logger.warning(
                            "Extra field '%s' (normalized='%s') not found in header. Skipping.",
                            raw_field,
                            nf,
                        )

            # Output the (possibly unchanged) header line
            yield line
            first_line = False
            continue

        # For data lines
        cols = line.split("\t")
        if gt_idx is None or gt_idx >= len(cols):
            # If we can't find GT or it's out of bounds, pass the line unchanged
            yield line
            continue

        gt_value = cols[gt_idx].strip()
        # Split subfields by snpsift_sep
        genotypes = gt_value.split(snpsift_sep) if gt_value else []

        # If we have extra fields, build a dictionary of field_name => [list of subfield values]
        extra_values_per_field = {}
        if append_extra_fields and extra_field_indices:
            for raw_field, col_idx in extra_field_indices.items():
                if col_idx < len(cols):
                    raw_val = cols[col_idx].strip()
                    extra_values_per_field[raw_field] = (
                        raw_val.split(snpsift_sep) if raw_val else []
                    )
                else:
                    extra_values_per_field[raw_field] = []

        new_gts = []
        for i, genotype_subfield in enumerate(genotypes):
            if i >= len(samples):
                logger.warning(
                    "Line %d => found %d genotype subfields but only %d samples. Truncating.",
                    line_idx,
                    len(genotypes),
                    len(samples),
                )
                break

            # Step 1: Replace phased genotypes
            genotype_subfield = genotype_subfield.replace("|", "/")

            # Step 2: Normalize genotypes with one missing allele (./N or N/.)
            if "/" in genotype_subfield:
                parts = genotype_subfield.split("/")
                if len(parts) == 2:
                    # Handle ./N case: parts[0] is '.', parts[1] consists of one or more digits
                    if parts[0] == "." and parts[1].isdigit():
                        genotype_subfield = "0/" + parts[1]
                    # Handle N/. case: parts[1] is '.', parts[0] consists of one or more digits
                    elif parts[1] == "." and parts[0].isdigit():
                        genotype_subfield = parts[0] + "/0"

            # Step 3: Apply numeric allele replacements
            if "/" in genotype_subfield:
                parts = genotype_subfield.split("/")
                processed_alleles = []
                if len(parts) == 2:  # Ensure it's diploid
                    a1, a2 = parts
                    for allele_part in [a1, a2]:
                        if (
                            allele_part == "." or not allele_part.isdigit()
                        ):  # Keep missing or non-numeric as is
                            processed_alleles.append(allele_part)
                            continue

                        current_allele_val = allele_part
                        # Apply configured regex patterns to the individual allele part
                        for regex_pat, replacement_val in compiled_patterns:
                            current_allele_val = regex_pat.sub(replacement_val, current_allele_val)

                        # After applying patterns, if it's a digit and not '0', make it '1'
                        # This ensures multi-digit numbers like '11' (from '12' via [2-9]->1 map) also become '1'
                        if current_allele_val.isdigit() and current_allele_val != "0":
                            processed_alleles.append("1")
                        elif current_allele_val == "0":  # Explicitly keep '0'
                            processed_alleles.append("0")
                        else:  # If it became non-digit (e.g. map to '.') or was already '.', keep it
                            processed_alleles.append(current_allele_val)
                    genotype_subfield = "/".join(processed_alleles)
                # else: if not 2 parts, do not attempt this logic, leave genotype_subfield as is from previous step
                # else: # Not a diploid genotype, apply regex to whole field (original behavior for non-diploid)
                # This part handles cases like single allele numbers or other formats if they reach here.
                # However, VCF genotypes are typically like N/N, N|N, or N.
                # If it's a single number (e.g. haploid), it should also be binarized.
                original_genotype_subfield_for_non_diploid = genotype_subfield
                for regex_pat, replacement_val in compiled_patterns:
                    if regex_pat.search(genotype_subfield):
                        genotype_subfield = regex_pat.sub(replacement_val, genotype_subfield)
                # If it was a digit, and after map it's still a digit and not 0, make it 1.
                if (
                    original_genotype_subfield_for_non_diploid.isdigit()
                    and genotype_subfield.isdigit()
                    and genotype_subfield != "0"
                ):
                    genotype_subfield = "1"
                # else, it remains as transformed by map, or original if not digit/not 0

            # Step 4: Skip "0/0" or "./. as "non-variant" after all transformations
            if genotype_subfield in ("0/0", "./."):
                continue

            # Build the genotype output
            sample_name = samples[i]
            if append_extra_fields and extra_field_indices:
                # Gather the i-th entry from each extra field
                bits = []
                for raw_field in extra_field_indices.keys():
                    ev_list = extra_values_per_field.get(raw_field, [])
                    val = ev_list[i] if i < len(ev_list) else ""
                    bits.append(val)

                if bits:
                    # e.g. "0/1:100:52,48"
                    genotype_plus_extra = f"{genotype_subfield}{extra_field_delimiter}{extra_field_delimiter.join(bits)}"
                else:
                    genotype_plus_extra = genotype_subfield

                # Format: sampleName(genotype_plus_extra)
                new_gts.append(f"{sample_name}({genotype_plus_extra})")
            else:
                # Just sample(genotype_subfield)
                new_gts.append(f"{sample_name}({genotype_subfield})")

        # Rejoin the new GT column
        if new_gts:
            cols[gt_idx] = separator.join(new_gts)
        else:
            cols[gt_idx] = ""

        yield "\t".join(cols)
