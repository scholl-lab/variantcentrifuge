# File: variantcentrifuge/replacer.py
# Location: variantcentrifuge/variantcentrifuge/replacer.py

"""
Genotype replacement module.

This module dynamically searches for a "GT" column in the header line of the
input TSV and applies configurable genotype replacement logic.

Logic:
- If `append_genotype` is True, output sample(genotype) for variant genotypes;
  else just sample.
- Non-variant genotypes ("0/0", "./.") are skipped for the final column.
- If `include_nocalls` is True, "./." still won't be output but is counted elsewhere.
- If no "GT" column is found, the lines are returned unchanged.
- Requires `cfg["sample_list"]` for sample names.

Optional Extra Fields:
- If `cfg["append_extra_sample_fields"]` is True, extra columns like DP, AD, etc.
  are appended to each genotype in parentheses, e.g. sample(0/1:DP:AD).
- We split each column for each sample using the same SnpSift subfield separator
  (cfg["extract_fields_separator"]).
"""

import logging
import re
from typing import Iterator, Dict, Any

logger = logging.getLogger("variantcentrifuge")


def replace_genotypes(lines: Iterator[str], cfg: Dict[str, Any]) -> Iterator[str]:
    """
    Replace genotypes based on user-defined logic and return updated lines.

    The replacement logic is controlled by `cfg["genotype_replacement_map"]`,
    e.g. {r"[2-9]": "1"} to replace any allele "2".."9" with "1".

    Phased genotypes (0|1, etc.) are normalized by replacing "|" with "/".

    If `cfg["append_extra_sample_fields"]` is True, we also append columns like DP,AD
    next to the genotype, separated by `cfg["extra_sample_field_delimiter"]`.
    E.g.: sample(0/1:100:52,48)

    Debug logs are added to help track potential issues with missing fields or columns.
    """
    append_genotype = cfg.get("append_genotype", True)
    sample_list_str = cfg.get("sample_list", "")
    proband_list_str = cfg.get("proband_list", "")
    control_list_str = cfg.get("control_list", "")
    include_nocalls = cfg.get("include_nocalls", False)
    list_samples = cfg.get("list_samples", False)
    separator = cfg.get("separator", ";")
    genotype_replacement_map = cfg.get("genotype_replacement_map", {r"[2-9]": "1"})

    # Extra fields config
    append_extra_fields = cfg.get("append_extra_sample_fields", False)
    extra_sample_fields = cfg.get("extra_sample_fields", [])
    extra_field_delimiter = cfg.get("extra_sample_field_delimiter", ":")

    # The real multi-sample subfield separator used by SnpSift for columns
    snpsift_sep = cfg.get("extract_fields_separator", ":")

    logger.debug("replace_genotypes called with config: %s", cfg)
    logger.debug("append_genotype=%s, extra_fields=%s, snpsift_sep='%s'",
                 append_genotype, extra_sample_fields, snpsift_sep)

    def parse_list_from_str(s: str) -> list:
        """Split a comma-separated sample_list or proband_list string into a list."""
        if not s:
            return []
        return [x.strip() for x in s.split(",") if x.strip()]

    samples = parse_list_from_str(sample_list_str)
    logger.debug("Parsed samples: %s", samples)
    if not samples:
        logger.warning("No samples found in cfg[sample_list]. Returning lines unchanged.")
        # Just yield lines and return
        for line in lines:
            yield line
        return

    # Optional usage for proband/control, though not strictly used here
    probands = parse_list_from_str(proband_list_str) if proband_list_str else samples[:]
    if control_list_str:
        controls = parse_list_from_str(control_list_str)
        logger.debug("Parsed controls: %s", controls)
    else:
        prob_set = set(probands)
        controls = [s for s in samples if s not in prob_set]
        logger.debug("Inferred controls: %s", controls)

    # Compile genotype replacement regex patterns
    compiled_patterns = []
    for pat, repl in genotype_replacement_map.items():
        compiled = re.compile(pat)
        compiled_patterns.append((compiled, repl))

    unique_samples = set()
    first_line = True
    gt_idx = None
    # We'll keep track of extra field column indices:
    extra_field_indices = {}

    for line_idx, line in enumerate(lines, start=1):
        line = line.rstrip("\n")
        if first_line:
            header_cols = line.split("\t")
            logger.debug("Header columns at line_idx=%d: %s", line_idx, header_cols)

            if "GT" not in header_cols:
                logger.warning("No GT column found in header. Returning lines unmodified.")
                yield line
                for rest in lines:
                    yield rest
                return

            gt_idx = header_cols.index("GT")
            logger.debug("Found GT column index: %d", gt_idx)

            # Identify requested extra fields
            if append_extra_fields and extra_sample_fields:
                for ef_name in extra_sample_fields:
                    if ef_name in header_cols:
                        idx_ = header_cols.index(ef_name)
                        extra_field_indices[ef_name] = idx_
                        logger.debug("Found extra sample field '%s' at col index %d", ef_name, idx_)
                    else:
                        logger.warning("Extra field '%s' not found in header. Skipping it.", ef_name)

            if list_samples:
                logger.debug("list_samples=True => not modifying header line.")
                yield line
            else:
                yield line
            first_line = False
            continue

        # Subsequent data lines
        cols = line.split("\t")
        if gt_idx >= len(cols):
            logger.warning("Line %d has no GT value (index out of range). Skipping line => %s", line_idx, line)
            if not list_samples:
                yield line
            continue

        gt_value = cols[gt_idx]
        # Split subfields by snpsift_sep (commonly ",")
        genotypes = gt_value.split(snpsift_sep)

        # Also split extra fields if applicable
        extra_values_per_field = {}
        if append_extra_fields and extra_field_indices:
            for ef_name, ef_col_idx in extra_field_indices.items():
                if ef_col_idx < len(cols):
                    raw_val = cols[ef_col_idx]
                    extra_values = raw_val.split(snpsift_sep)
                else:
                    extra_values = []
                extra_values_per_field[ef_name] = extra_values
                logger.debug("Line %d => Extra field '%s' => values: %s", line_idx, ef_name, extra_values)

        # Build the new GT column
        new_gts = []
        for i, genotype in enumerate(genotypes):
            if i >= len(samples):
                # More genotype subfields than sample names => skip extras
                logger.warning("Line %d => found %d genotype subfields but only %d samples. Truncating.",
                               line_idx, len(genotypes), len(samples))
                break

            genotype = genotype.replace("|", "/")
            # Apply regex replacements
            for (compiled_pat, repl_str) in compiled_patterns:
                if compiled_pat.search(genotype):
                    genotype = compiled_pat.sub(repl_str, genotype)

            # Skip "0/0" or "./."
            if genotype in ["0/0", "./."]:
                continue

            sample_name = samples[i]
            unique_samples.add(sample_name)
            # If we want "sample(0/1:DP:AD)"
            if append_genotype:
                if append_extra_fields and extra_values_per_field:
                    # gather the i-th entry from each extra field
                    bits = []
                    for ef_ in extra_sample_fields:
                        # e.g. bits for DP, AD, etc.
                        the_list = extra_values_per_field.get(ef_, [])
                        val = the_list[i] if i < len(the_list) else ""
                        bits.append(val)
                    if bits:
                        # e.g. "0/1:100:52,48"
                        genotype_with_extras = f"{genotype}{extra_field_delimiter}{extra_field_delimiter.join(bits)}"
                    else:
                        genotype_with_extras = genotype
                    new_gts.append(f"{sample_name}({genotype_with_extras})")
                else:
                    # Just sample(0/1)
                    new_gts.append(f"{sample_name}({genotype})")
            else:
                # If we do not want genotype appended, just sample
                new_gts.append(sample_name)

        if list_samples:
            # We won't modify columns, just keep track
            continue
        else:
            if new_gts:
                # Join with e.g. ";"
                cols[gt_idx] = separator.join(new_gts)
            else:
                # no variants => empty
                cols[gt_idx] = ""
            yield "\t".join(cols)

    if list_samples:
        # Output collected variant samples
        logger.debug("list_samples=True => final unique samples: %s", sorted(unique_samples))
        yield ",".join(sorted(unique_samples))
