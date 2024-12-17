# File: variantcentrifuge/replacer.py
# Location: variantcentrifuge/variantcentrifuge/replacer.py

"""
Genotype replacement module.

This module dynamically searches for a "GT" column in the header line of the input TSV
and applies configurable genotype replacement logic.

Logic:
- If `append_genotype` is True, output sample(genotype) for variant genotypes; else just sample.
- Only non-"0/0" genotypes ("0/1", "1/0", "1/1") are output as variants. "0/0" and "./." are skipped.
- If `include_nocalls` is True, "./." contributes to totals but is not output.
- If no variant genotypes are found, the GT field is left empty.
- If no "GT" column is found, the lines are returned unchanged.
- Requires `cfg["sample_list"]` for sample names; if not provided, returns lines unchanged.

Enhancement #14:
- Genotype replacement logic can be configured via `cfg["genotype_replacement_map"]`.
- Default: {r"[2-9]": "1"} replaces any allele >1 with '1'.
- Users can adjust this map or provide an empty map to skip replacements.

Relevant `cfg` parameters:
- "append_genotype": bool  
- "sample_list": str (comma-separated)  
- "proband_list": str (comma-separated)  
- "control_list": str (comma-separated)  
- "include_nocalls": bool  
- "count_genotypes": bool  
- "list_samples": bool  
- "separator": str  
- "genotype_replacement_map": dict (e.g., {r"[2-9]": "1"})  

If `list_samples` is True, instead of modifying lines, the module will collect and output unique variant samples at the end.
"""

import logging
import re
from typing import Iterator, Dict, Any

logger = logging.getLogger("variantcentrifuge")


def replace_genotypes(lines: Iterator[str], cfg: Dict[str, Any]) -> Iterator[str]:
    """
    Replace genotypes based on user-defined logic and return updated lines.

    The replacement logic is controlled by `cfg["genotype_replacement_map"]`,
    which should be a dictionary mapping regex patterns to replacement strings.

    Default logic (if not provided): {r"[2-9]": "1"}

    Parameters
    ----------
    lines : Iterator[str]
        An iterable of lines (strings) representing the extracted TSV data.
    cfg : dict
        Configuration dictionary. Relevant entries:
        - "append_genotype": bool
            Whether to append genotype (e.g., sample(0/1)) or just sample ID.
        - "sample_list": str
            Comma-separated list of samples in the order of genotypes.
        - "proband_list": str
            Comma-separated list of samples considered probands.
        - "control_list": str
            Comma-separated list of samples considered controls.
        - "include_nocalls": bool
            Whether to include no-calls (./.) in count totals.
        - "count_genotypes": bool
            Whether to append genotype count columns to the output.
        - "list_samples": bool
            Whether to list unique variant samples at the end instead of outputting variants.
        - "separator": str
            Separator to join replaced genotypes.
        - "genotype_replacement_map": dict
            Keys are regex patterns, values are replacement strings. Used to transform alleles.

    Yields
    ------
    str
        Modified lines with genotype replacements applied. If `list_samples` is True,
        yields a final line listing unique variant samples.
    """
    append_genotype = cfg.get("append_genotype", True)
    sample_list_str = cfg.get("sample_list", "")
    proband_list_str = cfg.get("proband_list", "")
    control_list_str = cfg.get("control_list", "")
    include_nocalls = cfg.get("include_nocalls", False)
    count_genotypes = cfg.get("count_genotypes", True)
    list_samples = cfg.get("list_samples", False)
    separator = cfg.get("separator", ";")
    genotype_replacement_map = cfg.get("genotype_replacement_map", {r"[2-9]": "1"})

    def parse_list_from_str(input_str: str):
        if not input_str:
            return []
        return [x.strip() for x in input_str.split(",") if x.strip()]

    samples = parse_list_from_str(sample_list_str)
    if not samples:
        # No sample list provided, just return lines unchanged
        for line in lines:
            yield line
        return

    probands = parse_list_from_str(proband_list_str) if proband_list_str else samples[:]
    if control_list_str:
        controls = parse_list_from_str(control_list_str)
    else:
        # Controls = samples not in probands
        proband_set = set(probands)
        controls = [s for s in samples if s not in proband_set]

    proband_map = {p: True for p in probands}
    control_map = {c: True for c in controls}

    unique_samples = set()
    gt_idx = None
    first_line = True

    # Compile regex patterns for performance
    compiled_patterns = [(re.compile(pat), repl) for pat, repl in genotype_replacement_map.items()]

    for line in lines:
        line = line.rstrip("\n")
        if first_line:
            header_fields = line.split("\t")
            # Find "GT" column
            if "GT" not in header_fields:
                # No GT column found, return lines unchanged
                yield line
                for l in lines:
                    yield l
                return

            gt_idx = header_fields.index("GT")

            if list_samples:
                # If listing samples, we don't modify the header,
                # counts might not be appended since we only output samples at the end.
                if count_genotypes:
                    # Still just yield the header as-is here.
                    pass
            else:
                # If not listing samples, output header (and append columns if count_genotypes is True)
                if count_genotypes:
                    if len(controls) > 0:
                        new_header = header_fields + [
                            "proband_count", "proband_variant_count", "proband_allele_count",
                            "control_count", "control_variant_count", "control_allele_count"
                        ]
                    else:
                        new_header = header_fields + [
                            "proband_count", "proband_variant_count", "proband_allele_count"
                        ]
                    yield "\t".join(new_header)
                else:
                    yield line
            first_line = False
            continue

        fields = line.split("\t")
        if gt_idx >= len(fields):
            # Invalid line structure
            if not list_samples:
                yield line
            continue

        gt_field = fields[gt_idx]
        genotypes = gt_field.split(",")

        het_count = 0
        hom_count = 0
        variant_count = 0
        total_count = len(probands)

        control_het_count = 0
        control_hom_count = 0
        control_variant_count = 0
        control_total_count = len(controls)

        new_gts = []

        for i, g in enumerate(genotypes):
            if i >= len(samples):
                # More genotypes than samples?
                continue
            sample = samples[i]

            g = g.replace("|", "/")

            # Apply replacement logic
            for pat, repl in compiled_patterns:
                if pat.search(g):
                    g = pat.sub(repl, g)

            # Only output non-0/0 variant genotypes
            if g not in ["0/0", "./."]:
                # Variant genotype
                unique_samples.add(sample)

                # Count stats for probands/controls
                if sample in proband_map:
                    if g in ["0/1", "1/0"]:
                        het_count += 1
                    elif g == "1/1":
                        hom_count += 1
                    variant_count += 1

                if sample in control_map:
                    if g in ["0/1", "1/0"]:
                        control_het_count += 1
                    elif g == "1/1":
                        control_hom_count += 1
                    control_variant_count += 1

                # Append sample or sample(genotype)
                if append_genotype:
                    new_gts.append(f"{sample}({g})")
                else:
                    new_gts.append(sample)
            else:
                # g is "0/0" or "./."
                if include_nocalls and g == "./.":
                    # No-call counts towards total but not output
                    if sample in proband_map:
                        total_count += 1
                    if sample in control_map:
                        control_total_count += 1

        if list_samples:
            # If listing samples, don't output line now
            pass
        else:
            # Replace GT field with only variant samples
            if new_gts:
                fields[gt_idx] = separator.join(new_gts)
            else:
                fields[gt_idx] = ""

            if count_genotypes:
                if len(controls) > 0:
                    fields += [
                        str(total_count),
                        str(variant_count),
                        str(het_count + (2 * hom_count)),
                        str(control_total_count),
                        str(control_variant_count),
                        str(control_het_count + (2 * control_hom_count))
                    ]
                else:
                    fields += [
                        str(total_count),
                        str(variant_count),
                        str(het_count + (2 * hom_count))
                    ]
            yield "\t".join(fields)

    if list_samples:
        # If listing samples, output them at the end
        yield ",".join(sorted(unique_samples))
