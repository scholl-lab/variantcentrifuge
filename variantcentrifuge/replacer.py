# File: variantcentrifuge/replacer.py
# Location: variantcentrifuge/variantcentrifuge/replacer.py

"""
Genotype replacement module.

This module dynamically searches for a "GT" column in the header line of the input TSV
and applies configurable genotype replacement logic.

Logic:
- If `append_genotype` is True, output sample(genotype) for variant genotypes; else just sample.
- Only non-"0/0" genotypes ("0/1", "1/0", "1/1") and their phased counterparts ("0|1", "1|0", "1|1") 
  are output as variants. "0/0" and "./." are skipped.
- If `include_nocalls` is True, "./." contributes to totals but is not output.
- If no variant genotypes are found, the GT field is left empty.
- If no "GT" column is found, the lines are returned unchanged.
- Requires `cfg["sample_list"]` for sample names; if not provided, returns lines unchanged.

Enhancement #14:
- Genotype replacement logic can be configured via `cfg["genotype_replacement_map"]`.
- Default: {r"[2-9]": "1"} replaces any allele >1 with '1'.
- Users can adjust this map or provide an empty map to skip replacements.

Enhancement #30 (bug fix):
- Previously, phased genotypes (0|1, 1|0, etc.) were not handled. Now, any "|" separators in genotypes 
  are replaced with "/" before applying the replacement logic. This ensures phased genotypes are treated 
  consistently and recognized as variants when appropriate.

Relevant `cfg` parameters:
- "append_genotype": bool
- "sample_list": str (comma-separated)
- "proband_list": str (comma-separated)
- "control_list": str (comma-separated)
- "include_nocalls": bool
- "list_samples": bool
- "separator": str
- "genotype_replacement_map": dict (e.g., {r"[2-9]": "1"})

If `list_samples` is True, instead of modifying lines, the module will collect and output unique variant samples at the end.

Note:
This module no longer performs any counting of genotypes. All counting logic has been removed.
Counting of variants, allele counts, and homozygous counts should be handled elsewhere (e.g., in `helpers.py`).
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

    Phased genotypes (e.g., 0|1, 1|0) are normalized to slash-based genotypes (0/1, 1/0)
    before applying the replacement logic. This ensures consistent handling of both 
    phased and unphased genotypes.

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
            Whether to include no-calls (./.) in processing totals 
            (no direct output impact here).
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
                # If listing samples, do not modify the header
                yield line
            else:
                # Just yield the original header line
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

        new_gts = []

        for i, g in enumerate(genotypes):
            if i >= len(samples):
                # More genotypes than samples?
                continue
            sample = samples[i]

            # Normalize phased genotypes by replacing '|' with '/'
            g = g.replace("|", "/")

            # Apply replacement logic
            for pat, repl in compiled_patterns:
                if pat.search(g):
                    g = pat.sub(repl, g)

            # Only output non-0/0 variant genotypes
            if g not in ["0/0", "./."]:
                # Variant genotype
                unique_samples.add(sample)
                if append_genotype:
                    new_gts.append(f"{sample}({g})")
                else:
                    new_gts.append(sample)
            else:
                # g is "0/0" or "./." - skip variant output
                pass

        if list_samples:
            # If listing samples, accumulate only
            pass
        else:
            # Replace GT field with only variant samples
            if new_gts:
                fields[gt_idx] = separator.join(new_gts)
            else:
                fields[gt_idx] = ""
            yield "\t".join(fields)

    if list_samples:
        # If listing samples, output them at the end
        yield ",".join(sorted(unique_samples))


#
# Tests for phased genotype handling
# These tests are embedded here for demonstration purposes.
# In a real scenario, place them in a separate test file.
#

def test_replace_genotypes_phased():
    """
    Test that phased genotypes (with '|') are correctly handled and normalized to '/'.
    """
    cfg = {
        "append_genotype": True,
        "sample_list": "sample1,sample2",
        "include_nocalls": False,
        "list_samples": False,
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"}
    }

    input_lines = [
        "CHROM\tPOS\tREF\tALT\tGT",
        "1\t100\tA\tT\tsample1(0|1),sample2(1|0)",
        "1\t200\tG\tC\tsample1(0|1),sample2(0/0)"
    ]

    # Expected:
    # For line 1: 0|1 and 1|0 become 0/1 and 1/0, both are variants -> "sample1(0/1);sample2(1/0)"
    # For line 2: 0|1 becomes 0/1 (variant), 0/0 is not output
    # -> "sample1(0/1)"

    expected_output = [
        "CHROM\tPOS\tREF\tALT\tGT",
        "1\t100\tA\tT\tsample1(0/1);sample2(1/0)",
        "1\t200\tG\tC\tsample1(0/1)"
    ]

    output = list(replace_genotypes(iter(input_lines), cfg))
    assert output == expected_output, f"Phased genotype handling failed. Got: {output}, Expected: {expected_output}"

def test_replace_genotypes_phased_with_replacements():
    """
    Test phased genotypes with replacements. For example, a genotype like 2|1 
    should be replaced using the provided map, turning '2' into '1'.
    """
    cfg = {
        "append_genotype": True,
        "sample_list": "sample1",
        "include_nocalls": False,
        "list_samples": False,
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"}
    }

    input_lines = [
        "CHROM\tPOS\tREF\tALT\tGT",
        "1\t300\tA\tT\tsample1(2|1)"
    ]

    # After normalization: 2|1 -> 2/1, after replacement: 2 -> 1, so becomes 1/1
    # 1/1 is a variant, so output "sample1(1/1)"

    expected_output = [
        "CHROM\tPOS\tREF\tALT\tGT",
        "1\t300\tA\tT\tsample1(1/1)"
    ]

    output = list(replace_genotypes(iter(input_lines), cfg))
    assert output == expected_output, f"Phased genotype replacement failed. Got: {output}, Expected: {expected_output}"
