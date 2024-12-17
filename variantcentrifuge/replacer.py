# File: variantcentrifuge/replacer.py
# Location: variantcentrifuge/variantcentrifuge/replacer.py

"""
Genotype replacement module.

Now we dynamically find the GT field by searching for a column named "GT" in the header line.
No need for a hardcoded gt_field_number.

Logic:
- If append_genotype is set, output sample(genotype) for variant genotypes.
- If not, just sample.
- Only non-0/0 genotypes (0/1, 1/0, 1/1) and variant output. Skip 0/0 and ./..
- If include_nocalls is true and ./., count in totals but don't output.
- If no variant genotypes found, GT field is empty.
- If we can't find a "GT" column in the header, return lines unchanged.

Requires cfg["sample_list"] for sample names.
If no sample_list, returns lines unchanged.

Enhancement #14:
- Make genotype replacement logic configurable:
  * Introduce a 'genotype_replacement_map' parameter in cfg, mapping regex patterns to their replacements.
  * Default: {r"[2-9]": "1"}
  * Document the default logic and provide options for users to adjust or turn off replacements by providing an empty map.
  * Increase transparency and flexibility in how genotypes are handled.
"""

import logging
import re

logger = logging.getLogger("variantcentrifuge")


def replace_genotypes(lines, cfg):
    """
    Replace genotypes based on user-defined logic and return updated lines.

    The replacement logic is controlled by cfg["genotype_replacement_map"],
    which should be a dictionary mapping regex patterns to replacement strings.

    Default logic (if not provided): {r"[2-9]": "1"}

    Parameters
    ----------
    lines : iterable
        An iterable of lines (strings) representing the extracted TSV data.
    cfg : dict
        Configuration dictionary. Relevant entries:
        - "append_genotype": bool, whether to append genotype (e.g., sample(0/1)) 
          or just sample ID.
        - "sample_list": str, comma-separated list of samples in the order of genotypes.
        - "proband_list": str, comma-separated list of proband samples.
        - "control_list": str, comma-separated list of control samples.
        - "include_nocalls": bool, whether to include no-calls (./.) in count totals.
        - "count_genotypes": bool, whether to append count columns to the output.
        - "list_samples": bool, whether to list unique samples at the end.
        - "separator": str, separator to join replaced genotypes.
        - "genotype_replacement_map": dict, keys are regex patterns, values are 
          replacement strings. Used to transform alleles in genotypes before processing.
          Example: {r"[2-9]": "1"} (default if not provided).

    Yields
    ------
    str
        Modified lines with genotype replacements applied.
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

    def parse_list_from_str(input_str):
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
    controls = []
    if control_list_str:
        controls = parse_list_from_str(control_list_str)
    else:
        # Controls = samples - probands
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
                # No GT column found, return lines unchanged from now on
                yield line
                for l in lines:
                    yield l
                return

            gt_idx = header_fields.index("GT")

            if list_samples:
                # If listing samples, we do not modify the header
                if count_genotypes:
                    # no changes to header, just yield if needed later after processing
                    pass
            else:
                # If not listing samples, print header and add columns if count_genotypes is True
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

            # Apply replacement logic based on genotype_replacement_map
            for pat, repl in compiled_patterns:
                if pat.search(g):
                    g = pat.sub(repl, g)

            # Only output non-0/0 variant genotypes
            if g not in ["0/0", "./."]:
                # This is a variant genotype
                unique_samples.add(sample)

                # Count stats
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
                # If ./., and include_nocalls is true, increment totals
                if include_nocalls and g == "./.":
                    if sample in proband_map:
                        total_count += 1
                    if sample in control_map:
                        control_total_count += 1
                # Do not append anything for these genotypes

        if list_samples:
            # If listing samples, no line output now
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
        # Print the collected unique samples
        yield ",".join(sorted(unique_samples))
