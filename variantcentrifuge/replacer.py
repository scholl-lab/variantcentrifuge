# File: variantcentrifuge/replacer.py

"""
Genotype replacement module.

This module implements the logic of the original replace_gt_with_sample.sh script.
It reads lines from a tab-delimited input, replaces non-"0/0" genotypes in a given
field with sample IDs, and optionally appends genotypes to the sample ID. It can also
output a list of unique samples or add genotype count columns for probands and controls.

If no sample list is provided, this function will simply return the lines unchanged.
"""

def replace_genotypes(lines, sample_file, cfg):
    """
    Replace genotype values with sample IDs and compute genotype counts, if sample info is available.

    Parameters
    ----------
    lines : iterator
        Iterator over lines with extracted fields (tab-delimited).
    sample_file : str
        Path or comma-separated list of samples. If None or empty, expects cfg["sample_list"].
    cfg : dict
        Configuration dictionary containing:
        - append_genotype (bool)
        - gt_field_number (int) (1-based field index)
        - sample_list (str, comma-separated) or sample_file (path)
        - proband_list (str or file)
        - control_list (str or file)
        - include_nocalls (bool)
        - count_genotypes (bool)
        - list_samples (bool) whether to output unique samples instead of modified lines
        - separator (str) the separator to use between multiple samples in GT field (default ";")

    Returns
    -------
    iterator
        Iterator over processed lines or original lines if no sample data is provided.
    """
    append_genotype = cfg.get("append_genotype", False)
    gt_field_number = cfg.get("gt_field_number", 10)  # 1-based index
    sample_list_str = cfg.get("sample_list", "")
    proband_list_str = cfg.get("proband_list", "")
    control_list_str = cfg.get("control_list", "")
    include_nocalls = cfg.get("include_nocalls", False)
    count_genotypes = cfg.get("count_genotypes", True)
    list_samples = cfg.get("list_samples", False)
    separator = cfg.get("separator", ";")

    def read_list(input_str_or_file):
        arr = []
        import os
        if not input_str_or_file:
            return arr
        if os.path.isfile(input_str_or_file):
            with open(input_str_or_file, "r", encoding="utf-8") as f:
                content = f.read().strip()
                if "," in content:
                    arr = [x.strip() for x in content.split(",") if x.strip()]
                else:
                    arr = [x.strip() for x in content.split("\n") if x.strip()]
        else:
            arr = [x.strip() for x in input_str_or_file.split(",") if x.strip()]
        return arr

    # Determine samples
    if sample_list_str and sample_file:
        samples = read_list(sample_list_str)
    elif sample_list_str:
        samples = read_list(sample_list_str)
    elif sample_file:
        samples = read_list(sample_file)
    else:
        # No sample list provided, just return lines unchanged
        for line in lines:
            yield line
        return

    # From here on, we know we have sample info and can proceed
    probands = read_list(proband_list_str) if proband_list_str else samples[:]
    controls = []
    if control_list_str:
        controls = read_list(control_list_str)
    else:
        # Controls = samples - probands
        proband_set = set(probands)
        controls = [s for s in samples if s not in proband_set]

    proband_map = {p: True for p in probands}
    control_map = {c: True for c in controls}

    gt_idx = gt_field_number - 1
    unique_samples = set()
    import sys

    first_line = True
    for line in lines:
        line = line.rstrip("\n")
        if first_line:
            header_fields = line.split("\t")
            if list_samples:
                # If listing samples, we just capture samples at the end
                # No direct line output now
                if count_genotypes:
                    # no changes in header needed here
                    pass
            else:
                # If not list mode, print header and add columns if count_genotypes is True
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

        # Data line
        if list_samples:
            # Collect unique samples info, no immediate output
            pass

        fields = line.split("\t")
        if gt_idx >= len(fields):
            # Invalid line structure, just yield as is if not listing samples
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
        import re
        for i, g in enumerate(genotypes):
            sample = samples[i] if i < len(samples) else None
            if sample is None:
                # More genotypes than samples? Just skip
                continue

            g = g.replace("|", "/")
            # Replace non-1 GT fields with 1 if found
            if re.search(r"[2-9]", g):
                g = re.sub(r"[2-9]", "1", g)
                print(f"Warning: Non-1 GT field detected and replaced with 1 in sample {sample} for genotype {g}", file=sys.stderr)

            if g != "0/0" or (include_nocalls and g == "./."):
                if g != "./.":
                    unique_samples.add(sample)
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
                else:
                    if sample in proband_map:
                        total_count += 1
                    if sample in control_map:
                        control_total_count += 1

            # Replace genotype if non-0/0 and non-./.
            if g not in ["0/0", "./."]:
                if append_genotype:
                    new_gts.append(f"{sample}({g})")
                else:
                    new_gts.append(sample)

        if list_samples:
            # no line output here
            pass
        else:
            # Replace GT field
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
        yield ",".join(sorted(unique_samples))
