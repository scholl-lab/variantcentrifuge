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
"""

def replace_genotypes(lines, cfg):
    append_genotype = cfg.get("append_genotype", True)
    sample_list_str = cfg.get("sample_list", "")
    proband_list_str = cfg.get("proband_list", "")
    control_list_str = cfg.get("control_list", "")
    include_nocalls = cfg.get("include_nocalls", False)
    count_genotypes = cfg.get("count_genotypes", True)
    list_samples = cfg.get("list_samples", False)
    separator = cfg.get("separator", ";")

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

    import re

    for line in lines:
        line = line.rstrip("\n")
        if first_line:
            header_fields = line.split("\t")
            # Find "GT" column
            # After extractFields and header fixes, we assume there's a column named "GT".
            # If not found, return lines unchanged.
            if "GT" not in header_fields:
                # No GT column found, return lines unchanged from now on
                yield line
                for l in lines:
                    yield l
                return

            gt_idx = header_fields.index("GT")

            if list_samples:
                # If listing samples, we just capture samples at the end
                # No direct line output now.
                if count_genotypes:
                    # no changes to header
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
            # Replace any digits 2-9 with '1'
            if re.search(r"[2-9]", g):
                g = re.sub(r"[2-9]", "1", g)
                print(f"Warning: Non-1 GT field detected and replaced with 1 in sample {sample} for genotype {g}", file=sys.stderr)

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
