import logging

from variantcentrifuge.vcf_eff_one_per_line import split_vcf_effects


def _info_map(vcf_line: str) -> dict[str, str | None]:
    info = vcf_line.split("\t")[7]
    result: dict[str, str | None] = {}
    for part in info.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            result[key] = value
        else:
            result[part] = None
    return result


def _split_rows_by_gene(vcf_line: str) -> dict[str, dict[str, str | None]]:
    rows = split_vcf_effects(vcf_line)
    by_gene = {}
    for row in rows:
        info = _info_map(row)
        ann_parts = (info["ANN"] or "").split("|")
        by_gene[ann_parts[3]] = info
    return by_gene


def test_split_ann_prunes_lof_nmd_to_matching_gene() -> None:
    line = (
        "chr1\t100\t.\tA\tT\t.\tPASS\t"
        "AC=1;"
        "ANN="
        "T|splice_donor_variant|HIGH|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||,"
        "T|downstream_gene_variant|MODIFIER|GENE_B|ENSGB|transcript|NM_b|protein_coding||||||||;"
        "LOF=(GENE_A|ENSGA|1|0.95);"
        "NMD=(GENE_A|ENSGA|1|0.91)"
        "\tGT\t0/1"
    )

    rows_by_gene = _split_rows_by_gene(line)

    assert rows_by_gene["GENE_A"]["LOF"] == "(GENE_A|ENSGA|1|0.95)"
    assert rows_by_gene["GENE_A"]["NMD"] == "(GENE_A|ENSGA|1|0.91)"
    assert "LOF" not in rows_by_gene["GENE_B"]
    assert "NMD" not in rows_by_gene["GENE_B"]
    assert rows_by_gene["GENE_A"]["AC"] == "1"
    assert rows_by_gene["GENE_B"]["AC"] == "1"


def test_split_ann_matches_lof_nmd_by_gene_id_when_symbol_differs() -> None:
    line = (
        "chr1\t100\t.\tA\tT\t.\tPASS\t"
        "ANN="
        "T|splice_donor_variant|HIGH|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||,"
        "T|missense_variant|MODERATE|GENE_B_ALIAS|ENSGB|transcript|NM_b|protein_coding||||||||;"
        "LOF=(GENE_B_CANONICAL|ENSGB|1|0.95);"
        "NMD=(GENE_B_CANONICAL|ENSGB|1|0.91)"
    )

    rows_by_gene = _split_rows_by_gene(line)

    assert "LOF" not in rows_by_gene["GENE_A"]
    assert "NMD" not in rows_by_gene["GENE_A"]
    assert rows_by_gene["GENE_B_ALIAS"]["LOF"] == "(GENE_B_CANONICAL|ENSGB|1|0.95)"
    assert rows_by_gene["GENE_B_ALIAS"]["NMD"] == "(GENE_B_CANONICAL|ENSGB|1|0.91)"


def test_split_ann_drops_malformed_lof_nmd_entries_and_warns(caplog) -> None:
    line = (
        "chr1\t100\t.\tA\tT\t.\tPASS\t"
        "ANN="
        "T|splice_donor_variant|HIGH|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||,"
        "T|downstream_gene_variant|MODIFIER|GENE_B|ENSGB|transcript|NM_b|protein_coding||||||||;"
        "LOF=(GENE_A|ENSGA|1|0.95),malformed;"
        "NMD=also_malformed"
    )

    with caplog.at_level(logging.WARNING, logger="variantcentrifuge.vcf_eff_one_per_line"):
        rows_by_gene = _split_rows_by_gene(line)

    assert rows_by_gene["GENE_A"]["LOF"] == "(GENE_A|ENSGA|1|0.95)"
    assert "NMD" not in rows_by_gene["GENE_A"]
    assert "LOF" not in rows_by_gene["GENE_B"]
    assert "NMD" not in rows_by_gene["GENE_B"]
    assert "Dropping malformed LOF entry" in caplog.text
    assert "Dropping malformed NMD entry" in caplog.text


def test_split_eff_preserves_lof_nmd_without_ann_pruning() -> None:
    line = (
        "chr1\t100\t.\tA\tT\t.\tPASS\t"
        "AC=1;"
        "EFF=STOP_GAINED(HIGH|MISSENSE|c.1A>T|p.Lys1*|1|GENE_A|protein_coding|CODING|NM_a|1),"
        "DOWNSTREAM(MODIFIER|||||GENE_B|protein_coding|CODING|NM_b|1);"
        "LOF=(GENE_A|ENSGA|1|0.95);"
        "NMD=(GENE_A|ENSGA|1|0.91)"
    )

    rows = split_vcf_effects(line)

    assert len(rows) == 2
    for row in rows:
        info = _info_map(row)
        assert info["LOF"] == "(GENE_A|ENSGA|1|0.95)"
        assert info["NMD"] == "(GENE_A|ENSGA|1|0.91)"


def test_multiallelic_somatic_record_decomposes_to_biallelic_rows_with_scalar_af_and_matching_ann(
    tmp_path,
) -> None:
    import pandas as pd
    import pytest

    from variantcentrifuge.extractor import extract_fields_bcftools, parse_ann_subfields
    from variantcentrifuge.filters import filter_dataframe_with_query

    line = (
        "chr1\t100\t.\tC\tA,G\t.\tPASS\t"
        "AC=1,2;AF=0.01,0.30;DP=121;"
        "ANN="
        "A|synonymous_variant|LOW|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||,"
        "G|missense_variant|MODERATE|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||"
        "\tGT:AD:AF:DP\t0/1:90,1,30:0.01,0.30:121\t0/0:100,0,0:0.00,0.00:100"
    )

    rows = split_vcf_effects(line)
    assert len(rows) == 2

    cols_a = rows[0].split("\t")
    cols_g = rows[1].split("\t")

    assert cols_a[4] == "A"
    assert cols_a[9] == "0/1:90,1:0.01:121"
    assert cols_a[10] == "0/0:100,0:0.00:100"
    assert _info_map(rows[0])["ANN"] == (
        "A|synonymous_variant|LOW|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||"
    )

    assert cols_g[4] == "G"
    assert cols_g[9] == "0/1:90,30:0.30:121"
    assert cols_g[10] == "0/0:100,0:0.00:100"
    assert _info_map(rows[1])["ANN"] == (
        "G|missense_variant|MODERATE|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||"
    )

    # Verify parse_ann_subfields + filter_dataframe_with_query on decomposed rows
    df = pd.DataFrame(
        [
            {
                "CHROM": cols[0],
                "POS": cols[1],
                "REF": cols[3],
                "ALT": cols[4],
                "ANN": _info_map(r)["ANN"],
                "GEN_0__AF": cols[9].split(":")[2],
            }
            for r, cols in ((rows[0], cols_a), (rows[1], cols_g))
        ]
    )
    df = parse_ann_subfields(df, ["ANN[0].ALLELE", "ANN[0].EFFECT", "ANN[0].IMPACT"])
    df = df.rename(columns=lambda c: c.replace("[", "_").replace("].", "__"))
    filtered = filter_dataframe_with_query(
        df, "(GEN_0__AF >= 0.05) & (ANN_0__IMPACT == 'MODERATE')"
    )
    assert len(filtered) == 1
    assert filtered.iloc[0]["ALT"] == "G"
    assert float(filtered.iloc[0]["GEN_0__AF"]) == 0.30
    assert filtered.iloc[0]["ANN_0__EFFECT"] == "missense_variant"

    # Also verify reject_multiallelic guard in extract_fields_bcftools
    raw_vcf = tmp_path / "multiallelic.vcf"
    raw_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n" + line + "\n"
    )
    with pytest.raises(ValueError, match="Unsplit multiallelic"):
        extract_fields_bcftools(
            str(raw_vcf),
            "CHROM POS REF ALT GEN[*].AF",
            {"reject_multiallelic": True},
            str(tmp_path / "out.tsv"),
        )
