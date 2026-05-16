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
