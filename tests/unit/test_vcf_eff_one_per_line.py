import gzip
import logging

from variantcentrifuge.vcf_eff_one_per_line import (
    restrict_vcf_ann_to_alt,
    split_vcf_effects,
)


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


def test_bcftools_norm_decomposed_rows_prune_shared_ann_and_extract_matching_allele(
    tmp_path, monkeypatch
) -> None:
    import pandas as pd

    from variantcentrifuge.extractor import parse_ann_subfields
    from variantcentrifuge.filters import extract_variants, filter_dataframe_with_query

    shared_ann = (
        "ANN="
        "A|synonymous_variant|LOW|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||,"
        "G|missense_variant|MODERATE|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||"
    )
    # Records as produced by `bcftools norm -m -any` from C -> A,G
    decomposed_a = (
        f"chr1\t100\t.\tC\tA\t.\tPASS\tAC=1;AF=0.01;DP=121;{shared_ann}"
        "\tGT:AD:AF:DP\t0/1:90,1:0.01:121\t0/0:100,0:0.00:100"
    )
    decomposed_g = (
        f"chr1\t100\t.\tC\tG\t.\tPASS\tAC=2;AF=0.30;DP=121;{shared_ann}"
        "\tGT:AD:AF:DP\t0/1:90,30:0.30:121\t0/0:100,0:0.00:100"
    )

    rows_a = split_vcf_effects(decomposed_a)
    rows_g = split_vcf_effects(decomposed_g)
    assert len(rows_a) == 1
    assert len(rows_g) == 1

    assert _info_map(rows_a[0])["ANN"] == (
        "A|synonymous_variant|LOW|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||"
    )
    assert _info_map(rows_g[0])["ANN"] == (
        "G|missense_variant|MODERATE|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||"
    )

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
            for r in (rows_a[0], rows_g[0])
            for cols in (r.split("\t"),)
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

    # Verify extract_variants -> apply_snpsift_filter with annotation splitting disabled
    # restricts ANN (and LOF/NMD) to the decomposed ALT before SnpSift evaluates HIGH
    import gzip
    from unittest.mock import Mock

    from variantcentrifuge.filters import apply_snpsift_filter

    shared_high_ann = (
        "ANN="
        "A|synonymous_variant|LOW|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||,"
        "G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||;"
        "LOF=(KCNJ5|ENSG1|1|1.0)"
    )
    raw_a = (
        f"chr1\t100\t.\tC\tA\t.\tPASS\tAC=1;AF=0.01;DP=121;{shared_high_ann}"
        "\tGT:AD:AF:DP\t0/1:90,1:0.01:121"
    )
    raw_g = (
        f"chr1\t100\t.\tC\tG\t.\tPASS\tAC=2;AF=0.30;DP=121;{shared_high_ann}"
        "\tGT:AD:AF:DP\t0/1:90,30:0.30:121"
    )

    cmds: list[list[str]] = []

    def fake_run(cmd: list[str], output_file: str | None = None, return_result: bool = False):
        cmds.append(cmd)
        if cmd[:2] == ["bcftools", "norm"] and "-o" in cmd:
            out_path = cmd[cmd.index("-o") + 1]
            with gzip.open(out_path, "wt") as fh:
                fh.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                )
                fh.write(raw_a + "\n" + raw_g + "\n")
        elif cmd[0] == "SnpSift" and output_file:
            in_path = cmd[-1]
            opener = gzip.open if in_path.endswith(".gz") else open
            kept: list[str] = []
            with opener(in_path, "rt") as fh:  # type: ignore[call-overload]
                for ln in fh:
                    if ln.startswith("#") or "|HIGH|" in ln:
                        kept.append(ln)
            with open(output_file, "w") as out_fh:
                out_fh.writelines(kept)
        elif cmd[0] == "bgzip" and output_file:
            with open(cmd[-1]) as in_fh, gzip.open(output_file, "wt") as out_fh:
                out_fh.write(in_fh.read())
        res = Mock()
        res.stderr = ""
        return res if return_result else None

    monkeypatch.setattr("variantcentrifuge.filters.run_command", fake_run)
    out_vcf = tmp_path / "extracted.vcf.gz"
    filtered_vcf = tmp_path / "filtered.vcf.gz"
    extract_variants("input.vcf.gz", "regions.bed", {"threads": 2}, str(out_vcf))
    assert len(cmds) == 2
    assert cmds[0][:2] == ["bcftools", "view"]
    assert cmds[1][:2] == ["bcftools", "norm"]
    assert "-W" in cmds[1] and "-m" in cmds[1] and "-any" in cmds[1]

    # Without splitting annotations, apply_snpsift_filter evaluates only the decomposed allele's ANN
    # and LOF/NMD is removed from the same-gene LOW allele (A) while retained on HIGH allele (G)
    with gzip.open(out_vcf, "rt") as fh:
        pruned_records = [ln for ln in fh.read().splitlines() if ln and not ln.startswith("#")]
    assert len(pruned_records) == 2
    assert "LOF=" not in pruned_records[0]
    assert "LOF=(KCNJ5|ENSG1|1|1.0)" in pruned_records[1]

    # Verify rewritten out_vcf is valid BGZF with rebuilt index supporting region queries
    import shutil
    import subprocess

    if shutil.which("tabix"):
        q = subprocess.run(
            ["tabix", str(out_vcf), "chr1:100-100"],
            capture_output=True,
            text=True,
            check=True,
        )
        assert len(q.stdout.strip().splitlines()) == 2

    apply_snpsift_filter(
        str(out_vcf),
        "(ANN[ANY].IMPACT has 'HIGH') | ((exists LOF[*].PERC) & (LOF[*].PERC > 0.9))",
        {"threads": 2},
        str(filtered_vcf),
    )
    with gzip.open(filtered_vcf, "rt") as fh:
        kept_records = [ln for ln in fh.read().splitlines() if ln and not ln.startswith("#")]
    assert len(kept_records) == 1
    assert kept_records[0].split("\t")[4] == "G"
    assert "A|synonymous_variant" not in kept_records[0]

    # Unannotated ALT record: ALT=A row carrying only foreign G stop_gained annotation
    from variantcentrifuge.vcf_eff_one_per_line import restrict_line_ann_to_alt

    foreign_only_a = (
        "chr1\t100\t.\tC\tA\t.\tPASS\tAC=1;"
        "ANN=G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||;"
        "LOF=(KCNJ5|ENSG1|1|1.0);NMD=(KCNJ5|ENSG1|1|1.0)\tGT\t0/1"
    )
    restricted_a = restrict_line_ann_to_alt(foreign_only_a)
    assert "ANN=" not in restricted_a
    assert "LOF=" not in restricted_a
    assert "NMD=" not in restricted_a

    split_foreign_a = split_vcf_effects(foreign_only_a)
    assert len(split_foreign_a) == 1
    assert "ANN=" not in split_foreign_a[0]
    assert "LOF=" not in split_foreign_a[0]
    assert "NMD=" not in split_foreign_a[0]

    unann_df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1"],
            "POS": ["100", "100"],
            "REF": ["C", "C"],
            "ALT": ["A", "G"],
            "ANN": [
                "G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||",
                "G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||",
            ],
        }
    )
    parsed_unann = parse_ann_subfields(unann_df, ["ANN[0].GENE", "ANN[0].EFFECT", "ANN[0].IMPACT"])
    assert parsed_unann.loc[0, "ANN[0].GENE"] == "NA"
    assert parsed_unann.loc[0, "ANN[0].EFFECT"] == "NA"
    assert parsed_unann.loc[0, "ANN[0].IMPACT"] == "NA"
    assert parsed_unann.loc[1, "ANN[0].GENE"] == "KCNJ5"
    assert parsed_unann.loc[1, "ANN[0].IMPACT"] == "HIGH"


def test_restrict_vcf_ann_streams_without_readlines(tmp_path, monkeypatch):
    """Verify that restrict_vcf_ann_to_alt works when readlines() is explicitly disabled."""
    from variantcentrifuge.vcf_eff_one_per_line import restrict_vcf_ann_to_alt

    vcf_file = tmp_path / "test.vcf"
    raw_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tC\tA\t.\tPASS\tANN=G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||\tGT\t0/1\n"
    )
    vcf_file.write_text(raw_content)

    import builtins

    orig_open = builtins.open

    class NoReadlinesFile:
        def __init__(self, path, mode="r", **kwargs):
            self._f = orig_open(path, mode, **kwargs)

        def __iter__(self):
            return iter(self._f)

        def readlines(self):
            raise AssertionError("readlines() must not be called when streaming records!")

        def read(self, *args, **kwargs):
            return self._f.read(*args, **kwargs)

        def write(self, *args, **kwargs):
            return self._f.write(*args, **kwargs)

        def close(self):
            self._f.close()

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self._f.close()

    monkeypatch.setattr("builtins.open", NoReadlinesFile)

    changed = restrict_vcf_ann_to_alt(str(vcf_file))
    assert changed is True
    content = vcf_file.read_text()
    assert "ANN=" not in content


def test_restrict_vcf_ann_preserves_unchanged_file_and_indexes(tmp_path):
    """Verify unchanged inputs do not alter the file or its index timestamps."""
    import time

    from variantcentrifuge.vcf_eff_one_per_line import restrict_vcf_ann_to_alt

    vcf_file = tmp_path / "unchanged.vcf"
    raw_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tC\tA\t.\tPASS\tAC=1\n"
    )
    vcf_file.write_text(raw_content)
    idx_file = tmp_path / "unchanged.vcf.csi"
    idx_file.write_text("dummy_index")

    orig_vcf_mtime = vcf_file.stat().st_mtime
    orig_idx_mtime = idx_file.stat().st_mtime

    time.sleep(0.01)
    changed = restrict_vcf_ann_to_alt(str(vcf_file))
    assert changed is False
    assert vcf_file.exists()
    assert idx_file.exists()
    assert vcf_file.stat().st_mtime == orig_vcf_mtime
    assert idx_file.stat().st_mtime == orig_idx_mtime


def test_restrict_vcf_ann_preserves_original_file_on_write_failure(tmp_path, monkeypatch):
    """Verify failed writes leave the original file and indexes untouched."""
    from variantcentrifuge.vcf_eff_one_per_line import restrict_vcf_ann_to_alt

    vcf_file = tmp_path / "fail.vcf"
    raw_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tC\tA\t.\tPASS\tANN=G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||\n"
    )
    vcf_file.write_text(raw_content)
    idx_file = tmp_path / "fail.vcf.tbi"
    idx_file.write_text("dummy_tbi")

    def boom(*args, **kwargs):
        raise OSError("Disk full simulation")

    monkeypatch.setattr("os.replace", boom)

    import pytest

    with pytest.raises(OSError, match="Disk full simulation"):
        restrict_vcf_ann_to_alt(str(vcf_file))

    assert vcf_file.read_text() == raw_content
    assert idx_file.exists()


def test_restrict_vcf_ann_gz_streaming_without_readlines(tmp_path, monkeypatch):
    """Verify that restrict_vcf_ann_to_alt streams .vcf.gz without calling readlines()."""
    vcf_file = tmp_path / "test.vcf.gz"
    raw_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tC\tA\t.\tPASS\tANN=G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||\tGT\t0/1\n"
    )
    with gzip.open(vcf_file, "wt") as f:
        f.write(raw_content)

    orig_gzip_open = gzip.open

    class NoReadlinesGzipFile:
        def __init__(self, *args, **kwargs):
            self._f = orig_gzip_open(*args, **kwargs)

        def __iter__(self):
            return iter(self._f)

        def readlines(self):
            raise AssertionError("readlines() must not be called when streaming .vcf.gz records!")

        def read(self, *args, **kwargs):
            return self._f.read(*args, **kwargs)

        def write(self, *args, **kwargs):
            return self._f.write(*args, **kwargs)

        def close(self):
            self._f.close()

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self._f.close()

    monkeypatch.setattr("gzip.open", NoReadlinesGzipFile)

    changed = restrict_vcf_ann_to_alt(str(vcf_file))
    assert changed is True
    with orig_gzip_open(vcf_file, "rt") as f:
        content = f.read()
    assert "ANN=" not in content


def test_restrict_vcf_ann_preserves_file_permissions(tmp_path):
    """Verify that restrict_vcf_ann_to_alt retains strict file permissions (e.g. 0o600)."""
    import os
    import stat

    vcf_file = tmp_path / "secure.vcf"
    raw_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tC\tA\t.\tPASS\tANN=G|stop_gained|HIGH|KCNJ5|ENSG1|transcript|ENST1|protein_coding||||||||\tGT\t0/1\n"
    )
    vcf_file.write_text(raw_content)
    os.chmod(vcf_file, 0o600)

    initial_mode = stat.S_IMODE(os.stat(vcf_file).st_mode)
    assert initial_mode == 0o600

    changed = restrict_vcf_ann_to_alt(str(vcf_file))
    assert changed is True

    final_mode = stat.S_IMODE(os.stat(vcf_file).st_mode)
    assert final_mode == 0o600
