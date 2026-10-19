"""Regression coverage for SnpEff split rows used in burden/association inputs."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.config import load_config
from variantcentrifuge.extractor import extract_fields_bcftools
from variantcentrifuge.filters import apply_snpsift_filter
from variantcentrifuge.gene_burden import perform_gene_burden_analysis
from variantcentrifuge.stages.processing_stages import ParallelCompleteProcessingStage
from variantcentrifuge.transcript_filter import filter_vcf_to_transcripts
from variantcentrifuge.vcf_eff_one_per_line import process_vcf_file as split_snpeff_annotations

pytestmark = pytest.mark.skipif(
    not all(shutil.which(tool) for tool in ("SnpSift", "bgzip", "bcftools", "tabix")),
    reason="SnpSift, bgzip, bcftools, and tabix are required",
)


FIELDS_TO_EXTRACT = "CHROM POS ANN[0].GENE ANN[0].IMPACT NMD[0].PERC GEN[*].GT"
VCF_SAMPLES = ["CASE", "CTRL"]


def _high_or_lof_or_nmd_filter() -> str:
    return load_config()["presets"]["high_or_lof_or_nmd"]


def _issue_106_record(*, pos: int, positive_for_gene_b: bool) -> str:
    if positive_for_gene_b:
        info = (
            "ANN="
            "C|downstream_gene_variant|MODIFIER|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||,"
            "C|missense_variant|MODERATE|GENE_B|ENSGB|transcript|NM_b|protein_coding||||||||;"
            "NMD=(GENE_B|ENSGB|1|0.95)"
        )
        return f"chr1\t{pos}\t.\tG\tC\t500\tPASS\t{info}\tGT\t0/1\t0/0"

    info = (
        "ANN="
        "T|splice_donor_variant|HIGH|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||,"
        "T|downstream_gene_variant|MODIFIER|GENE_B|ENSGB|transcript|NM_b|protein_coding||||||||;"
        "LOF=(GENE_A|ENSGA|1|0.95);"
        "NMD=(GENE_A|ENSGA|1|0.91)"
    )
    return f"chr1\t{pos}\t.\tA\tT\t500\tPASS\t{info}\tGT\t0/1\t0/0"


def _write_issue_106_vcf(
    path: Path,
    *,
    include_negative: bool = True,
    include_positive: bool = True,
) -> None:
    records = []
    if include_negative:
        records.append(_issue_106_record(pos=100, positive_for_gene_b=False))
    if include_positive:
        records.append(_issue_106_record(pos=200, positive_for_gene_b=True))

    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=1000>",
                (
                    '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional '
                    "annotations: Allele | Annotation | Annotation_Impact | Gene_Name | "
                    "Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | "
                    "HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | "
                    'AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO">'
                ),
                (
                    '##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of '
                    'function effects">'
                ),
                (
                    '##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense '
                    'mediated decay effects">'
                ),
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tCTRL",
                *records,
                "",
            ]
        ),
        encoding="utf-8",
    )


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _run_real_tool_flow(input_vcf: Path, tmp_path: Path) -> pd.DataFrame:
    tmp_path.mkdir(parents=True, exist_ok=True)
    split_vcf = tmp_path / "split.vcf.gz"
    filtered_vcf = tmp_path / "filtered.vcf.gz"
    transcript_vcf = tmp_path / "transcript_filtered.vcf.gz"
    extracted_tsv = tmp_path / "extracted.tsv.gz"

    split_snpeff_annotations(str(input_vcf), str(split_vcf))
    apply_snpsift_filter(
        str(split_vcf),
        _high_or_lof_or_nmd_filter(),
        {"threads": 1},
        str(filtered_vcf),
    )
    filter_vcf_to_transcripts(str(filtered_vcf), str(transcript_vcf), {"NM_b"})
    extract_fields_bcftools(
        variant_file=str(transcript_vcf),
        fields=FIELDS_TO_EXTRACT,
        cfg={"extract_fields_separator": ","},
        output_file=str(extracted_tsv),
        vcf_samples=VCF_SAMPLES,
    )

    return _read_tsv(extracted_tsv)


def _burden_input(df: pd.DataFrame) -> pd.DataFrame:
    renamed = df.rename(columns={"ANN[0].GENE": "GENE"})
    if "GENE" not in renamed.columns:
        renamed["GENE"] = pd.Series(dtype=str)
    return renamed


def test_split_filter_transcript_flow_excludes_cross_gene_issue_106_row(tmp_path: Path) -> None:
    input_vcf = tmp_path / "issue106.vcf"
    _write_issue_106_vcf(input_vcf)

    df = _run_real_tool_flow(input_vcf, tmp_path)

    assert df["POS"].tolist() == ["200"]
    assert df["ANN[0].GENE"].tolist() == ["GENE_B"]
    assert df["ANN[0].IMPACT"].tolist() == ["MODERATE"]
    assert df["NMD[0].PERC"].tolist() == ["0.95"]


def test_parallel_chunk_worker_before_filters_excludes_cross_gene_issue_106_row(
    tmp_path: Path,
) -> None:
    input_vcf = tmp_path / "issue106_negative.vcf"
    indexed_vcf = tmp_path / "issue106_negative.vcf.gz"
    chunk_bed = tmp_path / "chunk.bed"
    intermediate_dir = tmp_path / "worker"
    intermediate_dir.mkdir()

    _write_issue_106_vcf(input_vcf, include_positive=False)
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(indexed_vcf), str(input_vcf)], check=True)
    subprocess.run(["bcftools", "index", "-f", str(indexed_vcf)], check=True)
    chunk_bed.write_text("chr1\t0\t150\n", encoding="utf-8")

    chunk_tsv = ParallelCompleteProcessingStage()._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file=str(indexed_vcf),
        base_name="issue106",
        intermediate_dir=intermediate_dir,
        config={
            "threads_per_chunk": 1,
            "snpeff_splitting_mode": "before_filters",
            "filters": _high_or_lof_or_nmd_filter(),
            "fields_to_extract": FIELDS_TO_EXTRACT,
            "transcript_ids": {"NM_b"},
            "vcf_samples": VCF_SAMPLES,
            "gzip_intermediates": False,
        },
    )

    df = _read_tsv(chunk_tsv)
    assert df.empty
    assert "CHROM" in df.columns


def test_gene_burden_contract_uses_only_row_level_qualified_rows(tmp_path: Path) -> None:
    negative_vcf = tmp_path / "negative.vcf"
    positive_vcf = tmp_path / "positive.vcf"
    _write_issue_106_vcf(negative_vcf, include_positive=False)
    _write_issue_106_vcf(positive_vcf, include_negative=False)

    negative_df = _burden_input(_run_real_tool_flow(negative_vcf, tmp_path / "negative_flow"))
    positive_df = _burden_input(_run_real_tool_flow(positive_vcf, tmp_path / "positive_flow"))

    assert "GENE_B" not in set(negative_df["GENE"])

    burden_cfg = {
        "gene_burden_mode": "samples",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
    }
    negative_result = perform_gene_burden_analysis(
        negative_df,
        burden_cfg,
        case_samples={"CASE"},
        control_samples={"CTRL"},
        vcf_samples=VCF_SAMPLES,
    )
    positive_result = perform_gene_burden_analysis(
        positive_df,
        burden_cfg,
        case_samples={"CASE"},
        control_samples={"CTRL"},
        vcf_samples=VCF_SAMPLES,
    )

    assert negative_result.empty
    assert positive_result["GENE"].tolist() == ["GENE_B"]
    assert positive_result["proband_carrier_count"].tolist() == [1]
    assert positive_result["control_carrier_count"].tolist() == [0]


def test_transcript_filter_before_snpsift_excludes_unlisted_high_and_keeps_multi_isoform_rows(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.analysis_stages import VariantAnalysisStage

    input_vcf = tmp_path / "apa_isoforms.vcf"
    transcript_vcf = tmp_path / "transcript_filtered.vcf.gz"
    filtered_vcf = tmp_path / "filtered.vcf.gz"
    extracted_tsv = tmp_path / "extracted.tsv"

    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr3,length=200000000>",
                (
                    '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional '
                    "annotations: Allele | Annotation | Annotation_Impact | Gene_Name | "
                    "Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | "
                    "HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | "
                    'AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO">'
                ),
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL",
                # Record 1: HIGH on unlisted NM_999999.1, LOW on listed NM_000720.3
                (
                    "chr3\t100\t.\tA\tG\t500\tPASS\t"
                    "ANN=G|stop_gained|HIGH|CACNA1D|ENSG1|transcript|NM_999999.1|"
                    "protein_coding||||||||,"
                    "G|synonymous_variant|LOW|CACNA1D|ENSG1|transcript|NM_000720.3|"
                    "protein_coding||||||||\t"
                    "GT\t0/1\t0/0"
                ),
                # Record 2: Two listed CACNA1D isoforms (one version-skewed NM_000720.3 vs .4)
                (
                    "chr3\t200\t.\tG\tC\t500\tPASS\t"
                    "ANN=C|missense_variant|MODERATE|CACNA1D|ENSG1|transcript|NM_000720.3|"
                    "protein_coding|8/49|c.1207G>C|p.Gly403Arg|||||,"
                    "C|missense_variant|MODERATE|CACNA1D|ENSG1|transcript|NM_001128840.3|"
                    "protein_coding|8/48|c.1207G>C|p.Gly403Arg|||||\t"
                    "GT\t0/1\t0/0"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    with caplog.at_level("INFO", logger="variantcentrifuge.transcript_filter"):
        filter_vcf_to_transcripts(
            str(input_vcf),
            str(transcript_vcf),
            {"NM_000720.4", "NM_001128840.3", "NM_unmatched.1"},
        )

    assert "matched 2/3" in caplog.text
    assert "NM_unmatched.1" in caplog.text
    assert "NM_000720.4->NM_000720.3" in caplog.text

    apply_snpsift_filter(
        str(transcript_vcf),
        "(ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')",
        {"threads": 1},
        str(filtered_vcf),
    )
    extract_fields_bcftools(
        variant_file=str(filtered_vcf),
        fields="CHROM POS REF ALT ANN[0].GENE ANN[0].FEATUREID ANN[0].IMPACT GEN[*].GT",
        cfg={"extract_fields_separator": ","},
        output_file=str(extracted_tsv),
        vcf_samples=["TUMOR", "NORMAL"],
    )

    df = _read_tsv(extracted_tsv)
    assert len(df) == 2
    assert df["POS"].tolist() == ["200", "200"]
    assert df["ANN[0].FEATUREID"].tolist() == ["NM_000720.3", "NM_001128840.3"]

    df["GENE"] = df["ANN[0].GENE"]
    df["Extra_Info"] = ["iso1", "iso2"]
    workspace = Workspace(output_dir=tmp_path / "ws", base_name="test")
    workspace.intermediate_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        args=None,
        config={"no_stats": True, "gzip_intermediates": False},
        workspace=workspace,
        data=extracted_tsv,
    )
    ctx.vcf_samples = ["TUMOR", "NORMAL"]
    ctx.current_dataframe = df
    ctx.mark_complete("dataframe_loading")
    ctx.mark_complete("sample_config_loading")
    ctx = VariantAnalysisStage()._process(ctx)
    assert len(ctx.current_dataframe) == 2
    assert ctx.current_dataframe["ANN[0].FEATUREID"].tolist() == ["NM_000720.3", "NM_001128840.3"]
    assert ctx.current_dataframe["Extra_Info"].tolist() == ["iso1", "iso2"]


@pytest.mark.parametrize("reverse_order", [False, True])
def test_cross_gene_consequence_order_invariance(tmp_path: Path, reverse_order: bool) -> None:
    """Verify cross-gene consequence ordering produces identical qualified burden results."""
    ann_modifier = "T|intron_variant|MODIFIER|GENE_A|ENSGA|transcript|NM_a|protein_coding||||||||"
    ann_high = "T|stop_gained|HIGH|GENE_B|ENSGB|transcript|NM_b|protein_coding||||||||"
    ann_field = f"{ann_high},{ann_modifier}" if reverse_order else f"{ann_modifier},{ann_high}"

    vcf_content = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=1000>",
            (
                '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
                "Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
                "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | "
                "cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | "
                'Distance | ERRORS / WARNINGS / INFO">'
            ),
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tCTRL",
            f"chr1\t100\t.\tA\tT\t500\tPASS\tANN={ann_field}\tGT\t0/1\t0/0",
            "",
        ]
    )
    vcf_path = tmp_path / f"synthetic_{reverse_order}.vcf"
    vcf_path.write_text(vcf_content, encoding="utf-8")

    split_vcf = tmp_path / f"split_{reverse_order}.vcf.gz"
    filtered_vcf = tmp_path / f"filtered_{reverse_order}.vcf.gz"
    extracted_tsv = tmp_path / f"extracted_{reverse_order}.tsv"

    # Step 1: Split before filtering
    split_snpeff_annotations(str(vcf_path), str(split_vcf))

    # Step 2: Apply coding consequence filter
    coding_filter = "(ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')"
    apply_snpsift_filter(str(split_vcf), coding_filter, {"threads": 1}, str(filtered_vcf))

    # Step 3: Extract fields
    extract_fields_bcftools(
        variant_file=str(filtered_vcf),
        fields="CHROM POS ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
        cfg={"extract_fields_separator": ","},
        output_file=str(extracted_tsv),
        vcf_samples=["CASE", "CTRL"],
    )

    df = _read_tsv(extracted_tsv)
    assert len(df) == 1
    assert df["ANN[0].GENE"].tolist() == ["GENE_B"]
    assert df["ANN[0].IMPACT"].tolist() == ["HIGH"]

    # Step 4: Perform gene burden analysis
    burden_df = _burden_input(df)
    burden_cfg = {
        "gene_burden_mode": "samples",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
    }
    result = perform_gene_burden_analysis(
        burden_df,
        burden_cfg,
        case_samples={"CASE"},
        control_samples={"CTRL"},
        vcf_samples=["CASE", "CTRL"],
    )

    assert result["GENE"].tolist() == ["GENE_B"]
    assert result["proband_carrier_count"].tolist() == [1]
    assert result["control_carrier_count"].tolist() == [0]
