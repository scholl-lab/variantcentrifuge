"""Wildcard ANN[*]/EFF[*] extraction fields are rejected; ANN[0].GENE groups correctly."""

from __future__ import annotations

import shutil
import subprocess
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.extractor import (
    build_bcftools_format_string,
    extract_fields_bcftools,
    get_expected_extraction_columns,
    parse_ann_subfields,
    validate_extraction_fields,
)
from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage
from variantcentrifuge.stages.data_loading_stages import DataFrameLoadingStage
from variantcentrifuge.stages.setup_stages import ConfigurationLoadingStage

SAMPLES = [f"S{i}" for i in range(20)]


@pytest.mark.parametrize("field", ["ANN[*].GENE", "EFF[*].GENE", "ANN[*].EFFECT"])
def test_wildcard_annotation_fields_rejected(field):
    fields = ["CHROM", "POS", field, "GEN[*].GT"]
    with pytest.raises(ValueError, match=r"ANN\[0\]\.<FIELD>") as exc:
        build_bcftools_format_string(fields, SAMPLES)
    assert field in str(exc.value)
    with pytest.raises(ValueError, match="Unsupported extraction field"):
        get_expected_extraction_columns(" ".join(fields), SAMPLES)


def test_supported_fields_pass_validation():
    validate_extraction_fields(["CHROM", "ANN[0].GENE", "NMD[0].GENE", "GEN[*].GT", "dbNSFP_x"])


def test_configuration_loading_rejects_wildcard_fields_before_processing():
    context = PipelineContext(
        args=Namespace(config={"fields_to_extract": "CHROM POS ANN[*].GENE GEN[*].GT"}),
        config={},
        workspace=None,  # type: ignore[arg-type]
    )
    with pytest.raises(ValueError, match=r"ANN\[\*\]\.GENE"):
        ConfigurationLoadingStage()._process(context)


def _ann(gene: str, tx: str) -> str:
    return (
        f"C|missense_variant|MODERATE|{gene}|ID|transcript|{tx}|protein_coding|1/1|c.1A>C|p.X1Y"
        "|||||"
    )


def _association(tmp_path: Path, extracted: Path) -> pd.DataFrame:
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    context = PipelineContext(Namespace(), {"no_chunked_processing": True}, workspace)
    context.data = extracted
    context = DataFrameLoadingStage()._process(context)
    context.vcf_samples = SAMPLES
    context.config.update(
        perform_association=True,
        case_samples=SAMPLES[:10],
        control_samples=SAMPLES[10:],
        association_tests=["fisher"],
        output_dir=str(tmp_path),
        output_file_base="x",
        gzip_intermediates=False,
    )
    return AssociationAnalysisStage()._process(context).association_results


def _genotypes(v: int) -> list[str]:
    return ["0/1" if (i + v) % 4 == 0 else "0/0" for i in range(len(SAMPLES))]


def test_ann0_gene_groups_transcript_variants_into_one_gene(tmp_path):
    """Two TARGET variants on different transcripts form one gene with two variants."""
    raw = pd.DataFrame(
        {
            "CHROM": ["1", "1"],
            "POS": ["100", "200"],
            "REF": ["A", "A"],
            "ALT": ["C", "C"],
            "ANN": [_ann("TARGET", "TX_A"), _ann("TARGET", "TX_B")],
        }
    )
    df = parse_ann_subfields(raw, ["CHROM", "POS", "REF", "ALT", "ANN[0].GENE"])
    for j in range(len(SAMPLES)):
        df[f"GEN[{j}].GT"] = [_genotypes(0)[j], _genotypes(1)[j]]
    extracted = tmp_path / "extracted.tsv"
    df.to_csv(extracted, sep="\t", index=False)

    results = _association(tmp_path, extracted)
    assert results["gene"].tolist() == ["TARGET"]
    assert results["n_variants"].tolist() == [2]


@pytest.mark.skipif(not shutil.which("bcftools"), reason="requires bcftools on PATH")
def test_real_bcftools_ann0_extraction_aggregates_by_gene(tmp_path):
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=1000>",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
        "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | "
        "CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES),
    ]
    for v, (pos, tx) in enumerate([(100, "TX_A"), (200, "TX_B")]):
        lines.append(
            f"1\t{pos}\t.\tA\tC\t50\tPASS\tANN={_ann('TARGET', tx)}\tGT\t"
            + "\t".join(_genotypes(v))
        )
    vcf = tmp_path / "in.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    gz = tmp_path / "in.vcf.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(gz), str(vcf)], check=True)
    extracted = tmp_path / "extracted.tsv"
    extract_fields_bcftools(
        str(gz), "CHROM POS REF ALT ANN[0].GENE GEN[*].GT", {}, str(extracted), SAMPLES
    )

    results = _association(tmp_path, extracted)
    assert results["gene"].tolist() == ["TARGET"]
    assert results["n_variants"].tolist() == [2]
