"""Regression tests for SnpSift filter fail-closed behavior."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from variantcentrifuge.config import load_config
from variantcentrifuge.field_profile import resolve_profile
from variantcentrifuge.filters import apply_snpsift_filter


pytestmark = pytest.mark.skipif(
    not all(shutil.which(tool) for tool in ("SnpSift", "bgzip", "bcftools")),
    reason="SnpSift, bgzip, and bcftools are required",
)


def _write_issue_104_vcf(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO">',
                '##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects">',
                '##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects">',
                '##INFO=<ID=dbNSFP_gnomAD4.1_joint_AF,Number=A,Type=Float,Description="AF">',
                '##INFO=<ID=dbNSFP_gnomAD4.1_joint_AC,Number=A,Type=Integer,Description="AC">',
                '##INFO=<ID=dbNSFP_clinvar_clnsig,Number=.,Type=String,Description="ClinVar">',
                '##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="ClinVar">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr4\t78912028\t.\tC\tT\t500\tPASS\tANN=T|stop_gained|HIGH|BMP2K|ENSG00000138756|transcript|ENST00000264889|protein_coding|1/10|c.3481C>T|p.Gln1161Ter|1/100|1/100|1/33||;dbNSFP_gnomAD4.1_joint_AF=0.00221454;dbNSFP_gnomAD4.1_joint_AC=3546;dbNSFP_clinvar_clnsig=Benign/Likely_benign;ClinVar_CLNSIG=Benign/Likely_benign",
                "chr10\t94545923\t.\tT\tC\t500\tPASS\tANN=C|start_lost|HIGH|HELLS|ENSG00000119969|transcript|ENST00000371332|protein_coding|1/20|c.2T>C|p.Met1?|2/100|2/100|1/33||;dbNSFP_gnomAD4.1_joint_AF=0.000997683;dbNSFP_gnomAD4.1_joint_AC=1552;dbNSFP_clinvar_clnsig=Benign/Likely_benign;ClinVar_CLNSIG=Benign/Likely_benign",
                "chr1\t40307478\t.\tG\tA\t500\tPASS\tANN=A|stop_gained|HIGH|COL9A2|ENSG00000049089|transcript|ENST00000372748|protein_coding|1/30|c.1G>A|p.Trp1Ter|1/100|1/100|1/33||;dbNSFP_gnomAD4.1_joint_AF=0.0141844;dbNSFP_gnomAD4.1_joint_AC=22894;dbNSFP_clinvar_clnsig=Benign/Likely_benign;ClinVar_CLNSIG=Benign/Likely_benign",
                "chr4\t78913000\t.\tG\tA\t500\tPASS\tANN=A|stop_gained|HIGH|KEEPER|ENSGKEEP|transcript|ENSTKEEP|protein_coding|1/10|c.1G>A|p.Trp1Ter|1/100|1/100|1/33||;dbNSFP_gnomAD4.1_joint_AF=0.00001;dbNSFP_gnomAD4.1_joint_AC=1;dbNSFP_clinvar_clnsig=Pathogenic;ClinVar_CLNSIG=Pathogenic",
                "",
            ]
        )
    )


def _dbnsfp5_filter(*presets: str) -> str:
    cfg = load_config()
    cfg["field_profile"] = "dbnsfp5"
    resolve_profile(cfg)
    return " & ".join(f"({cfg['presets'][name]})" for name in presets)


def test_issue_104_high_common_benign_record_is_filtered_out(tmp_path):
    input_vcf = tmp_path / "issue104.vcf"
    output_vcf = tmp_path / "filtered.vcf.gz"
    _write_issue_104_vcf(input_vcf)

    filter_expr = _dbnsfp5_filter("high_or_lof_or_nmd", "rare", "not_benign")

    apply_snpsift_filter(str(input_vcf), filter_expr, {"threads": 1}, str(output_vcf))

    result = subprocess.run(
        ["bcftools", "view", "-H", str(output_vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = [line for line in result.stdout.splitlines() if line]
    assert len(rows) == 1
    assert "KEEPER" in rows[0]
    assert "BMP2K" not in rows[0]
    assert "HELLS" not in rows[0]
    assert "COL9A2" not in rows[0]


def test_snpsift_exit_zero_parser_diagnostics_raise(tmp_path):
    input_vcf = tmp_path / "issue104.vcf"
    output_vcf = tmp_path / "filtered.vcf.gz"
    _write_issue_104_vcf(input_vcf)

    malformed_expr = (
        "((((exists LOF[*].PERC) & (LOF[*].PERC > 0.9)) | "
        "((exists NMD[*].PERC) & (NMD[*].PERC > 0.9)) | "
        "(ANN[ANY].IMPACT has 'HIGH')))) & "
        "(((dbNSFP_gnomAD4.1_joint_AF[0] < 0.0001) | "
        "(na dbNSFP_gnomAD4.1_joint_AC[0])))"
    )

    with pytest.raises(RuntimeError, match="SnpSift filter reported parser diagnostics"):
        apply_snpsift_filter(str(input_vcf), malformed_expr, {"threads": 1}, str(output_vcf))
