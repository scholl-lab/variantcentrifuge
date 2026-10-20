"""cfg["temp_dir"] must reach every bulk temp writer through the real stage code (#213)."""

import contextlib
from pathlib import Path
from unittest.mock import patch

from tests.unit.stages.test_processing_stages_critical import (  # noqa: F401 (fixtures)
    base_context,
    mock_workspace,
)
from variantcentrifuge.stages.processing_stages import (
    FieldExtractionStage,
    ParallelCompleteProcessingStage,
    SnpSiftFilterStage,
    VariantExtractionStage,
)

TMP = "/scratch/vc_temp_x"


def test_snpsift_filter_stage_passes_temp_dir(base_context, tmp_path):  # noqa: F811
    vcf = tmp_path / "in.vcf.gz"
    vcf.touch()
    base_context.config = {"filter": "(AF > 0.1)", "temp_dir": TMP}
    base_context.data = vcf
    with patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter") as m:
        SnpSiftFilterStage()._process(base_context)
    assert m.call_args.args[2]["temp_dir"] == TMP


def test_field_extraction_stage_passes_temp_dir(base_context, tmp_path):  # noqa: F811
    vcf = tmp_path / "f.vcf.gz"
    vcf.touch()
    base_context.config = {"fields_to_extract": "CHROM POS", "temp_dir": TMP}
    base_context.filtered_vcf = vcf
    with patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools") as m:
        FieldExtractionStage()._process(base_context)
    assert m.call_args.kwargs["cfg"]["temp_dir"] == TMP


def test_variant_extraction_stage_passes_temp_dir(base_context, tmp_path):  # noqa: F811
    vcf = tmp_path / "in.vcf.gz"
    vcf.touch()
    base_context.config = {"vcf_file": str(vcf), "temp_dir": TMP}
    base_context.gene_bed_file = tmp_path / "g.bed"
    with (
        patch("variantcentrifuge.stages.processing_stages.extract_variants") as m,
        contextlib.suppress(Exception),  # only the call arguments matter here
    ):
        VariantExtractionStage()._process(base_context)
    assert m.call_args.kwargs["cfg"]["temp_dir"] == TMP


def test_parallel_chunk_passes_temp_dir(tmp_path: Path):
    bed = tmp_path / "c.bed"
    bed.write_text("1\t0\t100\n")
    cfg = {
        "threads_per_chunk": 1,
        "fields_to_extract": "CHROM POS",
        "filters": "(AF > 0.1)",
        "temp_dir": TMP,
        "vcf_samples": ["S1"],
    }
    p = "variantcentrifuge.stages.processing_stages."
    with (
        patch(p + "extract_variants") as ev,
        patch(p + "apply_snpsift_filter") as af,
        patch(p + "extract_fields_bcftools") as ef,
    ):
        ParallelCompleteProcessingStage()._process_single_chunk(
            0, bed, "in.vcf.gz", "s", tmp_path, cfg
        )
    assert ev.call_args.kwargs["cfg"]["temp_dir"] == TMP
    assert af.call_args.args[2]["temp_dir"] == TMP
    assert ef.call_args.kwargs["cfg"]["temp_dir"] == TMP
