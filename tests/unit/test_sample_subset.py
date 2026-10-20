"""Sample subsetting at field extraction for association runs (#212, decision D1)."""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.aggregation import _aggregate_gene_burden_from_columns
from variantcentrifuge.association.genotype_matrix import build_genotype_matrix
from variantcentrifuge.checkpoint_helpers import hash_configuration
from variantcentrifuge.extractor import extract_fields_bcftools
from variantcentrifuge.sample_subset import (
    resolve_sample_subset,
    sample_set_sha256,
)

VCF_SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]


def _ctx(**cfg):
    config = {
        "perform_association": True,
        "no_stats": True,
        "case_samples": ["S2", "S5", "ABSENT"],  # ABSENT: in the sample file, not in the VCF
        "control_samples": ["S1", "S4"],
        **cfg,
    }
    return SimpleNamespace(
        config=config, vcf_samples=list(VCF_SAMPLES), args=SimpleNamespace(), phenotype_data=None
    )


def test_subset_applied_in_vcf_order_and_recorded():
    ctx = _ctx()
    assert resolve_sample_subset(ctx) == ["S1", "S2", "S4", "S5"]
    assert ctx.vcf_samples == ["S1", "S2", "S4", "S5"]
    assert ctx.config["sample_list"] == "S1,S2,S4,S5"
    assert ctx.config["sample_subset_sha256"] == sample_set_sha256(["S5", "S4", "S2", "S1"])
    # ABSENT never reaches bcftools (no --force-samples needed); cohort lists are untouched
    assert ctx.config["case_samples"] == ["S2", "S5", "ABSENT"]
    # cached: later callers (chunks, checkpoint restore) get the same answer
    assert resolve_sample_subset(ctx) == ["S1", "S2", "S4", "S5"]
    # the metadata TSV dumps the config: only the hash goes there, not the ID list
    assert not [
        k for k, v in ctx.config.items() if k.startswith("_") or v == ["S1", "S2", "S4", "S5"]
    ]


@pytest.mark.parametrize(
    "cfg",
    [
        {"perform_association": False},
        {"no_sample_subset": True},
        {"case_samples": [], "control_samples": []},
        {"control_samples": []},
        {"case_samples": "S2,S5"},  # unparsed CLI string
        {"control_samples": ["S1", "S3", "S4", "S6"]},  # union covers every VCF sample
    ],
)
def test_no_subset_without_strict_analysed_set(cfg, caplog):
    ctx = _ctx(**cfg)
    with caplog.at_level(logging.INFO):
        assert resolve_sample_subset(ctx) is None
    assert ctx.vcf_samples == VCF_SAMPLES
    assert "sample_subset_sha256" not in ctx.config
    assert "skipped" not in caplog.text


@pytest.mark.parametrize(
    ("cfg", "reason"),
    [
        ({"calculate_inheritance": True}, "inheritance"),
        ({"ped": "fam.ped"}, "inheritance"),
        ({"sample_sheet": "sheet.tsv"}, "inheritance"),
        ({"append_extra_sample_fields": True}, "--append-extra-sample-fields"),
        ({"no_replacement": True}, "--no-replacement"),
        ({"igv_enabled": True}, "--igv"),
        ({"no_stats": False}, "statistics"),
        ({"late_filtering": True}, "--late-filtering"),
        ({"pseudonymize": True}, "--pseudonymize"),
        ({"tumor_sample_index": 0}, "tumor/normal sample index"),
        ({"normal_sample_index": 1}, "tumor/normal sample index"),
        ({"final_filter": "IMPACT == 'HIGH'"}, "--final-filter"),
    ],
)
def test_sample_consuming_feature_keeps_all_samples_and_logs_once(cfg, reason, caplog):
    ctx = _ctx(**cfg)
    with caplog.at_level(logging.INFO, logger="variantcentrifuge.sample_subset"):
        assert resolve_sample_subset(ctx) is None
        assert resolve_sample_subset(ctx) is None
    assert ctx.vcf_samples == VCF_SAMPLES
    skipped = [r for r in caplog.records if "skipped" in r.getMessage()]
    assert len(skipped) == 1 and reason in skipped[0].getMessage()


def test_scoring_that_needs_inheritance_blocks_subset():
    with patch(
        "variantcentrifuge.scoring_requirements.check_scoring_requires_inheritance",
        return_value=True,
    ):
        assert resolve_sample_subset(_ctx()) is None


def test_remove_sample_substring_maps_back_to_header_names():
    ctx = _ctx(remove_sample_substring="_X", vcf_file="in.vcf.gz")
    with patch(
        "variantcentrifuge.helpers.get_vcf_samples", return_value=[f"{s}_X" for s in VCF_SAMPLES]
    ):
        assert resolve_sample_subset(ctx) == ["S1_X", "S2_X", "S4_X", "S5_X"]
    assert ctx.vcf_samples == ["S1", "S2", "S4", "S5"]


def test_quantitative_trait_uses_samples_with_a_value():
    ctx = _ctx(
        trait_type="quantitative",
        phenotype_file="pheno.tsv",
        phenotype_value_column="bmi",
        case_samples=[],
        control_samples=[],
    )
    ctx.phenotype_data = {"S3": {"21.5"}, "S6": {"30"}, "S2": {"NA"}}
    assert resolve_sample_subset(ctx) == ["S3", "S6"]
    assert resolve_sample_subset(_ctx(trait_type="quantitative")) is None  # no phenotype file


def test_checkpoint_hash_tracks_sample_set(tmp_path):
    cases = tmp_path / "cases.txt"
    cases.write_text("S2\nS5\n")
    cfg = {"vcf_file": "in.vcf.gz", "case_samples_file": str(cases)}
    h1 = hash_configuration(cfg)
    cases.write_text("S5\nS2\n\n")  # same set, other order
    assert hash_configuration(cfg) == h1
    cases.write_text("S2\nS6\n")
    assert hash_configuration(cfg) != h1
    assert hash_configuration({**cfg, "no_sample_subset": True}) != hash_configuration(cfg)


def test_chunk_fingerprint_includes_subset(tmp_path):
    from variantcentrifuge.stages.parallel_processing_stage import (
        ParallelCompleteProcessingStage,
    )

    bed = tmp_path / "c.bed"
    bed.write_text("chr1\t0\t10\n")
    stage = ParallelCompleteProcessingStage()
    ctx = SimpleNamespace(config={"vcf_file": "in.vcf.gz"}, data=None)
    fp_all = stage._chunk_fingerprint(ctx, bed, "CHROM GEN[*].GT")
    ctx.config["sample_subset_sha256"] = "abc"
    assert stage._chunk_fingerprint(ctx, bed, "CHROM GEN[*].GT") != fp_all


def test_score_artifact_provenance_records_subset():
    from variantcentrifuge.stages.association_stage_helpers import _emit_score_artifact

    with (
        patch("variantcentrifuge.association.meta.artifact.emit_stratum_artifact") as emit,
        patch(
            "variantcentrifuge.association.meta.provenance.collect_provenance",
            return_value={"regions_bed_sha256": "none"},
        ),
    ):
        _emit_score_artifact({"sample_subset_sha256": "f00"}, "/x/S1.npz", [])
        _emit_score_artifact({}, "/x/S1.npz", [])
    first, second = (c.kwargs["provenance"] for c in emit.call_args_list)
    assert first == {"regions_bed_sha256": "none", "sample_subset_sha256": "f00"}
    assert second["sample_subset_sha256"] == "none"


def test_extract_passes_samples_file_to_bcftools(tmp_path):
    seen = {}

    def fake_run(cmd, stdout, **_kw):
        seen["cmd"] = cmd
        seen["samples"] = Path(cmd[cmd.index("-S") + 1]).read_text()
        return subprocess.CompletedProcess(cmd, 0, "", "")

    with patch("variantcentrifuge.extractor.subprocess.run", side_effect=fake_run):
        extract_fields_bcftools(
            "in.vcf.gz",
            "CHROM GEN[*].GT",
            {"bcftools_samples": ["S1", "S4"]},
            str(tmp_path / "o.tsv"),
            vcf_samples=["S1", "S4"],
        )
    assert seen["cmd"][-1] == "in.vcf.gz" and seen["samples"] == "S1\nS4\n"
    assert not list(tmp_path.glob("*.samples"))


# --- real bcftools on a synthetic VCF -----------------------------------------------------

_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">
##INFO=<ID=GENE,Number=1,Type=String,Description="gene">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
chr1\t10\t.\tA\tG\t.\tPASS\tAC=3;GENE=G1\tGT\t0/1\t0/0\t1/1\t./.\t0/0\t0/0
chr1\t20\t.\tC\tT\t.\tPASS\tAC=2;GENE=G1\tGT\t0/0\t0/1\t0/0\t0/0\t0/1\t0/0
chr1\t30\t.\tG\tA\t.\tPASS\tAC=1;GENE=G2\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1
chr1\t40\t.\tT\tC\t.\tPASS\tAC=2;GENE=G2\tGT\t1/1\t./.\t0/0\t0/1\t0/0\t0/0
"""


@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not installed")
def test_subset_extraction_gives_identical_association_inputs(tmp_path):
    vcf = tmp_path / "syn.vcf"
    vcf.write_text(_VCF)
    fields = "CHROM POS AC GENE GEN[*].GT"
    full_tsv, sub_tsv = tmp_path / "full.tsv.gz", tmp_path / "sub.tsv.gz"
    extract_fields_bcftools(str(vcf), fields, {}, str(full_tsv), vcf_samples=VCF_SAMPLES)
    ctx = _ctx()
    subset = resolve_sample_subset(ctx)
    extract_fields_bcftools(
        str(vcf), fields, {"bcftools_samples": subset}, str(sub_tsv), vcf_samples=ctx.vcf_samples
    )
    with gzip.open(full_tsv, "rt") as f:
        full = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    with gzip.open(sub_tsv, "rt") as f:
        sub = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)

    assert list(sub["AC"]) == list(full["AC"]) == ["3", "2", "1", "2"]  # INFO as stored
    assert len(sub.columns) == len(full.columns) - 2
    for i, s in enumerate(ctx.vcf_samples):  # per-sample columns follow VCF order
        assert list(sub[f"GEN[{i}].GT"]) == list(full[f"GEN[{VCF_SAMPLES.index(s)}].GT"])

    cases, controls = set(ctx.config["case_samples"]), set(ctx.config["control_samples"])

    def gt_cols(df):
        return [c for c in df.columns if c.startswith("GEN[")]

    agg_full = _aggregate_gene_burden_from_columns(
        full, cases, controls, VCF_SAMPLES, gt_cols(full)
    )
    agg_sub = _aggregate_gene_burden_from_columns(
        sub, cases, controls, ctx.vcf_samples, gt_cols(sub)
    )
    assert agg_sub == agg_full

    assigned = [s for s in VCF_SAMPLES if s in cases | controls]
    pheno = np.array([1.0 if s in cases else 0.0 for s in assigned])
    for gene in ("G1", "G2"):
        by_sample_full = dict(zip(VCF_SAMPLES, gt_cols(full), strict=True))
        by_sample_sub = dict(zip(ctx.vcf_samples, gt_cols(sub), strict=True))
        m_full = build_genotype_matrix(
            full[full["GENE"] == gene],
            assigned,
            [by_sample_full[s] for s in assigned],
            phenotype_vector=pheno,
        )
        m_sub = build_genotype_matrix(
            sub[sub["GENE"] == gene],
            assigned,
            [by_sample_sub[s] for s in assigned],
            phenotype_vector=pheno,
        )
        for a, b in zip(m_full, m_sub, strict=True):
            np.testing.assert_array_equal(np.asarray(a), np.asarray(b))


def test_threaded_extraction_runs_after_case_control_assignment():
    from variantcentrifuge.pipeline_core.runner_planner import ordered_stages
    from variantcentrifuge.stages.parallel_processing_stage import (
        ParallelCompleteProcessingStage,
    )
    from variantcentrifuge.stages.phenotype_case_control_stages import (
        PhenotypeCaseControlAssignmentStage,
    )
    from variantcentrifuge.stages.post_processing_stages import GenotypeReplacementStage

    names = [
        s.name
        for s in ordered_stages(
            [
                GenotypeReplacementStage(),
                ParallelCompleteProcessingStage(),
                PhenotypeCaseControlAssignmentStage(),
            ]
        )
    ]
    assert (
        names.index("phenotype_case_control_assignment")
        < names.index("parallel_complete_processing")
        < names.index("genotype_replacement")
    )


@pytest.mark.parametrize("threads", [1, 8])
@pytest.mark.parametrize(
    "toggle",
    [{"no_stats": False}, {"pseudonymize": True}, {"igv_enabled": True}, {"ped": "f.ped"}],
)
def test_resume_refused_when_a_blocker_changes_the_extracted_samples(tmp_path, threads, toggle):
    """Repro shape: --threads N --enable-checkpoint --no-stats, then --resume with a blocker."""
    from variantcentrifuge.checkpoint_state import PipelineState

    cases, controls = tmp_path / "cases.txt", tmp_path / "controls.txt"
    cases.write_text("S2\nS5\n")
    controls.write_text("S1\nS4\n")
    run1 = {
        "vcf_file": "in.vcf.gz",
        "threads": threads,
        "perform_association": True,
        "no_stats": True,
        "case_samples_file": str(cases),
        "control_samples_file": str(controls),
    }
    state = PipelineState(str(tmp_path))
    state.initialize(run1, "0.23")
    resumed = PipelineState(str(tmp_path))
    assert resumed.load()
    assert resumed.can_resume(dict(run1), "0.23", validate_outputs=False)
    assert not resumed.can_resume({**run1, **toggle}, "0.23", validate_outputs=False)
