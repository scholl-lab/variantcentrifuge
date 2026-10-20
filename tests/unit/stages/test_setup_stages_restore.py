"""Setup-stage normalization and checkpoint-restore paths (resume must rebuild state)."""

from __future__ import annotations

import logging
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from variantcentrifuge.pipeline_core import PipelineContext
from variantcentrifuge.pipeline_core.runner_planner import restore_completed_stage
from variantcentrifuge.stages import setup_stages
from variantcentrifuge.stages.setup_stages import (
    AnnotationConfigLoadingStage,
    ConfigurationLoadingStage,
    PhenotypeLoadingStage,
    ScoringConfigLoadingStage,
    load_terms_from_file,
)

REPO = Path(__file__).resolve().parents[3]
log = logging.getLogger(__name__)


def _ctx(args: Namespace | None = None, **config: object) -> PipelineContext:
    return PipelineContext(args=args or Namespace(), config=dict(config), workspace=None)  # type: ignore[arg-type]


def test_load_terms_from_file(tmp_path):
    assert load_terms_from_file(None, log) == []
    terms = tmp_path / "t.txt"
    terms.write_text("HP:1\n\n HP:2 \n")
    assert load_terms_from_file(str(terms), log) == ["HP:1", "HP:2"]
    (tmp_path / "empty.txt").write_text("\n")
    with pytest.raises(ValueError, match="empty or invalid"):
        load_terms_from_file(str(tmp_path / "empty.txt"), log)
    with pytest.raises(FileNotFoundError):
        load_terms_from_file(str(tmp_path / "missing.txt"), log)


@pytest.mark.parametrize("config_path", [None, "custom.json"])
def test_configuration_merges_cli_args_over_loaded_config(config_path):
    args = Namespace(
        config=config_path,
        reference="hg38",
        threads=None,
        _internal="x",
        phenotype_file="p.tsv",
        phenotype_sample_column="s",
        phenotype_value_column="v",
    )
    with patch.object(setup_stages, "load_config", return_value={"threads": 4}) as load:
        result = ConfigurationLoadingStage()._process(_ctx(args, keep="me"))
    load.assert_called_once_with(*([config_path] if config_path else []))
    cfg = result.config
    assert cfg["reference"] == "hg38" and cfg["threads"] == 4 and cfg["keep"] == "me"
    assert "_internal" not in cfg and cfg["phenotype_value_column"] == "v"
    assert "extract" not in cfg  # only derived when fields_to_extract is configured


def test_configuration_parses_lists_and_term_files(tmp_path):
    case_file = tmp_path / "case.txt"
    case_file.write_text("HP:9\n")
    ctrl_file = tmp_path / "ctrl.txt"
    ctrl_file.write_text("HP:8\n")
    merged = {
        "fields_to_extract": "CHROM POS",
        "case_phenotypes": ["HP:1"],
        "control_phenotypes": ["HP:2"],
        "case_phenotypes_file": str(case_file),
        "control_phenotypes_file": str(ctrl_file),
        "case_samples": ["A"],
        "control_samples": ["B"],
    }
    cfg = ConfigurationLoadingStage()._process(_ctx(Namespace(config=merged))).config
    assert cfg["extract"] == ["CHROM", "POS"]
    assert cfg["case_phenotypes"] == ["HP:1", "HP:9"]
    assert cfg["control_phenotypes"] == ["HP:2", "HP:8"]
    assert cfg["case_samples"] == ["A"] and cfg["control_samples"] == ["B"]


def _pheno(tmp_path: Path, rows: str) -> dict:
    path = tmp_path / "pheno.tsv"
    path.write_text("sample\tvalue\n" + rows)
    return {
        "phenotype_file": str(path),
        "phenotype_sample_column": "sample",
        "phenotype_value_column": "value",
    }


@pytest.mark.parametrize("restore", [False, True])
def test_phenotype_loading_and_restore(tmp_path, restore):
    stage = PhenotypeLoadingStage()
    run = stage._handle_checkpoint_skip if restore else stage._process
    ctx = run(_ctx(**_pheno(tmp_path, "A\t1.5\nB\t2\n")))
    assert ctx.phenotype_data == {"A": {"1.5"}, "B": {"2"}}
    assert ctx.config["phenotypes"] == ctx.phenotype_data
    with pytest.raises(ValueError, match="No phenotype data loaded"):
        run(_ctx(**_pheno(tmp_path, "")))
    assert run(_ctx()).phenotype_data is None

    sheet = SimpleNamespace(phenotype_data={"A": {"HP:1"}})
    with patch("variantcentrifuge.sample_sheet.load_sample_sheet", return_value=sheet):
        ctx = run(_ctx(sample_sheet="sheet.tsv"))
    assert ctx.phenotype_data == {"A": {"HP:1"}}
    assert stage.checkpoint_resume_policy == "restore"


@pytest.mark.parametrize("restore", [False, True])
def test_scoring_config_loading_and_restore(restore):
    stage = ScoringConfigLoadingStage()
    run = stage._handle_checkpoint_skip if restore else stage._process
    ctx = run(_ctx(scoring_config_path=str(REPO / "scoring" / "test_simple_score")))
    assert ctx.scoring_config["formulas"]
    with pytest.raises(ValueError, match="scoring configuration"):
        run(_ctx(scoring_config_path=str(REPO / "does-not-exist")))
    assert run(_ctx()).scoring_config is None
    assert stage.checkpoint_resume_policy == "restore"
    assert stage.parallel_safe is True


@pytest.mark.parametrize("restore", [False, True])
def test_annotation_config_loading_and_restore(restore):
    stage = AnnotationConfigLoadingStage()
    run = stage._handle_checkpoint_skip if restore else stage._process
    ctx = run(
        _ctx(
            annotate_bed_files="a.bed",
            annotate_gene_lists="genes.txt",
            annotate_json_genes="g.json",
            json_gene_mapping='{"identifier":"symbol"}',
        )
    )
    assert ctx.annotation_configs == {
        "bed_files": ["a.bed"],
        "gene_lists": ["genes.txt"],
        "json_genes": ["g.json"],
        "json_mapping": '{"identifier":"symbol"}',
    }
    empty = run(_ctx())
    assert not empty.annotation_configs
    assert stage.checkpoint_resume_policy == "restore"


def test_configuration_restore_is_complete_after_resume():
    merged = {"case_samples": "A,B"}
    stage = ConfigurationLoadingStage()
    ctx = restore_completed_stage(stage, _ctx(Namespace(config=merged)))
    assert ctx.config["case_samples"] == ["A", "B"]
    assert ctx.is_complete("configuration_loading")
