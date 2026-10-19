"""Regression tests for Issue #161:
- Guard sys.exit() calls in helper functions in validators.py and cli_config.py
- Raise ValueError / FileNotFoundError from helpers so importing/calling them is composable
- Provide optional exit_on_error=True for CLI entry points needing SystemExit
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pytest

from variantcentrifuge.cli_config import (
    populate_cli_config,
    resolve_presets_and_tumor_normal,
)
from variantcentrifuge.validators import (
    validate_cli_args,
    validate_igv_files,
    validate_mandatory_parameters,
    validate_phenotype_file,
    validate_vcf_file,
)


@pytest.mark.unit
class TestValidatorsGuardSysExit:
    """Test that validators raise exceptions rather than forcing sys.exit."""

    def test_validate_mandatory_parameters_raises_value_error(self):
        """validate_mandatory_parameters should raise ValueError when missing inputs."""
        with pytest.raises(ValueError, match="reference database must be specified"):
            validate_mandatory_parameters(None, "filter1", "field1")

        with pytest.raises(ValueError, match="No filters provided"):
            validate_mandatory_parameters("GRCh38", None, "field1")

        with pytest.raises(ValueError, match="No fields to extract provided"):
            validate_mandatory_parameters("GRCh38", "filter1", None)

    def test_validate_mandatory_parameters_exit_on_error(self):
        """validate_mandatory_parameters should exit with SystemExit(1) when exit_on_error=True."""
        with pytest.raises(SystemExit) as exc:
            validate_mandatory_parameters(None, "filter1", "field1", exit_on_error=True)
        assert exc.value.code == 1

    def test_validate_vcf_file_missing_raises_file_not_found(self, tmp_path: Path):
        """validate_vcf_file should raise FileNotFoundError or ValueError for missing/empty VCF."""
        logger = logging.getLogger("test_logger")
        missing_vcf = str(tmp_path / "nonexistent.vcf")
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            validate_vcf_file(missing_vcf, logger)

        with pytest.raises(ValueError, match="VCF file not found"):
            validate_vcf_file(None, logger)

    def test_validate_vcf_file_empty_raises_value_error(self, tmp_path: Path):
        """validate_vcf_file should raise ValueError for 0-byte VCF."""
        logger = logging.getLogger("test_logger")
        empty_vcf = tmp_path / "empty.vcf"
        empty_vcf.touch()
        with pytest.raises(ValueError, match="is empty"):
            validate_vcf_file(str(empty_vcf), logger)

    def test_validate_vcf_file_exit_on_error(self, tmp_path: Path):
        """validate_vcf_file should exit with SystemExit(1) when exit_on_error=True."""
        logger = logging.getLogger("test_logger")
        missing_vcf = str(tmp_path / "nonexistent.vcf")
        with pytest.raises(SystemExit) as exc:
            validate_vcf_file(missing_vcf, logger, exit_on_error=True)
        assert exc.value.code == 1

    def test_validate_phenotype_file_missing_raises_file_not_found(self, tmp_path: Path):
        """validate_phenotype_file should raise FileNotFoundError when file is missing."""
        missing = str(tmp_path / "missing_pheno.tsv")
        with pytest.raises(FileNotFoundError, match="Phenotype file not found"):
            validate_phenotype_file(missing, "sample", "pheno")

    def test_validate_phenotype_file_empty_or_bad_header(self, tmp_path: Path):
        """validate_phenotype_file should raise ValueError when file is empty or missing columns."""
        empty = tmp_path / "empty.tsv"
        empty.touch()
        with pytest.raises(ValueError, match="is empty"):
            validate_phenotype_file(str(empty), "sample", "pheno")

        no_col = tmp_path / "bad.tsv"
        no_col.write_text("sample\twrong_col\ns1\t1\n")
        with pytest.raises(ValueError, match="Phenotype value column 'pheno' not found"):
            validate_phenotype_file(str(no_col), "sample", "pheno")

    def test_validate_igv_files_missing_raises(self, tmp_path: Path):
        """validate_igv_files should raise FileNotFoundError when fasta or fai is missing."""
        missing_fa = str(tmp_path / "genome.fa")
        with pytest.raises(FileNotFoundError, match="Local FASTA file not found"):
            validate_igv_files(missing_fa, None)

        fa = tmp_path / "genome.fa"
        fa.write_text(">chr1\nACGT\n")
        with pytest.raises(FileNotFoundError, match="FASTA index file not found"):
            validate_igv_files(str(fa), None)


@pytest.mark.unit
class TestCliConfigHelpersGuardSysExit:
    """Test that cli_config helpers raise ValueError instead of sys.exit."""

    def test_resolve_presets_missing_preset_raises_value_error(self):
        """resolve_presets_and_tumor_normal should raise ValueError if preset not in config."""
        args = argparse.Namespace(preset=["unknown_preset"], final_filter=None)
        cfg = {"presets": {"known": "AF < 0.01"}}
        with pytest.raises(ValueError, match="Preset 'unknown_preset' not found"):
            resolve_presets_and_tumor_normal(args, cfg, filters=None)

    def test_populate_cli_config_bad_json_raises(self):
        """populate_cli_config should raise ValueError when parser is None and JSON is invalid."""
        args = argparse.Namespace(
            variant_weight_params="not-json",
            diagnostics_output=None,
            regions_bed=None,
            pca=None,
            pca_components=10,
            pca_sites=None,
            coast_weights=None,
            association_workers=0,
            skat_method="SKAT",
            gene_prior_weights=None,
            fdr_method="bh",
            interaction_term=None,
            interaction_reference="dominant",
            permutation_test=False,
            permutations=1000,
            perm_adaptive=False,
            perm_alpha=0.05,
            perm_batch_size=1000,
            perm_seed=None,
            null_model_file=None,
            save_null_model=None,
            association_summary_tsv=None,
            sample_group=None,
            sample_group_column=None,
            sample_group_delimiter=None,
            sample_group_min_size=1,
            sample_group_reference=None,
            igv=False,
            annotate_json_genes=False,
            json_gene_mapping=None,
            annotate_gene_list=None,
            association_tests=None,
            perform_association=False,
            late_filtering=False,
            final_filter=None,
        )
        cfg: dict = {}
        with pytest.raises(ValueError, match="invalid JSON"):
            populate_cli_config(
                args, cfg, parser=None, reference="GRCh38", filters=None, fields=None
            )

    def test_populate_cli_config_bad_coast_weights_raises(self):
        """populate_cli_config should raise ValueError when coast_weights count != 3."""
        args = argparse.Namespace(
            variant_weight_params=None,
            diagnostics_output=None,
            regions_bed=None,
            pca=None,
            pca_components=10,
            pca_sites=None,
            coast_weights="0.1,0.2",  # Only 2 weights instead of 3
            association_workers=0,
            skat_method="SKAT",
            gene_prior_weights=None,
            fdr_method="bh",
            interaction_term=None,
            interaction_reference="dominant",
            permutation_test=False,
            permutations=1000,
            perm_adaptive=False,
            perm_alpha=0.05,
            perm_batch_size=1000,
            perm_seed=None,
            null_model_file=None,
            save_null_model=None,
            association_summary_tsv=None,
            sample_group=None,
            sample_group_column=None,
            sample_group_delimiter=None,
            sample_group_min_size=1,
            sample_group_reference=None,
            igv=False,
            annotate_json_genes=False,
            json_gene_mapping=None,
            annotate_gene_list=None,
            association_tests=None,
            perform_association=False,
            late_filtering=False,
            final_filter=None,
        )
        cfg: dict = {}
        with pytest.raises(ValueError, match="must have exactly 3 values"):
            populate_cli_config(
                args, cfg, parser=None, reference="GRCh38", filters=None, fields=None
            )


@pytest.mark.unit
class TestValidateCliArgsHelper:
    """Test the composable validate_cli_args helper function."""

    def test_validate_cli_args_missing_fields_file_raises(self, tmp_path: Path):
        """validate_cli_args should raise FileNotFoundError when fields_file is missing."""
        args = argparse.Namespace(
            fields_file=str(tmp_path / "missing_fields.txt"),
            transcript_file=None,
            phenotype_file=None,
            phenotype_sample_column=None,
            phenotype_value_column=None,
            resume=False,
            resume_from=None,
            interactive_resume=False,
            enable_checkpoint=False,
            pseudonymize=False,
            pseudonymize_table=None,
        )
        cfg: dict = {}
        with pytest.raises(FileNotFoundError, match="Fields file not found"):
            validate_cli_args(args, cfg)

    def test_validate_cli_args_incompatible_checkpoint_flags(self):
        """validate_cli_args should raise ValueError for incompatible checkpoint flags."""
        args = argparse.Namespace(
            fields_file=None,
            transcript_file=None,
            phenotype_file=None,
            phenotype_sample_column=None,
            phenotype_value_column=None,
            resume=True,
            resume_from=None,
            interactive_resume=False,
            enable_checkpoint=False,
            pseudonymize=False,
            pseudonymize_table=None,
        )
        cfg: dict = {}
        with pytest.raises(ValueError, match="--resume requires --enable-checkpoint"):
            validate_cli_args(args, cfg)


@pytest.mark.unit
class TestListingCommandsWithoutVcf:
    """Test that listing commands work without -v/--vcf-file or pipeline inputs."""

    def test_list_stages_without_vcf_succeeds(self, monkeypatch: pytest.MonkeyPatch):
        """--list-stages succeeds without -v and never invokes pipeline execution."""
        from unittest.mock import MagicMock, patch

        from variantcentrifuge.cli import main

        mock_run = MagicMock()
        monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run)

        with patch("sys.argv", ["variantcentrifuge", "--list-stages"]):
            with pytest.raises(SystemExit) as exc:
                main()
            assert exc.value.code == 0

        assert not mock_run.called

    def test_list_checkpoints_without_vcf_succeeds(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ):
        """--list-checkpoints succeeds without -v and never invokes pipeline execution."""
        from unittest.mock import MagicMock, patch

        from variantcentrifuge.cli import main

        mock_run = MagicMock()
        monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run)

        with patch(
            "sys.argv", ["variantcentrifuge", "--list-checkpoints", "--output-dir", str(tmp_path)]
        ):
            with pytest.raises(SystemExit) as exc:
                main()
            assert exc.value.code == 0

        assert not mock_run.called

    def test_ordinary_execution_requires_vcf(self, monkeypatch: pytest.MonkeyPatch):
        """Ordinary execution without -v/--vcf-file fails with error."""
        from unittest.mock import MagicMock, patch

        from variantcentrifuge.cli import main

        mock_run = MagicMock()
        monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run)

        with patch("sys.argv", ["variantcentrifuge", "--gene-name", "BRCA1"]):
            with pytest.raises(SystemExit) as exc:
                main()
            assert exc.value.code == 2

        assert not mock_run.called


@pytest.mark.unit
class TestRequiredToolsValidation:
    """Test that external prerequisites are validated upfront according to selected mode."""

    def test_missing_bcftools_raises_runtime_error(self, monkeypatch: pytest.MonkeyPatch):
        """Missing bcftools raises RuntimeError before pipeline execution."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="BRCA1", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None if t == "bcftools" else "/usr/bin/" + t)

        with pytest.raises(RuntimeError, match=r"Required external tool.*bcftools"):
            validate_required_tools(args)

    def test_missing_snpsift_raises_runtime_error(self, monkeypatch: pytest.MonkeyPatch):
        """Missing SnpSift raises RuntimeError before pipeline execution."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="BRCA1", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None if t == "SnpSift" else "/usr/bin/" + t)

        with pytest.raises(RuntimeError, match=r"Required external tool.*SnpSift"):
            validate_required_tools(args)

    def test_snpeff_bypassed_with_gene_all_and_regions_bed(self, monkeypatch: pytest.MonkeyPatch):
        """snpEff is not required when gene_name is 'all' and regions_bed is provided."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="all", regions_bed="regions.bed")
        # snpEff missing, but bcftools and SnpSift present
        monkeypatch.setattr("shutil.which", lambda t: None if t == "snpEff" else "/usr/bin/" + t)

        # Should NOT raise because snpEff is bypassed
        validate_required_tools(args)

    def test_snpeff_required_when_specific_gene_named(self, monkeypatch: pytest.MonkeyPatch):
        """snpEff is required when a specific gene name is requested."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="BRCA1", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None if t == "snpEff" else "/usr/bin/" + t)

        with pytest.raises(RuntimeError, match=r"Required external tool.*snpEff"):
            validate_required_tools(args)


@pytest.mark.unit
class TestBgzfAndIndexValidation:
    """Test upfront validation for BGZF compression and tabix/csi indexes."""

    def test_plain_vcf_fails_validation(self, tmp_path: Path):
        """Plain uncompressed VCF fails with ValueError."""
        plain_vcf = tmp_path / "plain.vcf"
        plain_vcf.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

        with pytest.raises(ValueError, match="uncompressed plain text"):
            validate_vcf_file(str(plain_vcf))

    def test_standard_gzip_fails_validation(self, tmp_path: Path):
        """Standard gzip VCF without BGZF extra field fails with ValueError."""
        import gzip

        gz_vcf = tmp_path / "standard.vcf.gz"
        with gzip.open(gz_vcf, "wt") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        with pytest.raises(ValueError, match="standard gzip compressed, not BGZF"):
            validate_vcf_file(str(gz_vcf))

    def test_missing_index_fails_validation(self, tmp_path: Path):
        """Valid BGZF VCF without index fails with FileNotFoundError."""
        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())

        with pytest.raises(FileNotFoundError, match="VCF index file not found"):
            validate_vcf_file(str(dest_vcf))

    def test_empty_index_fails_validation(self, tmp_path: Path):
        """Valid BGZF VCF with empty (0-byte) index fails with ValueError."""
        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())
        dest_tbi = tmp_path / "test.vcf.gz.tbi"
        dest_tbi.touch()

        with pytest.raises(ValueError, match="is empty"):
            validate_vcf_file(str(dest_vcf))

    def test_valid_bgzf_with_tbi_passes(self, tmp_path: Path):
        """Valid BGZF VCF with valid .tbi index passes validation."""
        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        fixture_tbi = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz.tbi")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_tbi = tmp_path / "test.vcf.gz.tbi"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())
        dest_tbi.write_bytes(fixture_tbi.read_bytes())

        # Should pass without error
        validate_vcf_file(str(dest_vcf))

    def test_valid_bgzf_with_csi_passes(self, tmp_path: Path):
        """Valid BGZF VCF with valid .csi index passes validation."""
        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        fixture_tbi = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz.tbi")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_csi = tmp_path / "test.vcf.gz.csi"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())
        dest_csi.write_bytes(fixture_tbi.read_bytes())

        # Should pass without error
        validate_vcf_file(str(dest_vcf))


@pytest.mark.unit
class TestFinalFilterValidation:
    """Test final filter expression syntax and missing-column error handling."""

    def test_final_filter_missing_column_fails(self):
        """Filtering with a missing column name raises ValueError and aborts."""
        import pandas as pd

        from variantcentrifuge.filters import filter_dataframe_with_query

        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["100"], "score": ["0.5"]})
        with pytest.raises(ValueError, match="name 'NONEXISTENT_COL' is not defined"):
            filter_dataframe_with_query(df, "NONEXISTENT_COL > 10")

    def test_final_filter_invalid_syntax_fails(self):
        """Filtering with invalid syntax raises ValueError and aborts."""
        import pandas as pd

        from variantcentrifuge.filters import filter_dataframe_with_query

        df = pd.DataFrame({"CHROM": ["chr1"], "score": ["0.5"]})
        with pytest.raises(ValueError, match="Invalid final filter expression"):
            filter_dataframe_with_query(df, "score === 0.5")

    def test_valid_final_filter_retains_expected_rows(self):
        """Valid final filter expression retains the expected matching rows."""
        import pandas as pd

        from variantcentrifuge.filters import filter_dataframe_with_query

        df = pd.DataFrame(
            {
                "GENE": ["BRCA1", "TP53", "EGFR"],
                "score": [0.8, 0.4, 0.9],
            }
        )
        res = filter_dataframe_with_query(df, "score > 0.5")
        assert len(res) == 2
        assert list(res["GENE"]) == ["BRCA1", "EGFR"]
