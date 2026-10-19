"""Regression tests for Issue #161:
- Guard sys.exit() calls in helper functions in validators.py and cli_config.py
- Raise ValueError / FileNotFoundError from helpers so importing/calling them is composable
- Provide optional exit_on_error=True for CLI entry points needing SystemExit
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from unittest.mock import Mock, patch

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
        """snpEff, bedtools, and sortBed are not required when gene_name is 'all'
        and regions_bed is provided.
        """
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="all", regions_bed="regions.bed")
        # snpEff, bedtools, sortBed missing, but bcftools and SnpSift present
        monkeypatch.setattr(
            "shutil.which",
            lambda t: None if t in ("snpEff", "bedtools", "sortBed") else "/usr/bin/" + t,
        )

        # Should NOT raise because regions_bed is used directly without get_gene_bed
        validate_required_tools(args)

    def test_snpeff_bedtools_sortbed_required_with_gene_all_without_regions_bed(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """When gene_name is 'all' without regions_bed, get_gene_bed runs and
        requires snpEff, bedtools, sortBed.
        """
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="all", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None if t == "snpEff" else "/usr/bin/" + t)
        with pytest.raises(RuntimeError, match=r"Required external tool.*snpEff"):
            validate_required_tools(args)

        monkeypatch.setattr("shutil.which", lambda t: None if t == "bedtools" else "/usr/bin/" + t)
        with pytest.raises(RuntimeError, match=r"Required external tool.*bedtools"):
            validate_required_tools(args)

        monkeypatch.setattr("shutil.which", lambda t: None if t == "sortBed" else "/usr/bin/" + t)
        with pytest.raises(RuntimeError, match=r"Required external tool.*sortBed"):
            validate_required_tools(args)

    def test_snpeff_required_when_specific_gene_named(self, monkeypatch: pytest.MonkeyPatch):
        """snpEff is required when a specific gene name is requested."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="BRCA1", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None if t == "snpEff" else "/usr/bin/" + t)

        with pytest.raises(RuntimeError, match=r"Required external tool.*snpEff"):
            validate_required_tools(args)

    def test_bedtools_and_sortbed_required_for_interval_extraction(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """bedtools and sortBed are required when specific gene extraction runs."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="BRCA1", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None if t == "bedtools" else "/usr/bin/" + t)
        with pytest.raises(RuntimeError, match=r"Required external tool.*bedtools"):
            validate_required_tools(args)

        monkeypatch.setattr("shutil.which", lambda t: None if t == "sortBed" else "/usr/bin/" + t)
        with pytest.raises(RuntimeError, match=r"Required external tool.*sortBed"):
            validate_required_tools(args)

    def test_akt_required_when_pca_akt(self, monkeypatch: pytest.MonkeyPatch):
        """akt is required when PCA tool is akt."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="all", regions_bed="regions.bed", pca="akt")
        monkeypatch.setattr("shutil.which", lambda t: None if t == "akt" else "/usr/bin/" + t)
        with pytest.raises(RuntimeError, match=r"Required external tool.*akt"):
            validate_required_tools(args)

    def test_snpsift_not_required_when_late_filtering_without_early_filters(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """SnpSift is not required when late filtering is enabled and no early filters run."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(
            gene_name="all",
            regions_bed="regions.bed",
            late_filtering=True,
            filters=None,
        )
        monkeypatch.setattr("shutil.which", lambda t: None if t == "SnpSift" else "/usr/bin/" + t)
        # Should not raise because SnpSift is not required
        validate_required_tools(args)

    def test_snpsift_not_required_when_late_filtering_with_filters(
        self, monkeypatch: pytest.MonkeyPatch
    ):
        """SnpSift is not required when late filtering is enabled, even if filters are specified."""
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(
            gene_name="all",
            regions_bed="regions.bed",
            late_filtering=True,
            filters="QUAL > 10",
        )
        monkeypatch.setattr("shutil.which", lambda t: None if t == "SnpSift" else "/usr/bin/" + t)
        # Should not raise because SnpSift is not required when late_filtering=True
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
        import gzip

        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_csi = tmp_path / "test.vcf.gz.csi"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())
        dest_csi.write_bytes(gzip.compress(b"CSI\x01" + b"\x00" * 30))

        # Should pass without error
        validate_vcf_file(str(dest_vcf))

    def test_corrupted_csi_fails_validation(self, tmp_path: Path):
        """Valid BGZF VCF with invalid magic bytes in .csi fails with ValueError."""
        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_csi = tmp_path / "test.vcf.gz.csi"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())
        dest_csi.write_bytes(b"NOT_A_VALID_GZIP_OR_CSI")

        with pytest.raises(ValueError, match="corrupted or not a valid gzip/BGZF index"):
            validate_vcf_file(str(dest_vcf))

    def test_corrupted_tbi_fails_validation(self, tmp_path: Path):
        """Valid BGZF VCF with invalid magic bytes in .tbi fails with ValueError."""
        import gzip

        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        dest_vcf = tmp_path / "test.vcf.gz"
        dest_tbi = tmp_path / "test.vcf.gz.tbi"
        dest_vcf.write_bytes(fixture_vcf.read_bytes())
        dest_tbi.write_bytes(gzip.compress(b"BAD\x01" + b"\x00" * 30))

        with pytest.raises(ValueError, match="invalid format or magic header"):
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


@pytest.mark.unit
class TestValidatorsComprehensive:
    """Comprehensive branch coverage for variantcentrifuge.validators."""

    def test_is_bgzf_edge_cases(self):
        from variantcentrifuge.validators import _is_bgzf

        assert _is_bgzf(b"short") is False
        assert _is_bgzf(b"\x1f\x8b\x00" + b"\x00" * 20) is False  # Wrong magic
        assert _is_bgzf(b"\x1f\x8b\x08\x00" + b"\x00" * 20) is False  # No FEXTRA
        # FEXTRA set but xlen extends beyond buffer
        assert _is_bgzf(b"\x1f\x8b\x08\x04" + b"\x00" * 6 + b"\xff\xff" + b"\x00" * 10) is False
        # FEXTRA with other subfield
        header = b"\x1f\x8b\x08\x04" + b"\x00" * 6 + b"\x06\x00" + b"XY\x02\x00\x00\x00"
        assert _is_bgzf(header) is False

    def test_validate_vcf_file_missing_and_exit(self, tmp_path: Path):
        from variantcentrifuge.validators import validate_vcf_file

        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            validate_vcf_file(str(tmp_path / "nonexistent.vcf.gz"))

        with pytest.raises(SystemExit):
            validate_vcf_file(str(tmp_path / "nonexistent.vcf.gz"), exit_on_error=True)

    def test_validate_required_tools_exit_on_error(self, monkeypatch: pytest.MonkeyPatch):
        from variantcentrifuge.validators import validate_required_tools

        args = argparse.Namespace(gene_name="all", regions_bed=None)
        monkeypatch.setattr("shutil.which", lambda t: None)
        with pytest.raises(SystemExit):
            validate_required_tools(args, exit_on_error=True)

    def test_validate_phenotype_file_branches(self, tmp_path: Path):
        from variantcentrifuge.validators import validate_phenotype_file

        # None is no-op
        validate_phenotype_file(None, "sample", "val")

        # Missing file
        with pytest.raises(FileNotFoundError):
            validate_phenotype_file(str(tmp_path / "missing.tsv"), "sample", "val")
        with pytest.raises(SystemExit):
            validate_phenotype_file(
                str(tmp_path / "missing.tsv"), "sample", "val", exit_on_error=True
            )

        # Empty file
        empty_file = tmp_path / "empty.tsv"
        empty_file.touch()
        with pytest.raises(ValueError, match="is empty"):
            validate_phenotype_file(str(empty_file), "sample", "val")

        # No header
        no_hdr = tmp_path / "no_hdr.tsv"
        no_hdr.write_text("\n\n")
        with pytest.raises(ValueError, match="has no header"):
            validate_phenotype_file(str(no_hdr), "sample", "val")

        # Missing sample col
        bad_sample = tmp_path / "bad_sample.tsv"
        bad_sample.write_text("other\tval\nS1\t1\n")
        with pytest.raises(ValueError, match="sample column 'sample' not found"):
            validate_phenotype_file(str(bad_sample), "sample", "val")

        # Missing val col
        bad_val = tmp_path / "bad_val.tsv"
        bad_val.write_text("sample\tother\nS1\t1\n")
        with pytest.raises(ValueError, match="value column 'val' not found"):
            validate_phenotype_file(str(bad_val), "sample", "val")

        # Header only no data
        hdr_only = tmp_path / "hdr_only.tsv"
        hdr_only.write_text("sample\tval\n\n")
        with pytest.raises(ValueError, match="contains only a header"):
            validate_phenotype_file(str(hdr_only), "sample", "val")

        # Valid phenotype file
        good_pheno = tmp_path / "good.tsv"
        good_pheno.write_text("sample\tval\nS1\t1\n")
        validate_phenotype_file(str(good_pheno), "sample", "val")

    def test_validate_mandatory_parameters_branches(self):
        from variantcentrifuge.validators import validate_mandatory_parameters

        with pytest.raises(ValueError, match="reference database must be specified"):
            validate_mandatory_parameters(None, "filter", "field")
        with pytest.raises(SystemExit):
            validate_mandatory_parameters(None, "filter", "field", exit_on_error=True)

        with pytest.raises(ValueError, match="No filters provided"):
            validate_mandatory_parameters("ref", None, "field")
        with pytest.raises(SystemExit):
            validate_mandatory_parameters("ref", None, "field", exit_on_error=True)

        with pytest.raises(ValueError, match="No fields to extract provided"):
            validate_mandatory_parameters("ref", "filter", None)
        with pytest.raises(SystemExit):
            validate_mandatory_parameters("ref", "filter", None, exit_on_error=True)

        # Valid
        validate_mandatory_parameters("ref", "filter", "field")

    def test_validate_igv_files_branches(self, tmp_path: Path):
        from variantcentrifuge.validators import validate_igv_files

        # Missing fasta
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            validate_igv_files(str(tmp_path / "missing.fa"), None)

        # Empty fasta
        empty_fa = tmp_path / "empty.fa"
        empty_fa.touch()
        with pytest.raises(ValueError, match="FASTA file is empty"):
            validate_igv_files(str(empty_fa), None)

        # Missing fasta index
        fa = tmp_path / "ref.fa"
        fa.write_text(">chr1\nACGT\n")
        with pytest.raises(FileNotFoundError, match="FASTA index file not found"):
            validate_igv_files(str(fa), None)

        # Empty index
        fai = tmp_path / "ref.fa.fai"
        fai.touch()
        with pytest.raises(ValueError, match="FASTA index file is empty"):
            validate_igv_files(str(fa), None)

        # Valid fasta + fai
        fai.write_text("chr1\t4\t6\t4\t5\n")
        validate_igv_files(str(fa), None)

        # Missing ideogram
        with pytest.raises(FileNotFoundError, match="Ideogram file not found"):
            validate_igv_files(str(fa), str(tmp_path / "missing.ideo"))

        # Empty ideogram
        empty_ideo = tmp_path / "empty.ideo"
        empty_ideo.touch()
        with pytest.raises(ValueError, match="Ideogram file is empty"):
            validate_igv_files(str(fa), str(empty_ideo))

        # Valid ideogram
        good_ideo = tmp_path / "good.ideo"
        good_ideo.write_text("chr1\t0\t100\tp36.33\tgneg\n")
        validate_igv_files(str(fa), str(good_ideo))

    def test_validate_cli_args_branches(self, tmp_path: Path):
        from variantcentrifuge.validators import validate_cli_args

        # Missing fields_file
        args = argparse.Namespace(fields_file=str(tmp_path / "missing_fields.txt"), fields=None)
        with pytest.raises(FileNotFoundError, match="Fields file not found"):
            validate_cli_args(args, {})

        # Missing transcript_file
        args = argparse.Namespace(transcript_file=str(tmp_path / "missing_tx.txt"))
        with pytest.raises(FileNotFoundError, match="Transcript file not found"):
            validate_cli_args(args, {})

        # Phenotype file without sample/val cols
        args = argparse.Namespace(
            phenotype_file="pheno.tsv", phenotype_sample_column=None, phenotype_value_column=None
        )
        with pytest.raises(ValueError, match="both --phenotype-sample-column"):
            validate_cli_args(args, {})

        # Resume without enable_checkpoint
        args = argparse.Namespace(resume=True, enable_checkpoint=False)
        with pytest.raises(ValueError, match="--resume requires --enable-checkpoint"):
            validate_cli_args(args, {})

        args = argparse.Namespace(resume_from="s1", enable_checkpoint=False)
        with pytest.raises(ValueError, match="--resume-from requires --enable-checkpoint"):
            validate_cli_args(args, {})

        args = argparse.Namespace(interactive_resume=True, enable_checkpoint=False)
        with pytest.raises(ValueError, match="--interactive-resume requires --enable-checkpoint"):
            validate_cli_args(args, {})

        # Multiple resume options
        args = argparse.Namespace(
            resume=True,
            resume_from="s1",
            interactive_resume=False,
            enable_checkpoint=True,
        )
        with pytest.raises(ValueError, match="Cannot use multiple resume options"):
            validate_cli_args(args, {})

        # Pseudonymize without table
        args = argparse.Namespace(pseudonymize=True, pseudonymize_table=None)
        with pytest.raises(ValueError, match="--pseudonymize-table is required"):
            validate_cli_args(args, {})

        # Pseudonymize with cfg pseudonymize_table
        args = argparse.Namespace(pseudonymize=True, pseudonymize_table=None)
        validate_cli_args(args, {"pseudonymize_table": "table.tsv"})

        # Parser error paths
        mock_parser = Mock()
        mock_parser.error = Mock(side_effect=SystemExit(2))

        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(fields_file="bad.txt", fields=None), {}, parser=mock_parser
            )
        with pytest.raises(SystemExit):
            validate_cli_args(argparse.Namespace(transcript_file="bad.txt"), {}, parser=mock_parser)
        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(
                    phenotype_file="p.tsv",
                    phenotype_sample_column=None,
                    phenotype_value_column=None,
                ),
                {},
                parser=mock_parser,
            )
        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(resume=True, enable_checkpoint=False), {}, parser=mock_parser
            )
        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(resume_from="s", enable_checkpoint=False), {}, parser=mock_parser
            )
        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(interactive_resume=True, enable_checkpoint=False),
                {},
                parser=mock_parser,
            )
        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(
                    resume=True,
                    resume_from="s",
                    interactive_resume=False,
                    enable_checkpoint=True,
                ),
                {},
                parser=mock_parser,
            )
        with pytest.raises(SystemExit):
            validate_cli_args(
                argparse.Namespace(pseudonymize=True, pseudonymize_table=None),
                {},
                parser=mock_parser,
            )

    def test_validators_with_mock_logger(self, tmp_path: Path):
        from variantcentrifuge.validators import (
            validate_phenotype_file,
            validate_required_tools,
            validate_vcf_file,
        )

        mock_logger = Mock()

        # Missing vcf
        with pytest.raises(FileNotFoundError):
            validate_vcf_file("missing.vcf.gz", logger=mock_logger)
        assert mock_logger.error.called

        # Empty index
        vcf = tmp_path / "test.vcf.gz"
        fixture_vcf = Path("tests/fixtures/test_variants_scoring.GRCh37.annotated.vcf.gz")
        vcf.write_bytes(fixture_vcf.read_bytes())
        tbi = tmp_path / "test.vcf.gz.tbi"
        tbi.touch()
        mock_logger.reset_mock()
        with pytest.raises(ValueError):
            validate_vcf_file(str(vcf), logger=mock_logger)
        assert mock_logger.error.called

        # Phenotype errors with logger
        mock_logger.reset_mock()
        with pytest.raises(FileNotFoundError):
            validate_phenotype_file("missing.tsv", "s", "v", logger=mock_logger)
        assert mock_logger.error.called

        empty_tsv = tmp_path / "empty_l.tsv"
        empty_tsv.touch()
        mock_logger.reset_mock()
        with pytest.raises(ValueError):
            validate_phenotype_file(str(empty_tsv), "s", "v", logger=mock_logger)
        assert mock_logger.error.called

        # Missing tools with logger
        args = argparse.Namespace(gene_name="BRCA1", regions_bed=None)
        mock_logger.reset_mock()
        with (
            patch("shutil.which", return_value=None),
            pytest.raises(RuntimeError),
        ):
            validate_required_tools(args, logger=mock_logger)
        assert mock_logger.error.called

    def test_validate_vcf_file_bcf_and_exit_paths(self, tmp_path: Path):
        import gzip

        from variantcentrifuge.validators import validate_igv_files, validate_vcf_file

        # BCF file with .csi adjacent index
        bcf = tmp_path / "test.bcf"
        # BCF header starts with BCF\x02
        bcf.write_bytes(b"BCF\x02\x01\x00\x00\x00" + b"\x00" * 30)
        csi = tmp_path / "test.csi"
        csi.write_bytes(gzip.compress(b"CSI\x01" + b"\x00" * 30))

        validate_vcf_file(str(bcf))

        # validate_igv_files exit_on_error paths
        with pytest.raises(SystemExit):
            validate_igv_files("missing.fa", None, exit_on_error=True)

        empty_fa = tmp_path / "igv_empty.fa"
        empty_fa.touch()
        with pytest.raises(SystemExit):
            validate_igv_files(str(empty_fa), None, exit_on_error=True)

        fa = tmp_path / "igv_ok.fa"
        fa.write_text(">chr1\nACGT\n")
        with pytest.raises(SystemExit):
            validate_igv_files(str(fa), None, exit_on_error=True)

        empty_fai = tmp_path / "igv_ok.fa.fai"
        empty_fai.touch()
        with pytest.raises(SystemExit):
            validate_igv_files(str(fa), None, exit_on_error=True)

        empty_fai.write_text("chr1\t4\t6\t4\t5\n")
        with pytest.raises(SystemExit):
            validate_igv_files(str(fa), "missing.ideo", exit_on_error=True)

        empty_ideo = tmp_path / "igv_empty.ideo"
        empty_ideo.touch()
        with pytest.raises(SystemExit):
            validate_igv_files(str(fa), str(empty_ideo), exit_on_error=True)
