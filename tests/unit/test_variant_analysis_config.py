"""Unit tests for VariantAnalysisStage case/control configuration passing."""

# import tempfile  # Unused import
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.analysis_stages import VariantAnalysisStage


class TestVariantAnalysisConfig:
    """Test VariantAnalysisStage case/control configuration passing."""

    @pytest.fixture
    def context(self):
        """Create test context with case/control assignments."""
        context = create_test_context(
            gene_name="BRCA1",
            vcf_file="/tmp/test.vcf",
            config_overrides={
                "case_samples": ["S001", "S002"],
                "control_samples": ["S003", "S004"],
                "case_phenotypes": ["HP:0000001"],
                "control_phenotypes": [],
            },
        )

        # Mock VCF samples
        context.vcf_samples = ["S001", "S002", "S003", "S004"]

        # Create a mock DataFrame with test data
        test_df = pd.DataFrame(
            {
                "CHROM": ["1", "1", "2"],
                "POS": [1000, 2000, 3000],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "A"],
                "GT": ["S001(0/1)", "S002(1/1)", "S003(0/1)"],
            }
        )
        context.current_dataframe = test_df

        # Mock workspace for temporary files
        mock_workspace = MagicMock()
        mock_workspace.get_intermediate_path.return_value = Path("/tmp/temp_analysis.tsv")
        mock_workspace.base_name = "test_analysis"
        context.workspace = mock_workspace

        return context

    @patch("variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df")
    def test_case_control_config_passed_to_analyze_variants(self, mock_analyze_variants, context):
        """Test that case/control samples are passed to analyze_variants_df function."""
        mock_analyze_variants.return_value = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "proband_count": [2, 2],
                "control_count": [2, 2],
            }
        )

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Verify analyze_variants_df was called
        assert mock_analyze_variants.call_count == 1

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]  # Second argument is the config

        # Verify case/control samples were passed
        assert "case_samples" in analysis_config
        assert "control_samples" in analysis_config
        assert analysis_config["case_samples"] == ["S001", "S002"]
        assert analysis_config["control_samples"] == ["S003", "S004"]

        # Verify phenotype terms were also passed for legacy compatibility
        assert "case_phenotypes" in analysis_config
        assert "control_phenotypes" in analysis_config
        assert analysis_config["case_phenotypes"] == ["HP:0000001"]
        assert analysis_config["control_phenotypes"] == []

    @patch("variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df")
    def test_empty_case_control_config_passed(self, mock_analyze_variants, context):
        """Test that empty case/control samples are passed correctly."""
        # Set empty case/control samples
        context.config["case_samples"] = []
        context.config["control_samples"] = []

        mock_analyze_variants.return_value = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
                "proband_count": [0],
                "control_count": [4],
            }
        )

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]

        # Verify empty lists were passed
        assert analysis_config["case_samples"] == []
        assert analysis_config["control_samples"] == []

    @patch("variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df")
    def test_missing_case_control_config_defaults(self, mock_analyze_variants, context):
        """Test that missing case/control samples default to empty lists."""
        # Remove case/control samples from config
        if "case_samples" in context.config:
            del context.config["case_samples"]
        if "control_samples" in context.config:
            del context.config["control_samples"]

        mock_analyze_variants.return_value = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
            }
        )

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]

        # Verify defaults to empty lists
        assert analysis_config["case_samples"] == []
        assert analysis_config["control_samples"] == []

    @patch("variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df")
    def test_vcf_samples_config_passed(self, mock_analyze_variants, context):
        """Test that VCF samples and other config are passed correctly."""
        mock_analyze_variants.return_value = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
            }
        )

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]

        # Verify VCF samples and other required config
        assert analysis_config["vcf_sample_names"] == ["S001", "S002", "S003", "S004"]
        assert analysis_config["sample_list"] == "S001,S002,S003,S004"
        assert analysis_config["base_name"] == "test_analysis"
        assert "reference" in analysis_config

    def test_analyze_variants_df_direct(self):
        """Test analyze_variants_df directly in-memory."""
        from variantcentrifuge.analyze_variants import analyze_variants_df

        df = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["BRCA1", "BRCA1"],
                "GT": ["S001(0/1);S002(1/1)", "S003(0/1)"],
            }
        )
        cfg = {
            "sample_list": "S001,S002,S003,S004",
            "case_samples": ["S001", "S002"],
            "control_samples": ["S003", "S004"],
            "no_stats": True,
        }
        res = analyze_variants_df(df, cfg)
        assert "proband_count" in res.columns
        assert "control_count" in res.columns
        assert res["proband_count"].iloc[0] == 2
        assert res["control_count"].iloc[0] == 2
        assert res["proband_variant_count"].iloc[0] == 2
        assert res["control_variant_count"].iloc[0] == 0

    def test_variant_analysis_stage_in_memory(self, context, tmp_path):
        """Test VariantAnalysisStage runs in-memory without temp files or duplicate gene burden."""
        context.workspace.get_intermediate_path.return_value = tmp_path / "temp_analysis.tsv"
        context.current_dataframe["GENE"] = ["BRCA1", "BRCA1", "TP53"]
        context.current_dataframe["Inheritance_Pattern"] = ["AD", "AR", "XL"]
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        out_context = stage(context)

        # Temp file must not have been created on disk
        assert not (tmp_path / "temp_analysis.tsv").exists()
        assert not (tmp_path / "temp_analysis.tsv.gz").exists()

        # Output DataFrame has counts and preserved inheritance
        df = out_context.current_dataframe
        assert "proband_count" in df.columns
        assert "control_count" in df.columns
        assert "Inheritance_Pattern" in df.columns
        assert list(df["Inheritance_Pattern"]) == ["AD", "AR", "XL"]

    def test_variant_analysis_stage_none_dataframe(self, context):
        """Test VariantAnalysisStage gracefully handles None current_dataframe."""
        context.current_dataframe = None
        context.mark_complete("dataframe_loading")
        stage = VariantAnalysisStage()
        out = stage(context)
        assert out.current_dataframe is None

    def test_variant_analysis_stage_gt_reconstruction_and_renaming(self, context):
        """Test GT reconstruction when per-sample GT columns exist and rename map is used."""
        df = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
                "GENE": ["BRCA1"],
                "GEN_0__GT": ["0/1"],
                "GEN_1__GT": ["0/0"],
                "Extra_Col": ["extra_val"],
                "VAR_ID": ["1:1000:A:T"],
                "Custom_Annotation": ["ann"],
            }
        )
        context.current_dataframe = df
        context.vcf_samples = ["S001", "S002"]
        context.column_rename_map = {"ANN[0].FEATUREID": "ANN_0__FEATUREID"}
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        out = stage(context)
        assert out.current_dataframe is not None
        assert "Extra_Col" in out.current_dataframe.columns
        assert "GT" in out.current_dataframe.columns

    @patch(
        "variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df", return_value=None
    )
    def test_variant_analysis_stage_analysis_df_none_and_properties(self, mock_analyze, context):
        """Test handling when analyze_variants_df returns None and check stage properties."""
        context.mark_complete("dataframe_loading")
        stage = VariantAnalysisStage()
        assert stage.name == "variant_analysis"
        assert stage.description == "Perform variant-level analysis"
        assert stage.dependencies == {"dataframe_loading"}
        assert "inheritance_analysis" in stage.soft_dependencies
        assert stage.get_input_files(context) == []
        assert stage.get_output_files(context) == []

        out = stage(context)
        assert out.current_dataframe is not None

    @patch("variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df")
    def test_variant_analysis_stage_merge_length_mismatch(self, mock_analyze, context):
        """Test preservation merge when analysis_df has different length than input."""
        mock_analyze.return_value = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
                "proband_count": [1],
            }
        )
        context.current_dataframe["Extra"] = ["e1", "e2", "e3"]
        context.mark_complete("dataframe_loading")
        stage = VariantAnalysisStage()
        out = stage(context)
        assert len(out.current_dataframe) == 1
        assert "Extra" in out.current_dataframe.columns

    @patch("variantcentrifuge.stages.variant_analysis_stages.analyze_variants_df")
    def test_variant_analysis_stage_preserves_inheritance_columns(self, mock_analyze, context):
        """Test preserving inheritance columns when original_only_cols has Inheritance_Pattern."""
        mock_analyze.return_value = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "proband_count": [1, 2],
            }
        )
        context.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "Inheritance_Pattern": ["AD", "AR"],
                "VAR_ID": ["v1", "v2"],
                "Custom_Annotation": ["c1", "c2"],
                "GT": ["g1", "g2"],
            }
        )
        context.mark_complete("dataframe_loading")
        stage = VariantAnalysisStage()
        out = stage(context)
        assert "Inheritance_Pattern" in out.current_dataframe.columns
        assert list(out.current_dataframe["Inheritance_Pattern"]) == ["AD", "AR"]
        assert "VAR_ID" in out.current_dataframe.columns

    def test_variant_analysis_buffer_sharing(self, context):
        """Test that analysis stages and helpers maintain buffer sharing via shallow copies."""
        import numpy as np

        from variantcentrifuge.analyze_variants import analyze_variants_df
        from variantcentrifuge.gt_reconstruction import reconstruct_gt_column

        # 1. Test reconstruct_gt_column buffer sharing
        df_gt = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GEN[0].GT": ["0/1", "1/1"],
                "GEN[1].GT": ["0/0", "0/1"],
            }
        )
        reconstructed = reconstruct_gt_column(df_gt, ["S001", "S002"])
        assert np.shares_memory(df_gt["POS"].to_numpy(), reconstructed["POS"].to_numpy())
        assert "GT" in reconstructed.columns

        # 2. Test analyze_variants_df buffer sharing
        df_analysis = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["BRCA1", "BRCA1"],
                "GT": ["S001(0/1)", "S001(1/1);S002(0/1)"],
            }
        )
        cfg = {
            "sample_list": "S001,S002",
            "case_samples": ["S001"],
            "control_samples": ["S002"],
            "no_stats": True,
        }
        analyzed = analyze_variants_df(df_analysis, cfg)
        assert np.shares_memory(df_analysis["POS"].to_numpy(), analyzed["POS"].to_numpy())

        # 3. Test VariantAnalysisStage preserves buffer sharing for original columns
        context.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["BRCA1", "BRCA1"],
                "GEN_0__GT": ["0/1", "1/1"],
                "GEN_1__GT": ["0/0", "0/1"],
                "Extra_Col": [10, 20],
            }
        )
        context.vcf_samples = ["S001", "S002"]
        context.column_rename_map = {"GEN[0].GT": "GEN_0__GT", "GEN[1].GT": "GEN_1__GT"}
        context.mark_complete("dataframe_loading")
        orig_extra = context.current_dataframe["Extra_Col"].to_numpy()

        stage = VariantAnalysisStage()
        out = stage(context)
        assert out.current_dataframe is not None
        assert np.shares_memory(orig_extra, out.current_dataframe["Extra_Col"].to_numpy())

    def test_analyze_variants_burden_contract(self, tmp_path):
        """Test analyze_variants gene burden yielding and sidecar output contract."""
        from variantcentrifuge.analyze_variants import analyze_variants

        tsv_data = (
            "CHROM\tPOS\tREF\tALT\tGENE\tGT\n"
            "1\t1000\tA\tT\tBRCA1\tS001(0/1)\n"
            "1\t2000\tG\tC\tBRCA1\tS001(1/1);S002(0/1)\n"
        )
        burden_out = tmp_path / "custom_burden.tsv"
        cfg = {
            "sample_list": "S001,S002",
            "case_samples": ["S001"],
            "control_samples": ["S002"],
            "no_stats": True,
            "perform_gene_burden": True,
            "gene_burden_output_file": str(burden_out),
            "gene_burden_mode": "samples",
        }

        output_lines = list(analyze_variants(iter(tsv_data.splitlines(keepends=True)), cfg))
        assert any("GENE" in line and "p_value" in line for line in output_lines)
        assert not any("===BEGIN_GENE_BURDEN===" in line for line in output_lines)
        assert burden_out.exists()

        cfg["perform_gene_burden"] = False
        var_lines = list(analyze_variants(iter(tsv_data.splitlines(keepends=True)), cfg))
        assert var_lines[0].startswith("CHROM\tPOS\tREF\tALT")
        assert "proband_count" in var_lines[0]

    def test_analyze_variants_with_inheritance_and_scoring(self):
        """Test analyze_variants computes inheritance and scoring columns when configured."""
        from variantcentrifuge.analyze_variants import analyze_variants

        tsv_data = (
            "CHROM\tPOS\tREF\tALT\tGENE\tGT\tfather\tmother\tchild\n"
            "1\t1000\tA\tT\tBRCA1\tchild(0/1)\t0/0\t0/0\t0/1\n"
        )
        trio_pedigree = {
            "father": {
                "family_id": "FAM1",
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
            "mother": {
                "family_id": "FAM1",
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "child": {
                "family_id": "FAM1",
                "sample_id": "child",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "1",
                "affected_status": "2",
            },
        }
        scoring_cfg = {
            "variables": {"GENE": "g|default:unknown"},
            "formulas": [{"gene_score": "((g) == 'BRCA1') * 42"}],
        }
        cfg = {
            "sample_list": "father,mother,child",
            "case_samples": ["child"],
            "control_samples": ["father", "mother"],
            "no_stats": True,
            "calculate_inheritance": True,
            "pedigree_data": trio_pedigree,
            "scoring_config": scoring_cfg,
        }

        lines = list(analyze_variants(iter(tsv_data.splitlines(keepends=True)), cfg))
        header = lines[0].split("\t")
        assert "Inheritance_Pattern" in header
        assert "Inheritance_Details" in header
        assert "gene_score" in header

        data_row = lines[1].split("\t")
        pattern_idx = header.index("Inheritance_Pattern")
        score_idx = header.index("gene_score")
        assert "de_novo" in data_row[pattern_idx].lower() or "ad" in data_row[pattern_idx].lower()
        assert float(data_row[score_idx]) == 42.0
