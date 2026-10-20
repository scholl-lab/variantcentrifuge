"""Regression tests for Issue #169: validated --from-tsv mode for rescoring and reporting."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from variantcentrifuge import cli as cli_module
from variantcentrifuge.cli import main, parse_args
from variantcentrifuge.stages import processing_stages, sample_and_pedigree_stages


class TestFromTsvExecutionAndIsolation:
    """Verify --from-tsv rescoring, final filtering, reporting, and VCF stage bypass."""

    def test_from_tsv_bypasses_vcf_and_rescores_filters_and_reports(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Fail if any VCF validation or extraction helper is called; verify TSV mode succeeds."""

        def forbid_vcf(*args: Any, **kwargs: Any) -> Any:
            raise AssertionError("VCF preparation/extraction must not be called in --from-tsv mode")

        monkeypatch.setattr(cli_module, "validate_vcf_file", forbid_vcf)
        monkeypatch.setattr(cli_module, "validate_required_tools", forbid_vcf)
        monkeypatch.setattr(sample_and_pedigree_stages, "get_vcf_samples", forbid_vcf)
        monkeypatch.setattr(processing_stages, "get_gene_bed", forbid_vcf)
        monkeypatch.setattr(processing_stages, "extract_variants", forbid_vcf)
        monkeypatch.setattr(processing_stages, "extract_fields_bcftools", forbid_vcf)

        source_tsv = tmp_path / "cohort_final.tsv"
        source_content = (
            "CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\t"
            "dbNSFP_gnomAD_exomes_AF\tdbNSFP_gnomAD_genomes_AF\tdbNSFP_CADD_phred\tGT\n"
            "chr1\t100\tA\tG\tBRCA1\tmissense_variant\tMODERATE\t0.0001\t0.0001\t28.0\tS1(0/1);S2(0/0)\n"
            "chr1\t200\tC\tT\tTP53\tsynonymous_variant\tLOW\t0.25\t0.25\t2.0\tS1(0/0);S2(0/1)\n"
        )
        source_tsv.write_text(source_content, encoding="utf-8")
        source_sha_before = hashlib.sha256(source_tsv.read_bytes()).hexdigest()

        out_dir = tmp_path / "out_rescore"
        rc = main(
            [
                "--from-tsv",
                str(source_tsv),
                "--scoring-config-path",
                "scoring/test_simple_score",
                "--final-filter",
                "IMPACT == 'MODERATE' and test_simple_score > 0.4",
                "--html-report",
                "--xlsx",
                "--output-dir",
                str(out_dir),
                "-o",
                "rescored.tsv",
            ]
        )
        assert rc == 0

        # Source file must be untouched
        assert hashlib.sha256(source_tsv.read_bytes()).hexdigest() == source_sha_before

        out_tsv = out_dir / "rescored.tsv"
        assert out_tsv.is_file()
        df = pd.read_csv(out_tsv, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["GENE"] == "BRCA1"
        assert float(df.iloc[0]["test_simple_score"]) > 0.4

        # HTML and Excel reports generated
        assert (out_dir / "report" / "index.html").is_file()
        assert (out_dir / "report" / "variants.json").is_file()
        assert (out_dir / "cohort_final.xlsx").is_file()

    def test_from_tsv_recomputes_phenotypes_inheritance_and_scores(self, tmp_path: Path) -> None:
        """Changing phenotype definitions and PED recomputes counts, inheritance, and scores."""
        source_tsv = tmp_path / "variants.tsv"
        source_tsv.write_text(
            "CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\t"
            "dbNSFP_gnomAD_exomes_AF\tdbNSFP_gnomAD_genomes_AF\tdbNSFP_CADD_phred\tGT\n"
            "chr1\t100\tA\tG\tPKD1\tstop_gained\tHIGH\t0.0\t0.0\t35.0\t"
            "Child(0/1);Father(0/0);Mother(0/0)\n",
            encoding="utf-8",
        )

        # Run 1: Child has HP:0000107 (case), Father/Mother have HP:0000001 (control)
        pheno1 = tmp_path / "pheno1.tsv"
        pheno1.write_text(
            "sample_id\tphenotype\nChild\tHP:0000107\nFather\tHP:0000001\nMother\tHP:0000001\n",
            encoding="utf-8",
        )
        out1 = tmp_path / "run1"
        rc1 = main(
            [
                "--from-tsv",
                str(source_tsv),
                "--phenotype-file",
                str(pheno1),
                "--phenotype-sample-column",
                "sample_id",
                "--phenotype-value-column",
                "phenotype",
                "--case-phenotypes",
                "HP:0000107",
                "--control-phenotypes",
                "HP:0000001",
                "--output-dir",
                str(out1),
                "-o",
                "r1.tsv",
            ]
        )
        assert rc1 == 0
        df1 = pd.read_csv(out1 / "r1.tsv", sep="\t")
        assert int(df1.iloc[0]["proband_variant_count"]) == 1
        assert int(df1.iloc[0]["control_variant_count"]) == 0
        assert "HP:0000107" in str(df1.iloc[0]["Phenotypes"])

        # Inject stale inheritance columns into r1.tsv before Run 2
        df1_stale = pd.read_csv(out1 / "r1.tsv", sep="\t")
        df1_stale["Inheritance_Pattern"] = "stale_pattern"
        df1_stale["Inheritance_Details"] = "stale_details"
        df1_stale["Inheritance_Confidence"] = "stale_confidence"
        df1_stale["Inheritance_Description"] = "stale_description"
        df1_stale["Inheritance_Samples"] = "stale_samples"
        df1_stale.to_csv(out1 / "r1.tsv", sep="\t", index=False)

        # Run 2 on r1.tsv: swap case/control phenotype definitions and add trio PED + scoring
        pheno2 = tmp_path / "pheno2.tsv"
        pheno2.write_text(
            "sample_id\tphenotype\nChild\tHP:0009999\nFather\tHP:0000107\nMother\tHP:0000107\n",
            encoding="utf-8",
        )
        ped = tmp_path / "trio.ped"
        ped.write_text(
            "FAM1\tChild\tFather\tMother\t1\t2\n"
            "FAM1\tFather\t0\t0\t1\t1\n"
            "FAM1\tMother\t0\t0\t2\t1\n",
            encoding="utf-8",
        )
        out2 = tmp_path / "run2"
        rc2 = main(
            [
                "--from-tsv",
                str(out1 / "r1.tsv"),
                "--phenotype-file",
                str(pheno2),
                "--phenotype-sample-column",
                "sample_id",
                "--phenotype-value-column",
                "phenotype",
                "--case-phenotypes",
                "HP:0000107",
                "--control-phenotypes",
                "HP:0009999",
                "--ped",
                str(ped),
                "--inheritance-mode",
                "simple",
                "--scoring-config-path",
                "scoring/nephro_candidate_score",
                "--output-dir",
                str(out2),
                "-o",
                "r2.tsv",
            ]
        )
        assert rc2 == 0
        df2 = pd.read_csv(out2 / "r2.tsv", sep="\t")
        # Now Child is in control_phenotypes (HP:0009999), so proband_variant_count=0, control=1
        assert int(df2.iloc[0]["proband_variant_count"]) == 0
        assert int(df2.iloc[0]["control_variant_count"]) == 1
        assert "HP:0009999" in str(df2.iloc[0]["Phenotypes"])
        assert df2.iloc[0]["Inheritance_Pattern"] == "de_novo"
        for stale_col in (
            "Inheritance_Details",
            "Inheritance_Confidence",
            "Inheritance_Description",
            "Inheritance_Samples",
        ):
            expected_stale = f"stale_{stale_col.split('_')[1].lower()}"
            assert stale_col not in df2.columns or df2.iloc[0][stale_col] != expected_stale
        assert "nephro_candidate_score" in df2.columns
        # No duplicated link columns like ClinVar.1
        assert not any(c.endswith(".1") for c in df2.columns)

    def test_from_tsv_packed_gt_regression_association_tests(self, tmp_path: Path) -> None:
        """--from-tsv with packed GT supports regression association tests (logistic_burden)."""
        cases = [f"C{i}" for i in range(1, 11)]
        controls = [f"N{i}" for i in range(1, 11)]
        # Variant 1: 8 cases (0/1), 1 control (0/1), rest 0/0 (omitted or explicit)
        gt_v1_parts = (
            [f"{s}(0/1)" for s in cases[:8]]
            + [f"{controls[0]}(0/1)"]
            + [f"{s}(0/0)" for s in cases[8:] + controls[1:]]
        )
        # Variant 2: 6 cases (0/1), 1 control (0/1), non-carriers omitted from packed GT
        gt_v2_parts = [f"{s}(0/1)" for s in cases[2:8]] + [f"{controls[1]}(0/1)"]

        packed_tsv = tmp_path / "packed_assoc.tsv"
        packed_tsv.write_text(
            "CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tGT\n"
            f"chr1\t100\tA\tG\tBRCA1\tmissense_variant\tMODERATE\t{';'.join(gt_v1_parts)}\n"
            f"chr1\t200\tC\tT\tBRCA1\tstop_gained\tHIGH\t{';'.join(gt_v2_parts)}\n",
            encoding="utf-8",
        )

        out_dir = tmp_path / "assoc_out"
        rc = main(
            [
                "--from-tsv",
                str(packed_tsv),
                "--perform-association",
                "--association-tests",
                "fisher,logistic_burden",
                "--case-samples",
                ",".join(cases),
                "--control-samples",
                ",".join(controls),
                "--output-dir",
                str(out_dir),
                "-o",
                "assoc_variants.tsv",
                "--no-stats",
            ]
        )
        assert rc == 0
        assoc_files = list(out_dir.glob("*association*"))
        assert len(assoc_files) == 1
        assoc_df = pd.read_csv(assoc_files[0], sep="\t")
        assert len(assoc_df) == 1
        assert assoc_df.iloc[0]["gene"] == "BRCA1"
        assert int(assoc_df.iloc[0]["n_cases"]) == 10
        assert int(assoc_df.iloc[0]["n_controls"]) == 10
        assert "logistic_burden_pvalue" in assoc_df.columns
        assert pd.notna(assoc_df.iloc[0]["logistic_burden_pvalue"])
        assert 0.0 <= float(assoc_df.iloc[0]["logistic_burden_pvalue"]) <= 1.0


class TestFromTsvSchemasCompressionAndEdgeCases:
    """Cover extracted vs final TSV schemas, .tsv.gz, empty TSVs, and validation errors."""

    def test_extracted_schema_and_gzipped_input(self, tmp_path: Path) -> None:
        """Extracted TSV schema (ANN[0].GENE, GEN[0].GT) in .tsv.gz works with PED metadata."""
        extracted_gz = tmp_path / "extracted.tsv.gz"
        content = (
            "CHROM\tPOS\tREF\tALT\tANN[0].GENE\tANN[0].EFFECT\tANN[0].IMPACT\t"
            "GEN[0].GT\tGEN[1].GT\tGEN[2].GT\n"
            "chr1\t500\tG\tA\tCOL4A5\tmissense_variant\tMODERATE\t0/1\t0/0\t0/0\n"
        )
        with gzip.open(extracted_gz, "wt", encoding="utf-8") as f:
            f.write(content)

        ped = tmp_path / "trio.ped"
        ped.write_text(
            "FAM1\tProband\tDad\tMom\t1\t2\nFAM1\tDad\t0\t0\t1\t1\nFAM1\tMom\t0\t0\t2\t1\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out_extracted"
        rc = main(
            [
                "--from-tsv",
                str(extracted_gz),
                "--ped",
                str(ped),
                "--inheritance-mode",
                "simple",
                "--output-dir",
                str(out_dir),
                "-o",
                "from_extracted.tsv",
            ]
        )
        assert rc == 0
        df = pd.read_csv(out_dir / "from_extracted.tsv", sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["GENE"] == "COL4A5"
        assert "Proband(0/1)" in str(df.iloc[0]["GT"])
        assert df.iloc[0]["Inheritance_Pattern"] == "de_novo"

    def test_empty_tsv_and_fail_on_empty_in_from_tsv_mode(self, tmp_path: Path) -> None:
        """Empty (header-only) TSV succeeds normally and exits with code 2 under --fail-on-empty."""
        empty_tsv = tmp_path / "empty.tsv"
        empty_tsv.write_text(
            "CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tGT\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out_empty"
        rc_ok = main(
            [
                "--from-tsv",
                str(empty_tsv),
                "--output-dir",
                str(out_dir),
                "-o",
                "empty_out.tsv",
            ]
        )
        assert rc_ok == 0
        assert (out_dir / "empty_out.tsv").is_file()

        rc_fail = main(
            [
                "--from-tsv",
                str(empty_tsv),
                "--fail-on-empty",
                "--output-dir",
                str(tmp_path / "out_empty_fail"),
                "-o",
                "empty_fail.tsv",
            ]
        )
        assert rc_fail == 2

    def test_ambiguous_sample_mappings_and_invalid_inputs_fail_clearly(
        self, tmp_path: Path
    ) -> None:
        """Ambiguous sample mappings, missing columns, and VCF-only flags fail clearly."""
        # 1. Both -v and --from-tsv rejected by parser
        with pytest.raises(SystemExit):
            parse_args(["-v", "a.vcf.gz", "--from-tsv", "b.tsv"])

        # 2. Missing --from-tsv file
        assert main(["--from-tsv", str(tmp_path / "nonexistent.tsv")]) != 0

        # 3. Missing required columns (no GENE or CHROM)
        bad_cols_tsv = tmp_path / "bad_cols.tsv"
        bad_cols_tsv.write_text("FOO\tBAR\n1\t2\n", encoding="utf-8")
        assert main(["--from-tsv", str(bad_cols_tsv), "--output-dir", str(tmp_path / "o1")]) != 0

        # 4. Ambiguous anonymous GEN[i].GT columns when PED sample count mismatches GEN column count
        anon_tsv = tmp_path / "anon.tsv"
        anon_tsv.write_text(
            "CHROM\tPOS\tREF\tALT\tANN[0].GENE\tGEN[0].GT\tGEN[1].GT\n"
            "chr1\t100\tA\tG\tBRCA1\t0/1\t0/0\n",
            encoding="utf-8",
        )
        ped_3_samples = tmp_path / "three.ped"
        ped_3_samples.write_text(
            "FAM1\tS1\tS2\tS3\t1\t2\nFAM1\tS2\t0\t0\t1\t1\nFAM1\tS3\t0\t0\t2\t1\n",
            encoding="utf-8",
        )
        assert (
            main(
                [
                    "--from-tsv",
                    str(anon_tsv),
                    "--ped",
                    str(ped_3_samples),
                    "--output-dir",
                    str(tmp_path / "o2"),
                ]
            )
            != 0
        )

        # 5. Disjoint sample names between packed GT column and supplied PED file
        packed_tsv = tmp_path / "packed.tsv"
        packed_tsv.write_text(
            "CHROM\tPOS\tREF\tALT\tGENE\tGT\nchr1\t100\tA\tG\tBRCA1\tSampleA(0/1);SampleB(0/0)\n",
            encoding="utf-8",
        )
        disjoint_ped = tmp_path / "disjoint.ped"
        disjoint_ped.write_text("FAM1\tOtherX\t0\t0\t1\t2\n", encoding="utf-8")
        assert (
            main(
                [
                    "--from-tsv",
                    str(packed_tsv),
                    "--ped",
                    str(disjoint_ped),
                    "--output-dir",
                    str(tmp_path / "o3"),
                ]
            )
            != 0
        )

        # 6. Unsupported VCF-only flag (--bcftools-prefilter) with --from-tsv
        assert (
            main(
                [
                    "--from-tsv",
                    str(packed_tsv),
                    "--bcftools-prefilter",
                    "QUAL > 30",
                    "--output-dir",
                    str(tmp_path / "o4"),
                ]
            )
            != 0
        )

    def test_from_tsv_helper_and_validation_branches(self, tmp_path: Path) -> None:
        """Cover remaining validation and sample inference branches in from_tsv.py."""
        import argparse
        from types import SimpleNamespace

        from variantcentrifuge.from_tsv import (
            compute_tsv_base_name,
            infer_and_validate_from_tsv_samples,
            prepare_from_tsv_context,
            validate_from_tsv_args_and_file,
        )
        from variantcentrifuge.pipeline_core.context import PipelineContext
        from variantcentrifuge.pipeline_core.workspace import Workspace

        assert compute_tsv_base_name("/path/to/cohort.tsv.gz") == "cohort"
        assert compute_tsv_base_name("/path/to/cohort.tsv") == "cohort"

        # Anonymous GEN[i].GT without metadata succeeds when no sample-dependent stage is active
        df_anon = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": ["10"],
                "REF": ["A"],
                "ALT": ["T"],
                "GENE": ["G1"],
                "GEN[0].GT": ["0/1"],
            }
        )
        assert infer_and_validate_from_tsv_samples(df_anon, {}) == ["SAMPLE_0"]

        # Anonymous GEN[i].GT without metadata fails when inheritance is requested
        with pytest.raises(ValueError, match="require sample metadata"):
            infer_and_validate_from_tsv_samples(df_anon, {"calculate_inheritance": True})

        # Named sample columns with case_samples / control_samples_file
        ctrl_file = tmp_path / "controls.txt"
        ctrl_file.write_text("S2\n", encoding="utf-8")
        df_named = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": ["10"],
                "REF": ["A"],
                "ALT": ["T"],
                "GENE": ["G1"],
                "S1": ["0/1"],
                "S2": ["0/0"],
            }
        )
        samples_named = infer_and_validate_from_tsv_samples(
            df_named, {"case_samples": ["S1"], "control_samples_file": str(ctrl_file)}
        )
        assert samples_named == ["S1", "S2"]

        # Missing genotype column fails
        df_no_gt = pd.DataFrame(
            {"CHROM": ["1"], "POS": ["10"], "REF": ["A"], "ALT": ["T"], "GENE": ["G1"]}
        )
        with pytest.raises(ValueError, match="missing genotype column"):
            infer_and_validate_from_tsv_samples(df_no_gt, {})

        # 0-byte TSV, directory path, and --pca akt fail
        zero_tsv = tmp_path / "zero.tsv"
        zero_tsv.write_text("", encoding="utf-8")
        with pytest.raises(ValueError, match="empty"):
            validate_from_tsv_args_and_file(argparse.Namespace(from_tsv=str(zero_tsv)), {})

        with pytest.raises(ValueError, match="not a regular file"):
            validate_from_tsv_args_and_file(argparse.Namespace(from_tsv=str(tmp_path)), {})

        with pytest.raises(ValueError, match="--pca akt"):
            validate_from_tsv_args_and_file(
                argparse.Namespace(from_tsv=str(zero_tsv), pca="akt"), {}
            )

        # prepare_from_tsv_context with remove_sample_substring and checkpoint restore
        tsv_sub = tmp_path / "sub.tsv"
        tsv_sub.write_text(
            "CHROM\tPOS\tREF\tALT\tGENE\tInheritance_Pattern\tGT\n"
            "chr1\t10\tA\tG\tBRCA1\told\tS1_WGS(0/1)\n",
            encoding="utf-8",
        )
        ws = Workspace(tmp_path / "ws", "sub")
        ctx = PipelineContext(
            args=SimpleNamespace(from_tsv=str(tsv_sub)),  # type: ignore[arg-type]
            config={
                "from_tsv": str(tsv_sub),
                "remove_sample_substring": "_WGS",
                "calculate_inheritance": True,
                "gzip_intermediates": False,
            },
            workspace=ws,
        )
        prepare_from_tsv_context(ctx)
        assert ctx.vcf_samples == ["S1"]
        assert ctx.extracted_tsv is not None and ctx.extracted_tsv.is_file()
