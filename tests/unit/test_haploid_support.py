"""Regression tests for haploid genotype support (Issue #137)."""

import subprocess
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.genotype_matrix import (
    build_genotype_matrix,
    parse_gt_to_dosage,
)
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
from variantcentrifuge.extractor import extract_fields_bcftools
from variantcentrifuge.gene_burden import _aggregate_gene_burden_from_columns
from variantcentrifuge.generate_igv_report import generate_igv_report
from variantcentrifuge.gt_reconstruction import reconstruct_gt_column
from variantcentrifuge.inheritance.analyzer import analyze_inheritance
from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
)
from variantcentrifuge.inheritance.deducer import (
    check_x_linked_recessive,
    deduce_patterns_for_sample,
)
from variantcentrifuge.inheritance.vectorized_deducer import (
    vectorized_deduce_patterns,
)
from variantcentrifuge.phenotype_helpers import assign_case_control_counts


def test_parse_gt_to_dosage_haploid():
    """Test parse_gt_to_dosage with haploid calls."""
    assert parse_gt_to_dosage("0") == (0, False)
    assert parse_gt_to_dosage("1") == (1, False)
    assert parse_gt_to_dosage(".") == (None, False)


def test_genotype_matrix_haploid_calls_and_ploidy_imputation():
    """Test that matrix calls ['1','0/0'] have raw dosages [1,0] and freq 1/3."""
    gene_df = pd.DataFrame(
        {
            "CHROM": ["chrX"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["T"],
            "GEN_0__GT": ["1"],
            "GEN_1__GT": ["0/0"],
        }
    )
    vcf_samples = ["male_carrier", "female_ref"]
    gt_cols = ["GEN_0__GT", "GEN_1__GT"]

    geno, mafs, _, _ = build_genotype_matrix(
        gene_df=gene_df,
        vcf_samples=vcf_samples,
        gt_columns=gt_cols,
        is_binary=False,
    )

    # First site: ['1', '0/0']
    # Not excluded due to missingness
    assert len(mafs) == 1
    # Frequency: AC=1 / AN=(1 + 2) = 1/3
    assert pytest.approx(mafs[0], rel=1e-5) == 1.0 / 3.0

    # Dosages of site 0: male=1.0, female=0.0
    # Note: geno has shape (n_samples, n_kept_variants)
    assert geno[0, 0] == 1.0
    assert geno[1, 0] == 0.0

    # Imputation should respect ploidy
    impute_df = pd.DataFrame(
        {
            "CHROM": ["chrX"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["T"],
            "s_obs": ["0/1"],
            "s_hap_miss": ["."],
            "s_dip_miss": ["./."],
        }
    )
    geno_imp, _, _, _ = build_genotype_matrix(
        gene_df=impute_df,
        vcf_samples=["s_obs", "s_hap_miss", "s_dip_miss"],
        gt_columns=["s_obs", "s_hap_miss", "s_dip_miss"],
        is_binary=False,
        missing_site_threshold=0.8,
    )
    # Observed is 0/1: AC=1, AN=2 -> MAF = 0.5
    # s_hap_miss (ploidy 1): 1.0 * 0.5 = 0.5
    # s_dip_miss (ploidy 2): 2.0 * 0.5 = 1.0
    assert pytest.approx(geno_imp[0, 0]) == 1.0  # observed het
    assert pytest.approx(geno_imp[1, 0]) == 0.5  # haploid missing imputed to 0.5
    assert pytest.approx(geno_imp[2, 0]) == 1.0  # diploid missing imputed to 1.0


def test_packed_and_persample_counts_for_haploid():
    """Regression test: packed and per-sample counting paths report 1 carrier and 1 allele."""
    # 1. Packed GT path
    df_packed = pd.DataFrame(
        {
            "CHROM": ["chrX"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["T"],
            "GT": ["case_1(1);control_1(0)"],
        }
    )
    case_samples = {"case_1"}
    control_samples = {"control_1"}
    all_samples = case_samples | control_samples

    res_packed = assign_case_control_counts(df_packed, case_samples, control_samples, all_samples)
    assert res_packed["proband_variant_count"].iloc[0] == 1
    assert res_packed["proband_allele_count"].iloc[0] == 1
    assert res_packed["control_variant_count"].iloc[0] == 0
    assert res_packed["control_allele_count"].iloc[0] == 0

    # 2. Per-sample column path
    df_columns = pd.DataFrame(
        {
            "CHROM": ["chrX"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["GENE1"],
            "GEN_0__GT": ["1"],
            "GEN_1__GT": ["0"],
        }
    )
    vcf_samples = ["case_1", "control_1"]
    gt_columns = ["GEN_0__GT", "GEN_1__GT"]

    burden = _aggregate_gene_burden_from_columns(
        df=df_columns,
        case_samples=case_samples,
        control_samples=control_samples,
        vcf_samples=vcf_samples,
        gt_columns=gt_columns,
    )
    assert burden[0]["proband_carrier_count"] == 1
    assert burden[0]["proband_allele_count"] == 1
    assert burden[0]["control_carrier_count"] == 0
    assert burden[0]["control_allele_count"] == 0


def test_haploid_x_linked_trio_deduction():
    """Regression test: documented haploid X-linked trio yields x_linked_recessive."""
    pedigree_data = {
        "son": {"sex": "1", "affected_status": "2", "father_id": "father", "mother_id": "mother"},
        "father": {"sex": "1", "affected_status": "1", "father_id": "0", "mother_id": "0"},
        "mother": {"sex": "2", "affected_status": "1", "father_id": "0", "mother_id": "0"},
    }
    sample_list = ["son", "father", "mother"]

    variant_row = {
        "CHROM": "chrX",
        "POS": 10000000,
        "REF": "A",
        "ALT": "G",
        "son": "1",
        "father": "0",
        "mother": "0/1",
    }

    # 1. Scalar deduction
    scalar_xlr = check_x_linked_recessive("son", variant_row, pedigree_data, sample_list)
    assert scalar_xlr == "x_linked_recessive"

    scalar_patterns = deduce_patterns_for_sample("son", variant_row, pedigree_data, sample_list)
    assert "x_linked_recessive" in scalar_patterns

    # 2. Vectorized deduction
    df = pd.DataFrame([variant_row])
    vec_patterns = vectorized_deduce_patterns(df, pedigree_data, sample_list)
    assert "x_linked_recessive" in vec_patterns[0]

    # 3. analyze_inheritance
    analyzed_df = analyze_inheritance(df, pedigree_data, sample_list)
    assert analyzed_df["Inheritance_Pattern"].iloc[0] == "x_linked_recessive"


def test_haploid_mitochondrial_transmission():
    """Regression test: maternal haploid mitochondrial calls retain mitochondrial classification."""
    pedigree_data = {
        "child": {"sex": "1", "affected_status": "2", "father_id": "father", "mother_id": "mother"},
        "father": {"sex": "1", "affected_status": "1", "father_id": "0", "mother_id": "0"},
        "mother": {"sex": "2", "affected_status": "1", "father_id": "0", "mother_id": "0"},
    }
    sample_list = ["child", "father", "mother"]

    variant_row = {
        "CHROM": "chrM",
        "POS": 500,
        "REF": "C",
        "ALT": "T",
        "child": "1",
        "father": "0",
        "mother": "1",
    }

    df = pd.DataFrame([variant_row])
    analyzed_df = analyze_inheritance(df, pedigree_data, sample_list)
    assert analyzed_df["Inheritance_Pattern"].iloc[0] == "mitochondrial"


def test_haploid_alternate_never_compound_het_candidate():
    """Regression test: haploid alternate calls ('1') never become compound-het candidates."""
    pedigree_data = {
        "sample": {"sex": "1", "affected": True, "father": "0", "mother": "0"},
    }
    sample_list = ["sample"]

    # 2 variants in same gene with haploid '1'
    gene_df = pd.DataFrame(
        {
            "CHROM": ["chrX", "chrX"],
            "POS": [1000, 2000],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
            "sample": ["1", "1"],
            "Gene": ["GENE1", "GENE1"],
        }
    )

    results = analyze_gene_for_compound_het_vectorized(
        gene_df=gene_df,
        pedigree_data=pedigree_data,
        sample_list=sample_list,
    )
    # Should not produce any compound het candidate pairs
    assert len(results) == 0


@pytest.mark.unit
def test_fixed_alt_site_across_mixed_ploidy_skat():
    """Verify fixed-ALT site across mixed ploidy (haploid controls '1', diploid cases '1/1')
    is monomorphic (geno=0, MAF=0) and does not cause spurious SKAT associations.
    """
    n_ctrl = 50
    n_case = 50
    ctrl_samples = [f"ctrl_{i}" for i in range(n_ctrl)]
    case_samples = [f"case_{i}" for i in range(n_case)]
    vcf_samples = ctrl_samples + case_samples

    # Site 0: Fixed ALT across all samples ('1' in 50 haploids, '1/1' in 50 diploids)
    # Site 1: Balanced rare variant (1 carrier in ctrl, 1 carrier in case)
    # Site 2: Balanced rare variant (1 carrier in ctrl, 1 carrier in case)
    s0_ctrl = ["1"] * n_ctrl
    s0_case = ["1/1"] * n_case

    s1_ctrl = ["1"] + ["0"] * (n_ctrl - 1)
    s1_case = ["0/1"] + ["0/0"] * (n_case - 1)

    s2_ctrl = ["0", "1"] + ["0"] * (n_ctrl - 2)
    s2_case = ["0/0", "0/1"] + ["0/0"] * (n_case - 2)

    data = {
        "CHROM": ["chrX", "chrX", "chrX"],
        "POS": [1000, 2000, 3000],
        "REF": ["A", "C", "G"],
        "ALT": ["T", "T", "T"],
    }
    for s_idx, s in enumerate(ctrl_samples):
        data[s] = [s0_ctrl[s_idx], s1_ctrl[s_idx], s2_ctrl[s_idx]]
    for s_idx, s in enumerate(case_samples):
        data[s] = [s0_case[s_idx], s1_case[s_idx], s2_case[s_idx]]

    df = pd.DataFrame(data)

    geno, mafs, _, _ = build_genotype_matrix(
        gene_df=df,
        vcf_samples=vcf_samples,
        gt_columns=vcf_samples,
        is_binary=True,
    )

    # Site 0 must be zeroed out and have MAF 0.0
    assert np.unique(geno[:, 0]).tolist() == [0.0]
    assert mafs[0] == 0.0

    # Test with PurePythonSKATTest
    phenotype = np.array([0.0] * n_ctrl + [1.0] * n_case)
    config = AssociationConfig(
        trait_type="binary",
        skat_method="SKAT",
        variant_weights="beta:1,25",
    )

    test = PurePythonSKATTest()
    test.check_dependencies()
    test.prepare(1)

    # Run on all 3 sites
    contingency_data_3 = {
        "genotype_matrix": geno,
        "variant_mafs": mafs,
        "phenotype_vector": phenotype,
        "proband_count": n_case,
        "control_count": n_ctrl,
        "n_qualifying_variants": 3,
    }
    res_3 = test.run("GENE_TEST", contingency_data_3, config)

    # Run on the 2 polymorphic sites directly
    contingency_data_2 = {
        "genotype_matrix": geno[:, 1:],
        "variant_mafs": mafs[1:],
        "phenotype_vector": phenotype,
        "proband_count": n_case,
        "control_count": n_ctrl,
        "n_qualifying_variants": 2,
    }
    res_2 = test.run("GENE_TEST", contingency_data_2, config)
    test.finalize()

    assert res_3.p_value == 1.0
    assert res_3.extra["skat_n_marker_test"] == 2
    assert res_3.extra.get("acat_v_p") is None

    assert res_2.p_value == 1.0
    assert res_2.extra["skat_n_marker_test"] == 2
    assert res_2.extra.get("acat_v_p") is None


@pytest.mark.unit
def test_reconstruct_gt_custom_delimiter_schedules_igv_report(tmp_path):
    """Verify reconstruct_gt_column with delimiter=',' and ':' both schedule 1 report
    in generate_igv_report.
    """
    bam_file = tmp_path / "bams.txt"
    bam_file.write_text("S1,/path/s1.bam\n")

    for delim in [",", ":"]:
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "GEN_0__GT": ["0/1"],
                "GEN_0__DP": ["20"],
            }
        )
        reconstructed = reconstruct_gt_column(
            df,
            vcf_samples=["S1"],
            extra_sample_fields=["GEN[*].DP"],
            delimiter=delim,
        )
        assert f"S1(0/1{delim}20)" in reconstructed["GT"].iloc[0]

        tsv_file = tmp_path / f"variants_{delim.replace(':', 'colon')}.tsv"
        reconstructed.to_csv(tsv_file, sep="\t", index=False)

        with patch(
            "variantcentrifuge.generate_igv_report._generate_single_igv_report"
        ) as mock_report:
            mock_report.return_value = (True, "S1", "chr1", "100", "A", "T", "report.html")
            generate_igv_report(
                variants_tsv=str(tsv_file),
                output_dir=str(tmp_path),
                bam_mapping_file=str(bam_file),
            )
            assert len(mock_report.call_args_list) == 1
            assert mock_report.call_args_list[0][0][0][0] == "S1"


@pytest.mark.unit
def test_multiallelic_alternate_calls_consistency_across_modules():
    """Verify multiallelic non-ALT1 calls ('2', '0/2') evaluate as non-variant/0 alleles
    consistently across scalar deduction, vectorized deduction, inheritance analysis,
    gene burden, and phenotype helpers.
    """
    pedigree_data = {
        "child": {"sex": "1", "affected_status": "2", "father_id": "father", "mother_id": "mother"},
        "father": {"sex": "1", "affected_status": "1", "father_id": "0", "mother_id": "0"},
        "mother": {"sex": "2", "affected_status": "1", "father_id": "0", "mother_id": "0"},
    }
    sample_list = ["child", "father", "mother"]
    variant_row = {
        "CHROM": "chrX",
        "POS": 1000,
        "REF": "A",
        "ALT": "G",
        "child": "2",
        "father": "0",
        "mother": "0/2",
    }

    # 1. Scalar deduction
    scalar_patterns = deduce_patterns_for_sample("child", variant_row, pedigree_data, sample_list)
    assert scalar_patterns == []

    # 2. Vectorized deduction
    df = pd.DataFrame([variant_row])
    vec_patterns = vectorized_deduce_patterns(df, pedigree_data, sample_list)
    assert vec_patterns == [[]]

    # 3. Inheritance analysis
    analyzed_df = analyze_inheritance(df, pedigree_data, sample_list)
    assert analyzed_df["Inheritance_Pattern"].iloc[0] == "unknown"

    # 4. Phenotype helpers (packed GT)
    df_packed = pd.DataFrame(
        {
            "CHROM": ["chrX"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["G"],
            "GT": ["child(2);father(0);mother(0/2)"],
        }
    )
    counts = assign_case_control_counts(
        df_packed,
        case_samples={"child"},
        control_samples={"father", "mother"},
        all_samples={"child", "father", "mother"},
    )
    assert counts["proband_variant_count"].iloc[0] == 0
    assert counts["proband_allele_count"].iloc[0] == 0
    assert counts["control_variant_count"].iloc[0] == 0
    assert counts["control_allele_count"].iloc[0] == 0

    # 5. Gene burden
    df_cols = pd.DataFrame(
        {
            "CHROM": ["chrX"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["G"],
            "GENE": ["GENE1"],
            "GEN_0__GT": ["2"],
            "GEN_1__GT": ["0"],
            "GEN_2__GT": ["0/2"],
        }
    )
    burden = _aggregate_gene_burden_from_columns(
        df=df_cols,
        case_samples={"child"},
        control_samples={"father", "mother"},
        vcf_samples=["child", "father", "mother"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT"],
    )
    assert burden[0]["proband_carrier_count"] == 0
    assert burden[0]["proband_allele_count"] == 0
    assert burden[0]["control_carrier_count"] == 0
    assert burden[0]["control_allele_count"] == 0


@pytest.mark.unit
def test_extractor_to_build_genotype_matrix_haploid_support(tmp_path):
    """Verify extractor preserves haploid missing '.' and matrix imputes haploid missing."""
    out_file = str(tmp_path / "extracted.tsv")

    def fake_subprocess_run(cmd, stdout=None, **kwargs):
        # Format: CHROM POS REF ALT ID GEN[0].GT GEN[1].GT GEN[2].GT GEN[3].GT
        stdout.write("chrX\t1000\tA\tG\t.\t0\t1\t.\t\n")
        stdout.flush()
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    with patch("variantcentrifuge.extractor.subprocess.run", side_effect=fake_subprocess_run):
        extract_fields_bcftools(
            variant_file="dummy.vcf",
            fields="CHROM POS REF ALT ID GEN[*].GT",
            cfg={},
            output_file=out_file,
            vcf_samples=["s0", "s1", "s2", "s3"],
        )

    extracted_df = pd.read_csv(out_file, sep="\t", keep_default_na=False)

    # Verify ID was replaced from '.' to 'NA', but GEN[2].GT '.' was preserved
    assert extracted_df["ID"].iloc[0] == "NA"
    assert extracted_df["GEN[2].GT"].iloc[0] == "."
    assert extracted_df["GEN[3].GT"].iloc[0] == "NA"

    # In binary mode: missing haploid samples ('.' and 'NA') receive round(1.0 * 0.5) = 0.0
    geno_bin, mafs, _, _ = build_genotype_matrix(
        gene_df=extracted_df,
        vcf_samples=["s0", "s1", "s2", "s3"],
        gt_columns=["GEN[0].GT", "GEN[1].GT", "GEN[2].GT", "GEN[3].GT"],
        is_binary=True,
        missing_site_threshold=0.8,
    )
    assert mafs[0] == 0.5
    assert geno_bin[0, 0] == 0.0  # s0: obs 0
    assert geno_bin[1, 0] == 1.0  # s1: obs 1
    assert geno_bin[2, 0] == 0.0  # s2: haploid missing '.', imputed to 0.0
    assert geno_bin[3, 0] == 0.0  # s3: generic missing 'NA', imputed to 0.0

    # In quantitative mode: missing haploid samples receive 1.0 * 0.5 = 0.5 (not 1.0)
    geno_quant, _, _, _ = build_genotype_matrix(
        gene_df=extracted_df,
        vcf_samples=["s0", "s1", "s2", "s3"],
        gt_columns=["GEN[0].GT", "GEN[1].GT", "GEN[2].GT", "GEN[3].GT"],
        is_binary=False,
        missing_site_threshold=0.8,
    )
    assert geno_quant[0, 0] == 0.0
    assert geno_quant[1, 0] == 1.0
    assert geno_quant[2, 0] == 0.5
    assert geno_quant[3, 0] == 0.5


def test_coast_python_backend_and_test_haploid_variant_mafs():
    """Verify COAST python backend and PurePythonCOASTTest use ploidy-aware variant_mafs."""
    from variantcentrifuge.association.backends.coast_python import (
        PythonCOASTBackend,
        _compute_allelic_skat_weights,
    )
    from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest

    # Call with haploid+diploid samples ['1', '0/0'] -> raw dosages [1.0, 0.0]
    geno = np.array([[1.0], [0.0]])
    anno = np.array([1])  # BMV category
    coast_weights = [1.0, 2.0, 3.0]
    variant_mafs = np.array([1.0 / 3.0])

    # 1. Direct weight calculation check
    weights_with_mafs = _compute_allelic_skat_weights(
        geno, anno, coast_weights, variant_mafs=variant_mafs
    )
    expected_weight = 1.0 / np.sqrt((1.0 / 3.0) * (2.0 / 3.0))  # ~2.12132
    assert np.isclose(weights_with_mafs[0], expected_weight)

    # Naive without variant_mafs would be 1.0 / np.sqrt((1/4) * (3/4)) ~ 2.3094
    weights_naive = _compute_allelic_skat_weights(geno, anno, coast_weights, variant_mafs=None)
    assert not np.isclose(weights_naive[0], expected_weight)
    assert np.isclose(weights_naive[0], 1.0 / np.sqrt(0.25 * 0.75))

    # 2. PythonCOASTBackend.test_gene with variant_mafs passed
    backend = PythonCOASTBackend()
    backend._ensure_skat_backend()
    assert backend._skat_backend is not None
    phenotype = np.array([1, 0])
    null_model = backend._skat_backend.fit_null_model(
        phenotype=phenotype, covariates=None, trait_type="binary"
    )

    with patch(
        "variantcentrifuge.association.backends.coast_python._compute_allelic_skat_weights",
        wraps=_compute_allelic_skat_weights,
    ) as mock_weights:
        backend.test_gene(
            gene="TEST_GENE",
            geno_filtered=geno,
            anno_codes_filtered=anno,
            phenotype=phenotype,
            covariates=None,
            coast_weights=coast_weights,
            trait_type="binary",
            null_model=null_model,
            variant_mafs=variant_mafs,
        )
        mock_weights.assert_called_once()
        passed_mafs = mock_weights.call_args.kwargs.get("variant_mafs")
        assert passed_mafs is not None
        assert np.isclose(passed_mafs[0], 1.0 / 3.0)

    # 3. PurePythonCOASTTest end-to-end with contingency_data containing variant_mafs
    gene_df = pd.DataFrame(
        {
            "EFFECT": ["missense_variant"],
            "IMPACT": ["MODERATE"],
            "SIFT_pred": ["T"],
            "Polyphen2_HDIV_pred": ["B"],
        }
    )
    contingency_data = {
        "genotype_matrix": geno,
        "phenotype_vector": phenotype,
        "gene_df": gene_df,
        "variant_mafs": variant_mafs,
        "n_cases": 1,
        "n_controls": 1,
    }
    coast_test = PurePythonCOASTTest()
    coast_test.check_dependencies()
    with patch.object(
        coast_test._backend, "test_gene", wraps=coast_test._backend.test_gene
    ) as mock_test_gene:
        res = coast_test.run("TEST_GENE", contingency_data, AssociationConfig(trait_type="binary"))
        assert res.test_name == "coast"
        mock_test_gene.assert_called_once()
        passed_backend_mafs = mock_test_gene.call_args.kwargs.get("variant_mafs")
        assert passed_backend_mafs is not None
        assert np.isclose(passed_backend_mafs[0], 1.0 / 3.0)
