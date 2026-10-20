"""Association genotype-matrix boundary regressions (#175 multi-allelic rejection)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine_workers import run_sequential_genes
from variantcentrifuge.association.genotype_matrix import (
    _vectorized_parse_gt_columns,
    build_genotype_matrix,
    parse_gt_to_dosage,
)
from variantcentrifuge.association.meta.artifact import emit_stratum_artifact
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

UNSPLIT = ["0/2", "1/2", "2/2", "2|1", "1|2", "2", "0/3", "./2", "2/.", ".|2", "0/1/2"]


def _site_df(gt: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "GENE": ["GENE_A", "GENE_A"],
            "s0": ["0/1", "0/0"],
            "s1": ["0/0", gt],
            "s2": ["0/0", "0/0"],
        }
    )


def _builder(df: pd.DataFrame, samples: list[str]) -> _GenotypeMatrixBuilder:
    return _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=samples,
        gt_columns=samples,
        is_binary=True,
        missing_site_threshold=0.5,
        missing_sample_threshold=0.8,
        phenotype_vector=np.array([1.0, 0.0, 0.0]),
        covariate_matrix=None,
    )


@pytest.mark.unit
@pytest.mark.parametrize("gt", UNSPLIT)
def test_unsplit_multiallelic_rejected_with_site_and_hint(gt: str) -> None:
    """Allele indices >1 never become a dosage; the error names gene, site, column."""
    df = _site_df(gt)
    with pytest.raises(ValueError) as exc:
        build_genotype_matrix(df, ["s0", "s1", "s2"], ["s0", "s1", "s2"])
    msg = str(exc.value)
    assert repr(gt) in msg
    assert "gene GENE_A" in msg
    assert "chr1:200:C:T" in msg
    assert "'s1'" in msg
    assert "bcftools norm -m-" in msg


@pytest.mark.unit
@pytest.mark.parametrize("gt", ["./2", "2/.", ".|2", "0/1/2", "2|1"])
def test_scalar_parser_rejects_partially_missing_multiallelic(gt: str) -> None:
    """A missing allele does not hide an observed second-ALT allele."""
    with pytest.raises(ValueError, match="bcftools norm -m-"):
        parse_gt_to_dosage(gt)


@pytest.mark.unit
def test_partial_biallelic_and_polyploid_without_second_alt_keep_semantics() -> None:
    """'./1' stays missing; a polyploid call without index >1 stays unparsed (NaN)."""
    assert parse_gt_to_dosage("./1") == (None, False)
    assert parse_gt_to_dosage("1/1/1") == (None, False)
    df = pd.DataFrame({"a": ["./1", "1/1/1", "0/1"], "b": ["0/0", "0/0", "0/0"]})
    raw = _vectorized_parse_gt_columns(df, ["a", "b"])[0]
    np.testing.assert_array_equal(raw[:, 0], [np.nan, np.nan, 1.0])


@pytest.mark.unit
def test_unsplit_multiallelic_without_site_columns_names_row() -> None:
    """Without CHROM/POS/REF/ALT the error falls back to the row label."""
    df = pd.DataFrame({"a": ["0/0", "2/2"], "b": ["0/1", "0/0"]}, index=["x", "y"])
    with pytest.raises(ValueError, match=r"'2/2' \(variant row 'y', column 'a'\)"):
        build_genotype_matrix(df, ["a", "b"], ["a", "b"])


@pytest.mark.unit
def test_categorical_multiallelic_column_rejected() -> None:
    df = pd.DataFrame({"a": pd.Categorical(["0/1", "1/2"]), "b": ["0/0", "0/0"]})
    with pytest.raises(ValueError, match="bcftools norm"):
        build_genotype_matrix(df, ["a", "b"], ["a", "b"])


@pytest.mark.unit
def test_normalized_inputs_keep_current_results() -> None:
    """Haploid, phased and missing calls retain their dosage/ploidy semantics."""
    df = pd.DataFrame(
        {
            "a": ["0/1", "1", "0|1", "1/1"],
            "b": ["1|0", "0", "./.", "0/0"],
            "c": ["0/0", ".", ".|.", "0/0"],
            "d": ["0/0", "0", "0/0", "."],
        }
    )
    geno, mafs, mask, warns = build_genotype_matrix(
        df, list(df), list(df), is_binary=False, missing_site_threshold=1.0
    )
    assert warns == []
    assert mask == [True] * 4
    # Row 0: AC=2, AN=8 -> 0.25. Row 1 haploid: AC=1, AN=3 -> 1/3; '.' imputed 1/3.
    # Row 2: AC=1, AN=4 -> 0.25; './.', '.|.' imputed 0.5. Row 3: haploid '.' in a
    # diploid row keeps ploidy 1: AC=2, AN=6 -> 1/3, imputed 1/3.
    np.testing.assert_allclose(mafs, [0.25, 1 / 3, 0.25, 1 / 3])
    np.testing.assert_allclose(geno[:, 0], [1, 1, 0, 0])
    np.testing.assert_allclose(geno[:, 1], [1, 0, 1 / 3, 0])
    np.testing.assert_allclose(geno[:, 2], [1, 0.5, 0.5, 0])
    np.testing.assert_allclose(geno[:, 3], [2, 0, 0, 1 / 3])


@pytest.mark.unit
def test_multiallelic_error_surfaces_through_engine_and_artifact(tmp_path) -> None:
    """The lazy builder error is not swallowed per gene by the engine or score export."""
    samples = ["s0", "s1", "s2"]
    df = _site_df("1/2")
    gene_data = {"GENE": "GENE_A", "_genotype_matrix_builder": _builder(df, samples)}
    with pytest.raises(ValueError, match="bcftools norm"):
        run_sequential_genes([gene_data], {}, AssociationConfig(), {})
    with pytest.raises(ValueError, match="bcftools norm"):
        emit_stratum_artifact(
            str(tmp_path / "s"),
            {"GENE_A": _builder(df, samples)},
            stratum_id="s",
            genome_build="GRCh38",
            annotation_version="v",
            mask_id="m",
            mask_hash="h",
        )
    assert not (tmp_path / "s.npz").exists()


@pytest.mark.unit
@pytest.mark.parametrize(
    ("gt", "expected"),
    [("bad", None), ("0/bad", None), ("bad/.", None), ("+2", "raise"), ("0/+2", "raise")],
)
def test_scalar_parser_malformed_and_signed_alleles(gt: str, expected) -> None:
    """Non-integer alleles are missing; int()-parsable alleles >1 still raise."""
    if expected == "raise":
        with pytest.raises(ValueError, match="bcftools norm -m-"):
            parse_gt_to_dosage(gt)
    else:
        assert parse_gt_to_dosage(gt) == (expected, False)
