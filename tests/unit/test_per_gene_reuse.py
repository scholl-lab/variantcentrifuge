"""Per-gene work is done once (0.25 perf audit fixes 3 and 4).

One nuisance projection per distinct input array of a gene, and under ``--skat-robust``
one per-variant SPA per gene shared by the robust test and the score artifact.
"""

from __future__ import annotations

from unittest.mock import PropertyMock, patch

import numpy as np
import pandas as pd
import pytest

from tests.unit.test_artifact_single_pass import _genes
from variantcentrifuge.association import spa
from variantcentrifuge.association.backends import null_model as nm
from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine
from variantcentrifuge.association.meta.artifact import (
    _build_for_artifact,
    compute_gene_score,
    fit_stratum_null,
)
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
from variantcentrifuge.stages.association_stage_helpers import _score_artifact_hook
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

N = 400


def _gene(seed: int = 7) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    geno = rng.binomial(2, [0.02, 0.05, 0.1, 0.03, 0.08], size=(N, 5)).astype(np.float64)
    pheno = np.r_[np.ones(80), np.zeros(N - 80)]
    cov = rng.normal(size=(N, 2))
    return geno, pheno, cov


def _count_projections(fn) -> int:
    """Projections made by ``fn()`` (``_nuisance_basis`` runs once per ``_project``)."""
    real = nm._nuisance_basis
    calls = 0

    def counting(null_model):
        nonlocal calls
        calls += 1
        return real(null_model)

    with patch.object(nm, "_nuisance_basis", counting):
        fn()
    return calls


@pytest.mark.unit
def test_skato_gene_projects_each_distinct_array_once() -> None:
    """SKAT-O + ACAT-V: tested genotypes, weighted genotypes, burden sum -> 3 projections."""
    geno, pheno, cov = _gene()
    test = PurePythonSKATTest()
    test.check_dependencies()
    data = {
        "genotype_matrix": geno,
        "phenotype_vector": pheno,
        "covariate_matrix": cov,
        "proband_count": 80,
        "control_count": N - 80,
        "n_qualifying_variants": 5,
    }
    cfg = AssociationConfig(skat_method="SKATO")
    first = test.run("G", data, cfg)  # fits (and caches) the null model
    assert first.p_value is not None and first.extra["acat_v_p"] is not None
    result = None

    def run() -> None:
        nonlocal result
        result = test.run("G", data, cfg)

    assert _count_projections(run) == 3
    assert result is not None
    assert (result.p_value, result.extra["acat_v_p"]) == (first.p_value, first.extra["acat_v_p"])


@pytest.mark.unit
def test_artifact_gene_projects_once() -> None:
    geno, pheno, cov = _gene()
    null = fit_stratum_null(pheno, cov, "binary")
    keys = [f"1:{i}:A:G" for i in range(geno.shape[1])]
    assert _count_projections(lambda: compute_gene_score(geno, geno, null, keys)) == 1


@pytest.mark.unit
def test_robust_artifact_computes_variant_spa_once_per_gene() -> None:
    """--skat-robust with the score artifact: one per-variant SPA per tested gene."""
    genes = [g for i, g in enumerate(_genes()) if i not in (1, 3, 5, 7, 9)]  # same samples
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in genes}
    engine = AssociationEngine.from_names(
        ["skat_python"], AssociationConfig(association_workers=1, skat_robust=True)
    )
    real = spa.calibrated_variance  # runs exactly once per uncached variant_pvalues
    calls = 0

    def counting(*args):
        nonlocal calls
        calls += 1
        return real(*args)

    with patch.object(spa, "calibrated_variance", counting):
        engine.run_all(genes, gene_hook=_score_artifact_hook(builders, robust=True))
    assert calls == len(genes)


def _site_frame() -> tuple[pd.DataFrame, dict[str, _GenotypeMatrixBuilder]]:
    """One run-level frame; each gene's builder holds a slice of it and its row positions.

    Genes: mixed chrom/allele spelling with a multi-allelic site, a site failing QC,
    a duplicate canonical key (``chr1`` vs ``1``), and no variants.
    """
    rng = np.random.default_rng(5)
    samples = [f"S{i}" for i in range(N)]
    sites = {
        "A": [("chr1", 100, "a", "g"), ("1", 100, "A", "T"), ("CHR1", 120, "C", "t")],
        "B": [("2", 5, "G", "A"), ("2", 9, "T", "C"), ("2", 9, "T", "CA")],
        "DUP": [("chr1", 300, "A", "G"), ("1", 300, "a", "g")],
        "EMPTY": [],
    }
    frames, builders, start = [], {}, 0
    for gene, rows in sites.items():
        codes = rng.choice(["0/0", "0/1", "1/1"], size=(len(rows), N), p=[0.8, 0.17, 0.03])
        if gene == "B":
            codes[1, : N // 2] = "./."  # fails site QC: not kept
        part = pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT"])
        frames.append(pd.concat([part, pd.DataFrame(codes, columns=samples)], axis=1))
        builders[gene] = np.arange(start, start + len(rows))
        start += len(rows)
    df = pd.concat(frames, ignore_index=True)
    pheno = np.r_[np.ones(80), np.zeros(N - 80)]
    cov = np.random.default_rng(1).normal(size=(N, 2))
    return df, {
        gene: _GenotypeMatrixBuilder(
            gene_df=df.iloc[rows[0] : rows[-1] + 1] if len(rows) else df.iloc[0:0],
            vcf_samples=samples,
            gt_columns=samples,
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=pheno,
            covariate_matrix=cov,
            rows=rows,
        )
        for gene, rows in builders.items()
    }


@pytest.mark.unit
def test_artifact_site_index_gives_identical_records() -> None:
    """Keys from the run-level site columns equal the per-gene pandas path, bit for bit."""
    df, builders = _site_frame()
    indexed = _score_artifact_hook(builders, df=df)
    per_gene = _score_artifact_hook(builders)
    for gene in ("A", "B"):
        b = builders[gene]
        built = _build_for_artifact(b)
        no_loc = PropertyMock(side_effect=AssertionError("pandas .loc per gene"))
        with patch.object(pd.DataFrame, "loc", new=no_loc):
            sa, ka, ea = indexed(gene, b, built)
        sb, kb, eb = per_gene(gene, b, built)
        assert (sa, ka) == (sb, kb)
        assert ea.variant_keys == eb.variant_keys
        for field in ("u", "v", "ac", "an"):
            assert getattr(ea, field).tobytes() == getattr(eb, field).tobytes()
    assert indexed("A", builders["A"], _build_for_artifact(builders["A"]))[2].variant_keys == [
        "1:100:A:G",
        "1:100:A:T",
        "1:120:C:T",
    ]
    assert len(indexed("B", builders["B"], _build_for_artifact(builders["B"]))[2].variant_keys) == 2
    for hook in (indexed, per_gene):
        with pytest.raises(ValueError, match="unique"):
            hook("DUP", builders["DUP"], _build_for_artifact(builders["DUP"]))
        assert hook("EMPTY", builders["EMPTY"], {}) is None
