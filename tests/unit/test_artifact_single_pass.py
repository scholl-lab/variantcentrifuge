"""Score artifact computed in the engine's per-gene pass (#212).

The engine's workers compute each gene's artifact entry from the matrix they build for
the tests; the result must equal the artifact built gene by gene in the parent.
"""

from __future__ import annotations

import multiprocessing
from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine
from variantcentrifuge.association.meta import artifact as art_mod
from variantcentrifuge.association.meta.artifact import emit_stratum_artifact, load_artifact
from variantcentrifuge.stages.association_stage_helpers import _score_artifact_hook
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

N = 80
KW: dict[str, Any] = {
    "genome_build": "GRCh38",
    "annotation_version": "GRCh38.99",
    "mask_id": "m",
    "mask_hash": "h",
}


def _genes(n_genes: int = 12) -> list[dict[str, Any]]:
    """Synthetic genes: some samples mostly missing in a few genes (per-gene sample masks),
    one MAC<5 gene, one gene without variants."""
    rng = np.random.default_rng(3)
    samples = [f"S{i}" for i in range(N)]
    pheno = np.r_[np.ones(25), np.zeros(N - 25)]
    cov = np.random.default_rng(1).normal(size=(N, 2))
    genes = []
    for g in range(n_genes):
        n_var = 0 if g == 5 else 5
        codes = rng.choice(["0/0", "0/1", "1/1"], size=(n_var, N), p=[0.85, 0.12, 0.03])
        if g == 3:
            codes[:] = "0/0"
            codes[0, :2] = "0/1"  # MAC < 5: no test matrix, still in the artifact
        if g in (1, 7, 9):
            codes[:, g : g + 3] = "./."  # these samples fail sample QC in this gene only
        df = pd.DataFrame(
            {"CHROM": "1", "POS": np.arange(n_var) * 10 + 100 * g, "REF": "A", "ALT": "G"}
        )
        df = pd.concat([df, pd.DataFrame(codes, columns=samples)], axis=1)
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=samples,
            gt_columns=samples,
            is_binary=True,
            missing_site_threshold=0.5,
            missing_sample_threshold=0.5,
            phenotype_vector=pheno,
            covariate_matrix=cov,
        )
        genes.append(
            {
                "GENE": f"G{g:02d}",
                "proband_count": 25,
                "control_count": N - 25,
                "n_qualifying_variants": n_var,
                "_genotype_matrix_builder": builder,
            }
        )
    return genes


def _arrays(path: str) -> dict[str, np.ndarray]:
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def _assert_same_artifact(a: str, b: str) -> None:
    xa, xb = _arrays(a + ".npz"), _arrays(b + ".npz")
    assert list(xa) == list(xb)
    for k in xa:
        assert xa[k].tobytes() == xb[k].tobytes(), k
    ga, gb = load_artifact(a).genes, load_artifact(b).genes
    assert {g: (s.variant_keys, s.n_carriers) for g, s in ga.items()} == {
        g: (s.variant_keys, s.n_carriers) for g, s in gb.items()
    }


@pytest.fixture(autouse=True)
def _parallel_for_small_panels(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("variantcentrifuge.association.engine.MIN_GENES_FOR_PARALLEL", 2)


@pytest.mark.unit
@pytest.mark.parametrize("workers", [1, 3])
def test_engine_pass_artifact_equals_parent_built(tmp_path, workers: int) -> None:
    """Bit-identical arrays and gene order vs building every gene in the parent."""
    ref_builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in _genes()}
    emit_stratum_artifact(str(tmp_path / "ref"), ref_builders, stratum_id="ref", **KW)

    genes = _genes()
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in genes}
    engine = AssociationEngine.from_names(
        ["logistic_burden", "skat_python"],
        AssociationConfig(association_workers=workers, skat_method="SKATO"),
    )
    engine.run_all(genes, gene_hook=_score_artifact_hook(builders))
    assert set(engine.gene_hook_results) == set(builders)
    emit_stratum_artifact(
        str(tmp_path / "new"),
        builders,
        records=engine.gene_hook_results,
        stratum_id="ref",
        **KW,
    )
    _assert_same_artifact(str(tmp_path / "ref"), str(tmp_path / "new"))
    assert "G03" in load_artifact(str(tmp_path / "new")).genes  # MAC<5 gene kept


@pytest.mark.unit
def test_each_gene_is_parsed_once_with_the_artifact(tmp_path) -> None:
    """The engine pass builds each gene once; the artifact reuses that build."""
    calls = multiprocessing.get_context("fork").Value("i", 0)
    real = _GenotypeMatrixBuilder.gene_arrays

    def counting(self):
        with calls.get_lock():
            calls.value += 1
        return real(self)

    genes = _genes()
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in genes}
    engine = AssociationEngine.from_names(
        ["logistic_burden"], AssociationConfig(association_workers=3)
    )
    with patch.object(_GenotypeMatrixBuilder, "gene_arrays", counting):
        engine.run_all(genes, gene_hook=_score_artifact_hook(builders))
        emit_stratum_artifact(
            str(tmp_path / "s"), builders, records=engine.gene_hook_results, stratum_id="s", **KW
        )
    assert calls.value == len(genes) - 1  # the gene without variants is never parsed


@pytest.mark.unit
def test_full_sample_null_fitted_once_in_parent() -> None:
    """Without per-gene sample masks the artifact's null model is fitted once per run."""
    calls = multiprocessing.get_context("fork").Value("i", 0)
    real = art_mod.fit_stratum_null

    def counting(*args, **kwargs):
        with calls.get_lock():
            calls.value += 1
        return real(*args, **kwargs)

    genes = [g for i, g in enumerate(_genes()) if i not in (1, 7, 9)]
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in genes}
    engine = AssociationEngine.from_names(
        ["logistic_burden"], AssociationConfig(association_workers=3)
    )
    with patch.object(art_mod, "fit_stratum_null", counting):
        engine.run_all(genes, gene_hook=_score_artifact_hook(builders))
    assert calls.value == 1


@pytest.mark.unit
def test_null_cache_is_bounded_without_changing_results(tmp_path, monkeypatch) -> None:
    """Evicting cached per-mask null models only costs refits; output is unchanged."""
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in _genes()}
    emit_stratum_artifact(str(tmp_path / "a"), builders, stratum_id="s", **KW)
    sizes: list[int] = []
    real = art_mod.null_for_mask

    def tracking(nulls, *args):
        out = real(nulls, *args)
        sizes.append(len(nulls))
        return out

    monkeypatch.setattr(art_mod, "_MAX_CACHED_NULLS", 1)
    monkeypatch.setattr(art_mod, "null_for_mask", tracking)
    emit_stratum_artifact(str(tmp_path / "b"), builders, stratum_id="s", **KW)
    assert sizes and max(sizes) == 1
    _assert_same_artifact(str(tmp_path / "a"), str(tmp_path / "b"))
