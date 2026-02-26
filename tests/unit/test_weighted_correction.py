# File: tests/unit/test_weighted_correction.py
"""
Unit tests for weighted Benjamini-Hochberg FDR correction and gene weight loading.

Tests cover:
- load_gene_weights(): valid file, missing column errors, invalid weight errors
- apply_weighted_correction(): mathematical properties, edge cases, warnings
- Backward compatibility: apply_correction() unchanged, uniform weights reproduce it
"""

from __future__ import annotations

import logging
import os
import tempfile

import numpy as np
import pytest

from variantcentrifuge.association.correction import (
    apply_correction,
    apply_weighted_correction,
    load_gene_weights,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_weight_file(content: str) -> str:
    """Write content to a temp TSV file and return path."""
    fd, path = tempfile.mkstemp(suffix=".tsv", prefix="weights_")
    with os.fdopen(fd, "w") as fh:
        fh.write(content)
    return path


# ---------------------------------------------------------------------------
# Tests for load_gene_weights()
# ---------------------------------------------------------------------------


class TestLoadGeneWeights:
    def test_valid_file_two_columns(self):
        """Standard two-column TSV: gene, weight."""
        content = "gene\tweight\nGENEA\t1.5\nGENEB\t0.5\nGENEC\t2.0\n"
        path = _write_weight_file(content)
        try:
            weights = load_gene_weights(path)
            assert weights == {"GENEA": 1.5, "GENEB": 0.5, "GENEC": 2.0}
        finally:
            os.unlink(path)

    def test_valid_file_extra_columns(self):
        """Extra columns are allowed; only the named weight column is used."""
        content = "gene\tpLI\tweight\tnotes\nTP53\t0.95\t3.0\toncogene\nBRCA1\t0.80\t2.0\t\n"
        path = _write_weight_file(content)
        try:
            weights = load_gene_weights(path, weight_column="weight")
            assert weights["TP53"] == pytest.approx(3.0)
            assert weights["BRCA1"] == pytest.approx(2.0)
        finally:
            os.unlink(path)

    def test_custom_weight_column(self):
        """Non-default weight column name is resolved correctly."""
        content = "gene\tpLI_score\nTP53\t0.99\nBRCA1\t0.50\n"
        path = _write_weight_file(content)
        try:
            weights = load_gene_weights(path, weight_column="pLI_score")
            assert weights["TP53"] == pytest.approx(0.99)
        finally:
            os.unlink(path)

    def test_missing_weight_column_raises(self):
        """Missing weight column raises ValueError with helpful message."""
        content = "gene\tpLI\nTP53\t0.99\n"
        path = _write_weight_file(content)
        try:
            with pytest.raises(ValueError, match=r"weight.*not found"):
                load_gene_weights(path, weight_column="weight")
        finally:
            os.unlink(path)

    def test_zero_weight_raises(self):
        """Zero weight raises ValueError naming the offending gene."""
        content = "gene\tweight\nTP53\t1.0\nBRCA1\t0.0\n"
        path = _write_weight_file(content)
        try:
            with pytest.raises(ValueError, match="zero or negative"):
                load_gene_weights(path)
        finally:
            os.unlink(path)

    def test_negative_weight_raises(self):
        """Negative weight raises ValueError naming the offending gene."""
        content = "gene\tweight\nTP53\t1.0\nBRCA1\t-0.5\n"
        path = _write_weight_file(content)
        try:
            with pytest.raises(ValueError, match="zero or negative"):
                load_gene_weights(path)
        finally:
            os.unlink(path)

    def test_multiple_bad_weights_listed(self):
        """Multiple bad genes are listed in the error message."""
        lines = ["gene\tweight"]
        for i in range(12):
            lines.append(f"GENE{i}\t{-i if i > 0 else 0.0}")
        content = "\n".join(lines) + "\n"
        path = _write_weight_file(content)
        try:
            with pytest.raises(ValueError, match=r"and \d+ more"):
                load_gene_weights(path)
        finally:
            os.unlink(path)

    def test_logs_count(self, caplog):
        """load_gene_weights logs the number of loaded weights."""
        content = "gene\tweight\nGENEA\t1.0\nGENEB\t2.0\n"
        path = _write_weight_file(content)
        try:
            with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
                load_gene_weights(path)
            assert any("2 gene weights" in r.message for r in caplog.records)
        finally:
            os.unlink(path)


# ---------------------------------------------------------------------------
# Tests for apply_weighted_correction()
# ---------------------------------------------------------------------------


class TestApplyWeightedCorrection:
    """Mathematical and edge-case tests for apply_weighted_correction()."""

    @pytest.mark.parametrize("method", ["fdr", "bonferroni"])
    def test_uniform_weights_match_apply_correction(self, method):
        """Uniform weights (all 1.0) should give same result as apply_correction()."""
        pvals = [0.001, 0.01, 0.05, 0.1, 0.5]
        genes = ["G1", "G2", "G3", "G4", "G5"]
        weight_map = dict.fromkeys(genes, 1.0)

        corrected_unweighted = apply_correction(pvals, method=method)
        corrected_weighted, norm_w = apply_weighted_correction(
            pvals, genes, weight_map, method=method
        )

        np.testing.assert_allclose(corrected_weighted, corrected_unweighted, rtol=1e-10)
        # Normalized weights should all be 1.0 with uniform inputs
        np.testing.assert_allclose(norm_w, np.ones(len(pvals)), rtol=1e-10)

    def test_non_uniform_weights_differ_from_unweighted(self):
        """Non-uniform weights should produce different results from unweighted BH."""
        pvals = [0.001, 0.01, 0.05, 0.1, 0.5]
        genes = ["G1", "G2", "G3", "G4", "G5"]
        # G1 has very high weight (biologically important), G5 very low
        weight_map = {"G1": 5.0, "G2": 1.0, "G3": 1.0, "G4": 1.0, "G5": 0.2}

        corrected_unweighted = apply_correction(pvals, method="fdr")
        corrected_weighted, _norm_w = apply_weighted_correction(
            pvals, genes, weight_map, method="fdr"
        )

        # Results must differ (high weight on G1 should help it)
        assert not np.allclose(corrected_weighted, corrected_unweighted)

    def test_missing_genes_get_weight_one(self):
        """Genes not in weight_map receive weight=1.0."""
        pvals = [0.01, 0.05]
        genes = ["KNOWN", "UNKNOWN"]
        weight_map = {"KNOWN": 1.0}  # UNKNOWN not in map

        _, norm_w = apply_weighted_correction(pvals, genes, weight_map)
        # Both should be 1.0 because KNOWN=1.0 and UNKNOWN defaults to 1.0 → uniform
        np.testing.assert_allclose(norm_w, np.ones(2), rtol=1e-10)

    def test_single_gene_returned_unchanged(self):
        """Single gene input is returned unchanged (weighting has no effect)."""
        pvals = [0.03]
        genes = ["G1"]
        weight_map = {"G1": 5.0}

        corrected, norm_w = apply_weighted_correction(pvals, genes, weight_map)
        assert float(corrected[0]) == pytest.approx(0.03)
        assert float(norm_w[0]) == pytest.approx(1.0)

    def test_single_gene_logs_info(self, caplog):
        """Single gene path logs an info message."""
        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            apply_weighted_correction([0.05], ["G1"], {"G1": 2.0})
        assert any("Single gene" in r.message for r in caplog.records)

    def test_empty_pvals_returns_empty(self):
        """Empty input returns empty arrays."""
        corrected, norm_w = apply_weighted_correction([], [], {})
        assert len(corrected) == 0
        assert len(norm_w) == 0

    def test_coverage_warning_50_percent(self, caplog):
        """Warning emitted when >50% of genes are missing from weight map."""
        genes = [f"G{i}" for i in range(10)]
        pvals = [0.1] * 10
        # Only 4 genes in map → 6/10 = 60% missing → should warn at >50%
        weight_map = {f"G{i}": 1.0 for i in range(4)}

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            apply_weighted_correction(pvals, genes, weight_map)

        warning_messages = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("50%" in m for m in warning_messages)

    def test_coverage_warning_80_percent(self, caplog):
        """Stronger warning emitted when >80% of genes are missing from weight map."""
        genes = [f"G{i}" for i in range(10)]
        pvals = [0.1] * 10
        # Only 1 gene in map → 9/10 = 90% missing → should trigger >80% strong warning
        weight_map = {"G0": 1.0}

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            apply_weighted_correction(pvals, genes, weight_map)

        warning_messages = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("80%" in m or "STRONG" in m for m in warning_messages)

    def test_extreme_weight_ratio_warning(self, caplog):
        """Warning emitted when max/min normalized weight ratio exceeds 100."""
        genes = [f"G{i}" for i in range(5)]
        pvals = [0.1] * 5
        # G0 has very high weight vs others → extreme ratio after normalization
        weight_map = {"G0": 500.0, "G1": 1.0, "G2": 1.0, "G3": 1.0, "G4": 1.0}

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            apply_weighted_correction(pvals, genes, weight_map)

        warning_messages = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("ratio" in m.lower() or "extreme" in m.lower() for m in warning_messages)

    def test_weights_renormalized_to_mean_one(self):
        """Normalized weights must have mean=1.0."""
        pvals = [0.01, 0.05, 0.1, 0.2, 0.5]
        genes = ["G1", "G2", "G3", "G4", "G5"]
        weight_map = {"G1": 3.0, "G2": 1.0, "G3": 2.0, "G4": 0.5, "G5": 1.5}

        _, norm_w = apply_weighted_correction(pvals, genes, weight_map)
        assert float(norm_w.mean()) == pytest.approx(1.0, rel=1e-9)

    def test_corrected_pvals_in_unit_interval(self):
        """Corrected p-values must be in [0, 1]."""
        pvals = [0.001, 0.002, 0.005, 0.01]
        genes = ["G1", "G2", "G3", "G4"]
        # Extreme weights to stress-test clipping
        weight_map = {"G1": 100.0, "G2": 0.01, "G3": 1.0, "G4": 1.0}

        corrected, _ = apply_weighted_correction(pvals, genes, weight_map)
        assert float(corrected.min()) >= 0.0
        assert float(corrected.max()) <= 1.0

    def test_effective_number_of_tests_logged(self, caplog):
        """Effective number of tests is logged at INFO level."""
        genes = ["G1", "G2", "G3"]
        pvals = [0.01, 0.05, 0.1]
        weight_map = {"G1": 2.0, "G2": 1.0, "G3": 0.5}

        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            apply_weighted_correction(pvals, genes, weight_map)

        assert any("effective" in r.message.lower() for r in caplog.records)

    @pytest.mark.parametrize("method", ["fdr", "bonferroni"])
    def test_bonferroni_method_path(self, method):
        """Bonferroni path runs without error and returns valid p-values."""
        pvals = [0.001, 0.01, 0.05]
        genes = ["G1", "G2", "G3"]
        weight_map = {"G1": 2.0, "G2": 1.0, "G3": 0.5}

        corrected, norm_w = apply_weighted_correction(pvals, genes, weight_map, method=method)
        assert len(corrected) == 3
        assert len(norm_w) == 3
        assert all(0.0 <= v <= 1.0 for v in corrected)

    def test_return_types_are_numpy_arrays(self):
        """Return types should be numpy arrays."""
        pvals = [0.01, 0.05]
        genes = ["G1", "G2"]
        weight_map = {"G1": 1.5, "G2": 0.5}

        corrected, norm_w = apply_weighted_correction(pvals, genes, weight_map)
        assert isinstance(corrected, np.ndarray)
        assert isinstance(norm_w, np.ndarray)


# ---------------------------------------------------------------------------
# Backward compatibility tests
# ---------------------------------------------------------------------------


class TestApplyCorrectionBackwardCompat:
    """Verify apply_correction() signature and output are unchanged."""

    def test_signature_unchanged(self):
        """apply_correction() still accepts (pvals, method='fdr') and returns ndarray."""
        pvals = [0.01, 0.05, 0.1]
        result = apply_correction(pvals)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_empty_input(self):
        """Empty input returns empty array."""
        result = apply_correction([])
        assert isinstance(result, np.ndarray)
        assert len(result) == 0

    def test_fdr_method(self):
        """FDR BH correction produces monotone q-values."""
        pvals = sorted([0.001, 0.01, 0.05, 0.1, 0.5])
        result = apply_correction(pvals, method="fdr")
        # BH corrected q-values should all be >= raw p-values
        assert all(q >= p for q, p in zip(result, pvals, strict=True))

    def test_bonferroni_method(self):
        """Bonferroni correction multiplies by n (clipped to 1)."""
        pvals = [0.01, 0.05]
        result = apply_correction(pvals, method="bonferroni")
        # Bonferroni: p * n = p * 2
        assert float(result[0]) == pytest.approx(min(0.01 * 2, 1.0))
        assert float(result[1]) == pytest.approx(min(0.05 * 2, 1.0))


# ---------------------------------------------------------------------------
# write_fdr_weight_diagnostics tests
# ---------------------------------------------------------------------------


class TestWriteFdrWeightDiagnostics:
    """Tests for write_fdr_weight_diagnostics()."""

    def test_creates_tsv_file(self, tmp_path):
        """Diagnostics TSV is created in the specified directory."""
        from variantcentrifuge.association.diagnostics import write_fdr_weight_diagnostics

        genes = ["GENEA", "GENEB", "GENEC"]
        raw_weights = {"GENEA": 2.0, "GENEB": 1.0}  # GENEC missing
        norm_w = np.array([1.5, 0.75, 0.75])
        unw_q = np.array([0.03, 0.06, 0.15])
        wt_q = np.array([0.02, 0.08, 0.18])

        write_fdr_weight_diagnostics(
            genes=genes,
            raw_weights=raw_weights,
            normalized_weights=norm_w,
            unweighted_corrected=unw_q,
            weighted_corrected=wt_q,
            fdr_threshold=0.05,
            diagnostics_dir=tmp_path,
        )

        out_file = tmp_path / "fdr_weight_diagnostics.tsv"
        assert out_file.exists()

    def test_tsv_columns(self, tmp_path):
        """Output TSV has correct columns."""
        import pandas as pd

        from variantcentrifuge.association.diagnostics import write_fdr_weight_diagnostics

        genes = ["G1", "G2"]
        write_fdr_weight_diagnostics(
            genes=genes,
            raw_weights={"G1": 1.0, "G2": 2.0},
            normalized_weights=np.array([0.667, 1.333]),
            unweighted_corrected=np.array([0.1, 0.2]),
            weighted_corrected=np.array([0.08, 0.25]),
            fdr_threshold=0.05,
            diagnostics_dir=tmp_path,
        )

        df = pd.read_csv(tmp_path / "fdr_weight_diagnostics.tsv", sep="\t")
        expected_cols = {
            "gene",
            "raw_weight",
            "normalized_weight",
            "unweighted_q",
            "weighted_q",
            "significance_change",
        }
        assert expected_cols <= set(df.columns)

    def test_significance_gained(self, tmp_path):
        """Gene that crosses threshold gains significance = 'gained'."""
        import pandas as pd

        from variantcentrifuge.association.diagnostics import write_fdr_weight_diagnostics

        genes = ["G1"]
        write_fdr_weight_diagnostics(
            genes=genes,
            raw_weights={"G1": 3.0},
            normalized_weights=np.array([1.0]),
            unweighted_corrected=np.array([0.06]),  # was NOT significant
            weighted_corrected=np.array([0.04]),  # IS significant with weight
            fdr_threshold=0.05,
            diagnostics_dir=tmp_path,
        )

        df = pd.read_csv(tmp_path / "fdr_weight_diagnostics.tsv", sep="\t")
        assert df.loc[df["gene"] == "G1", "significance_change"].values[0] == "gained"

    def test_significance_lost(self, tmp_path):
        """Gene that drops below threshold loses significance = 'lost'."""
        import pandas as pd

        from variantcentrifuge.association.diagnostics import write_fdr_weight_diagnostics

        genes = ["G1"]
        write_fdr_weight_diagnostics(
            genes=genes,
            raw_weights={},
            normalized_weights=np.array([1.0]),
            unweighted_corrected=np.array([0.04]),  # WAS significant
            weighted_corrected=np.array([0.06]),  # NOT significant with weight
            fdr_threshold=0.05,
            diagnostics_dir=tmp_path,
        )

        df = pd.read_csv(tmp_path / "fdr_weight_diagnostics.tsv", sep="\t")
        assert df.loc[df["gene"] == "G1", "significance_change"].values[0] == "lost"

    def test_no_significance_change(self, tmp_path):
        """Gene unchanged in significance has empty significance_change."""
        import pandas as pd

        from variantcentrifuge.association.diagnostics import write_fdr_weight_diagnostics

        genes = ["G1"]
        write_fdr_weight_diagnostics(
            genes=genes,
            raw_weights={"G1": 1.0},
            normalized_weights=np.array([1.0]),
            unweighted_corrected=np.array([0.08]),
            weighted_corrected=np.array([0.09]),
            fdr_threshold=0.05,
            diagnostics_dir=tmp_path,
        )

        df = pd.read_csv(tmp_path / "fdr_weight_diagnostics.tsv", sep="\t")
        val = df.loc[df["gene"] == "G1", "significance_change"].values[0]
        # Empty string or NaN (from fillna) is acceptable
        assert val == "" or (isinstance(val, float) and np.isnan(val))

    def test_logs_impact_summary(self, tmp_path, caplog):
        """Summary of gained/lost significance is logged."""
        from variantcentrifuge.association.diagnostics import write_fdr_weight_diagnostics

        genes = ["G1", "G2"]
        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            write_fdr_weight_diagnostics(
                genes=genes,
                raw_weights={"G1": 3.0},
                normalized_weights=np.array([1.5, 0.5]),
                unweighted_corrected=np.array([0.06, 0.04]),
                weighted_corrected=np.array([0.04, 0.06]),
                fdr_threshold=0.05,
                diagnostics_dir=tmp_path,
            )

        messages = [r.message for r in caplog.records]
        assert any("gained" in m or "lost" in m for m in messages)
