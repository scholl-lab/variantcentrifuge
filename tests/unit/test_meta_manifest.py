"""Manifest reading + comparability/safety validation (acceptance #4, #5)."""

from __future__ import annotations

import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.manifest import (
    MetaCompatibilityError,
    read_study_manifest,
    validate_comparability,
)

KEYS = {"GENEA": [f"1:{100 + j}:A:G" for j in range(6)]}


def _pair(**overrides_b):
    a = make_stratum("A", seed=11, gene_keys=KEYS)
    b = make_stratum("B", seed=22, gene_keys=KEYS, **overrides_b)
    return a, b


def test_comparable_strata_pass():
    a, b = _pair()
    validate_comparability([a, b])  # no raise


@pytest.mark.parametrize(
    "override",
    [
        {"genome_build": "GRCh37"},
        {"annotation_version": "vep-110"},
        {"mask_hash": "different-hash"},
    ],
)
def test_mismatched_metadata_refused(override):
    """Acceptance #4: build/annotation/mask mismatch → explicit error, no pooling."""
    a, b = _pair(**override)
    with pytest.raises(MetaCompatibilityError, match="disagree on"):
        validate_comparability([a, b])


def test_sample_overlap_refused():
    """Acceptance #5: shared samples across strata → refuse (independence violated)."""
    a = make_stratum("A", seed=1, gene_keys=KEYS, sample_prefix="shared")
    b = make_stratum("B", seed=2, gene_keys=KEYS, sample_prefix="shared")
    with pytest.raises(MetaCompatibilityError, match="Sample overlap"):
        validate_comparability([a, b])


def test_mixed_trait_type_refused():
    a, b = _pair()
    object.__setattr__(b, "trait_type", "quantitative")
    with pytest.raises(MetaCompatibilityError, match="trait_type"):
        validate_comparability([a, b])


def test_manifest_requires_two_strata(tmp_path):
    m = tmp_path / "manifest.tsv"
    m.write_text("stratum_id\tartifact_path\nA\t/tmp/a.npz\n")
    with pytest.raises(MetaCompatibilityError, match="at least 2"):
        read_study_manifest(str(m))


def test_manifest_missing_columns(tmp_path):
    m = tmp_path / "bad.tsv"
    m.write_text("id\tpath\nA\tx\nB\ty\n")
    with pytest.raises(MetaCompatibilityError, match="must have columns"):
        read_study_manifest(str(m))
