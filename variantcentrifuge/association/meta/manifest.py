# File: variantcentrifuge/association/meta/manifest.py
"""
Study-manifest reading and cross-stratum comparability validation.

The manifest is a small TSV: one row per stratum with columns ``stratum_id`` and
``artifact_path``. Before any combination, strata must be comparable — the meta
step refuses (hard error, no silent fallback) when strata disagree on genome
build, annotation version, mask definition, or trait type, or when their
(hashed) sample sets overlap (overlap violates the additive-covariance
independence assumption; issue #112 forbids shared controls).

Sample-hash non-overlap does NOT by itself establish independence for related or
shared-control populations — that remains an explicit caller assumption.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass

from variantcentrifuge.association.meta.artifact import ScoreCovarianceArtifact, load_artifact

logger = logging.getLogger("variantcentrifuge")


class MetaCompatibilityError(ValueError):
    """Raised when strata are not comparable and must not be pooled."""


@dataclass
class StratumRef:
    """One manifest row: a stratum id and the path to its artifact."""

    stratum_id: str
    artifact_path: str


def read_study_manifest(path: str) -> list[StratumRef]:
    """Read a study manifest TSV (columns: stratum_id, artifact_path)."""
    refs: list[StratumRef] = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"stratum_id", "artifact_path"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise MetaCompatibilityError(
                f"Manifest '{path}' must have columns {sorted(required)}; found {reader.fieldnames}"
            )
        for row in reader:
            sid = (row.get("stratum_id") or "").strip()
            apath = (row.get("artifact_path") or "").strip()
            if not sid or not apath:
                raise MetaCompatibilityError(
                    f"Manifest '{path}' has a row with empty fields: {row}"
                )
            refs.append(StratumRef(stratum_id=sid, artifact_path=apath))
    if len(refs) < 2:
        raise MetaCompatibilityError(
            f"Manifest '{path}' lists {len(refs)} stratum; meta-analysis needs at least 2."
        )
    seen: set[str] = set()
    for ref in refs:
        if ref.stratum_id in seen:
            raise MetaCompatibilityError(f"Duplicate stratum_id '{ref.stratum_id}' in manifest.")
        seen.add(ref.stratum_id)
    return refs


def load_manifest_artifacts(refs: list[StratumRef]) -> list[ScoreCovarianceArtifact]:
    """Load every artifact referenced by the manifest (clear error on missing file)."""
    artifacts: list[ScoreCovarianceArtifact] = []
    for ref in refs:
        try:
            artifacts.append(load_artifact(ref.artifact_path))
        except FileNotFoundError as exc:
            raise MetaCompatibilityError(
                f"Artifact for stratum '{ref.stratum_id}' not found: {ref.artifact_path}"
            ) from exc
    return artifacts


def validate_comparability(artifacts: list[ScoreCovarianceArtifact]) -> None:
    """Raise MetaCompatibilityError unless all strata are comparable.

    Checks: identical genome_build, annotation_version, mask_hash, trait_type,
    and disjoint sample_id_hashes across strata.
    """
    if len(artifacts) < 2:
        raise MetaCompatibilityError("Need at least 2 strata to validate comparability.")

    def _uniform(attr: str) -> None:
        values = {getattr(a, attr) for a in artifacts}
        if len(values) > 1:
            raise MetaCompatibilityError(
                f"Strata disagree on '{attr}': {sorted(map(str, values))}. "
                "Refusing to pool non-comparable strata (no silent fallback)."
            )

    for attr in ("genome_build", "annotation_version", "mask_hash", "trait_type"):
        _uniform(attr)

    # Sample-hash overlap check (independence of covariances).
    seen: dict[str, str] = {}
    for art in artifacts:
        for h in art.sample_id_hashes:
            if h in seen and seen[h] != art.stratum_id:
                raise MetaCompatibilityError(
                    f"Sample overlap between strata '{seen[h]}' and '{art.stratum_id}'. "
                    "Shared samples/controls violate the additive-covariance assumption."
                )
            seen[h] = art.stratum_id


__all__ = [
    "MetaCompatibilityError",
    "StratumRef",
    "load_manifest_artifacts",
    "read_study_manifest",
    "validate_comparability",
]
