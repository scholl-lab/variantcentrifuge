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
import os
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
            # Relative paths are relative to the manifest, not the working directory.
            apath = os.path.join(os.path.dirname(os.path.abspath(path)), apath)
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
    """Load every artifact referenced by the manifest (clear error on missing file).

    The manifest's ``stratum_id`` names the stratum in all outputs.
    """
    artifacts: list[ScoreCovarianceArtifact] = []
    for ref in refs:
        try:
            art = load_artifact(ref.artifact_path)
        except FileNotFoundError as exc:
            raise MetaCompatibilityError(
                f"Artifact for stratum '{ref.stratum_id}' not found: {ref.artifact_path}"
            ) from exc
        art.stratum_id = ref.stratum_id
        artifacts.append(art)
    return artifacts


def _check_missing_policy(artifacts: list[ScoreCovarianceArtifact]) -> None:
    """Refuse mixed missing-genotype imputation; warn on differing QC thresholds (#200)."""
    by_stratum = {a.stratum_id: a.missing_policy for a in artifacts}
    if len({p["imputation"] for p in by_stratum.values()}) > 1:
        listing = ", ".join(f"{s}={p['imputation']}" for s, p in by_stratum.items())
        raise MetaCompatibilityError(
            f"Strata use different missing-genotype imputation ({listing}). Re-run the "
            "strata with the same --missing-as-ref setting before pooling."
        )
    for key in ("site_threshold", "sample_threshold"):
        if len({p.get(key) for p in by_stratum.values()}) > 1:
            listing = ", ".join(f"{s}={p.get(key)}" for s, p in by_stratum.items())
            logger.warning("Strata use different missingness %s (%s).", key, listing)


def validate_comparability(artifacts: list[ScoreCovarianceArtifact]) -> None:
    """Raise MetaCompatibilityError unless all strata are comparable.

    Checks: identical genome_build, annotation_version, mask_hash, trait_type,
    missing-genotype imputation, and disjoint sample_id_hashes across strata.
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
    _check_missing_policy(artifacts)

    # Sample-hash overlap check (independence of covariances).
    seen: dict[str, str] = {}
    for art in artifacts:
        if not art.sample_id_hashes:
            logger.warning(
                "Stratum '%s' has no sample hashes; sample overlap cannot be checked.",
                art.stratum_id,
            )
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
