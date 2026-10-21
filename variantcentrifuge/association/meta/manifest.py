# File: variantcentrifuge/association/meta/manifest.py
"""
Study-manifest reading and cross-stratum comparability validation.

The manifest is a small TSV: one row per stratum with columns ``stratum_id`` and
``artifact_path``, plus an optional ``mask_id`` column for several masks in one run
(one row per mask and stratum; each mask is validated and combined separately, #196).
Before any combination, strata must be comparable — the meta step refuses (hard
error, no silent fallback) when strata disagree on genome build, annotation
version, mask definition, or trait type, or when their (hashed) sample sets
overlap (overlap violates the additive-covariance independence assumption;
issue #112 forbids shared controls).

Sample-hash non-overlap does NOT by itself establish independence for related or
shared-control populations — that remains an explicit caller assumption.
"""

from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass

from variantcentrifuge.association.meta.artifact import ScoreCovarianceArtifact, load_artifact
from variantcentrifuge.association.meta.provenance import check_provenance

logger = logging.getLogger("variantcentrifuge")


class MetaCompatibilityError(ValueError):
    """Raised when strata are not comparable and must not be pooled."""


@dataclass
class StratumRef:
    """One manifest row: a stratum id and the path to its artifact."""

    stratum_id: str
    artifact_path: str
    mask_id: str | None = None


def read_study_manifest(path: str) -> list[StratumRef]:
    """Read a study manifest TSV (columns: stratum_id, artifact_path[, mask_id])."""
    refs: list[StratumRef] = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"stratum_id", "artifact_path"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise MetaCompatibilityError(
                f"Manifest '{path}' must have columns {sorted(required)}; found {reader.fieldnames}"
            )
        has_mask = "mask_id" in reader.fieldnames
        for row in reader:
            sid = (row.get("stratum_id") or "").strip()
            apath = (row.get("artifact_path") or "").strip()
            mask = (row.get("mask_id") or "").strip() if has_mask else None
            if not sid or not apath or mask == "":
                raise MetaCompatibilityError(
                    f"Manifest '{path}' has a row with empty fields: {row}"
                )
            # mask_id names a diagnostics subdirectory: keep it a plain file name
            if mask and (not re.fullmatch(r"[A-Za-z0-9._-]+", mask) or mask in (".", "..")):
                raise MetaCompatibilityError(
                    f"Manifest '{path}': mask_id '{mask}' may only use letters, digits, "
                    "'.', '_' and '-' (and must not be '.' or '..')."
                )
            # Relative paths are relative to the manifest, not the working directory.
            apath = os.path.join(os.path.dirname(os.path.abspath(path)), apath)
            refs.append(StratumRef(stratum_id=sid, artifact_path=apath, mask_id=mask))
    if not refs:
        raise MetaCompatibilityError(
            f"Manifest '{path}' lists no strata; meta-analysis needs at least 2."
        )
    for mask, group in group_by_mask(refs).items():
        what = f"mask '{mask}'" if mask is not None else "it"
        if len(group) < 2:
            raise MetaCompatibilityError(
                f"Manifest '{path}' lists {len(group)} stratum for {what}; "
                "meta-analysis needs at least 2."
            )
        ids = [r.stratum_id for r in group]
        dup = sorted({s for s in ids if ids.count(s) > 1})
        if dup:
            raise MetaCompatibilityError(f"Duplicate stratum_id {dup} for {what} in manifest.")
    return refs


def group_by_mask(refs: list[StratumRef]) -> dict[str | None, list[StratumRef]]:
    """Manifest rows per mask, in first-appearance order (one group without ``mask_id``)."""
    groups: dict[str | None, list[StratumRef]] = {}
    for ref in refs:
        groups.setdefault(ref.mask_id, []).append(ref)
    return groups


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


def _drop_mixed_robust(artifacts: list[ScoreCovarianceArtifact]) -> None:
    """Mixed --skat-robust and default strata: warn once, use the uncalibrated V (#223)."""
    robust = {
        a.stratum_id: any(g.v_robust is not None for g in a.genes.values()) for a in artifacts
    }
    if len(set(robust.values())) > 1:
        listing = ", ".join(f"{s}={'robust' if r else 'default'}" for s, r in robust.items())
        logger.warning(
            "Strata mix --skat-robust and default score artifacts (%s); meta-analysis uses "
            "the uncalibrated covariance V for all strata.",
            listing,
        )
        for art in artifacts:
            for gs in art.genes.values():
                gs.v_robust = None


def validate_comparability(
    artifacts: list[ScoreCovarianceArtifact], allow_annotation_mismatch: bool = False
) -> None:
    """Raise MetaCompatibilityError unless all strata are comparable.

    Checks: identical genome_build, annotation_version, mask_hash, trait_type,
    missing-genotype imputation, annotation provenance (a mismatch only warns with
    ``allow_annotation_mismatch``; unknown values warn), and disjoint
    sample_id_hashes across strata. Calibrated variances are dropped (one warning)
    unless every stratum has them.
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
    _drop_mixed_robust(artifacts)
    check_provenance(artifacts, allow_annotation_mismatch, MetaCompatibilityError)

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
    "group_by_mask",
    "load_manifest_artifacts",
    "read_study_manifest",
    "validate_comparability",
]
