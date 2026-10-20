"""Quantitative outcome alignment for association analysis."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping

import numpy as np

logger = logging.getLogger(__name__)

# Tokens that mean "no outcome recorded"; anything else must parse as a finite number.
_MISSING_TOKENS = {"", "na", "nan", "."}


def build_quantitative_outcome(
    phenotype_data: Mapping[str, set[str]] | None,
    vcf_samples: list[str],
    value_column: str | None,
    remove_substring: str | None = None,
) -> tuple[list[str], np.ndarray]:
    """Return VCF-ordered samples with a numeric outcome and the aligned outcome vector.

    ``phenotype_data`` is the sample -> values map loaded by ``PhenotypeLoadingStage`` from
    ``--phenotype-value-column``. Samples without a value are excluded (and counted);
    non-numeric, non-finite or conflicting values abort, since fitting a different
    outcome than the one requested would be silently wrong. ``vcf_samples`` are already
    stripped of ``remove_substring``; phenotype IDs get the same removal, as in the binary
    phenotype-based case/control matching.
    """
    if not phenotype_data:
        raise ValueError(
            "Quantitative association needs the numeric outcome: provide --phenotype-file "
            "with --phenotype-sample-column and --phenotype-value-column."
        )
    if remove_substring and remove_substring.strip():
        merged: dict[str, set[str]] = {}
        for sample, raw in phenotype_data.items():
            merged.setdefault(sample.replace(remove_substring, ""), set()).update(raw)
        phenotype_data = merged
    samples: list[str] = []
    values: list[float] = []
    invalid: list[str] = []
    for sample in vcf_samples:
        present = {
            v for v in phenotype_data.get(sample, set()) if v.strip().lower() not in _MISSING_TOKENS
        }
        if not present:
            continue
        if len(present) > 1:
            raise ValueError(
                f"Sample '{sample}' has conflicting quantitative outcome values "
                f"{sorted(present)} in column '{value_column}'."
            )
        text = present.pop()
        try:
            value = float(text)
        except ValueError:
            value = math.nan
        if not math.isfinite(value):
            invalid.append(f"{sample}={text!r}")
            continue
        samples.append(sample)
        values.append(value)

    if invalid:
        raise ValueError(
            f"Quantitative outcome column '{value_column}' has {len(invalid)} non-numeric or "
            f"non-finite value(s), e.g. {', '.join(invalid[:5])}."
        )
    if not samples:
        raise ValueError(
            f"No VCF sample has a quantitative outcome in column '{value_column}'; "
            "check that phenotype sample IDs match the VCF."
        )
    n_missing = len(vcf_samples) - len(samples)
    if n_missing:
        logger.warning(
            "Association analysis: %d of %d VCF sample(s) have no '%s' outcome and are excluded",
            n_missing,
            len(vcf_samples),
            value_column,
        )
    return samples, np.asarray(values, dtype=float)
