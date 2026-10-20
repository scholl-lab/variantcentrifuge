"""Picklable genotype matrix builder and column alias utilities."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class _GenotypeMatrixBuilder:
    """Picklable lazy builder for per-gene genotype matrices (PERF-06).

    Instead of pre-building all matrices, this callable is stored in gene_data
    and invoked by the engine just before running tests. The matrix is discarded
    after test execution, keeping peak memory at O(1 gene).
    """

    gene_df: pd.DataFrame
    vcf_samples: list[str]
    gt_columns: list[str]
    is_binary: bool
    missing_site_threshold: float
    missing_sample_threshold: float
    phenotype_vector: np.ndarray | None
    covariate_matrix: np.ndarray | None
    score_column: str | None = None
    cadd_column: str | None = None
    revel_column: str | None = None
    effect_column: str | None = None

    def _aligned_annotation_payload(self, keep_variants_mask: np.ndarray) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if self.score_column is not None:
            payload["score_values"] = self.gene_df[self.score_column].to_numpy(dtype=object)[
                keep_variants_mask
            ]
            payload["variant_weight_column"] = self.score_column
        if self.cadd_column is not None:
            payload["cadd_scores"] = self.gene_df[self.cadd_column].to_numpy(dtype=object)[
                keep_variants_mask
            ]
        if self.revel_column is not None:
            payload["revel_scores"] = self.gene_df[self.revel_column].to_numpy(dtype=object)[
                keep_variants_mask
            ]
        if self.effect_column is not None:
            payload["variant_effects"] = self.gene_df[self.effect_column].to_numpy(dtype=object)[
                keep_variants_mask
            ]
        return payload

    def _empty_annotation_payload(self) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if self.score_column is not None:
            payload["score_values"] = np.asarray([], dtype=object)
            payload["variant_weight_column"] = self.score_column
        if self.cadd_column is not None:
            payload["cadd_scores"] = np.asarray([], dtype=object)
        if self.revel_column is not None:
            payload["revel_scores"] = np.asarray([], dtype=object)
        if self.effect_column is not None:
            payload["variant_effects"] = np.asarray([], dtype=object)
        return payload

    def __call__(self) -> dict[str, Any]:
        """Build genotype matrix and apply sample mask + MAC check.

        Returns dict with keys: genotype_matrix, variant_mafs,
        phenotype_vector, covariate_matrix, gt_warnings, mac_filtered.
        """
        from ..association.genotype_matrix import build_genotype_matrix

        if self.gene_df.empty:
            n_samples = len(self.vcf_samples)
            result: dict[str, Any] = {
                "genotype_matrix": np.zeros((n_samples, 0), dtype=float),
                "variant_mafs": np.zeros(0, dtype=float),
                "phenotype_vector": self.phenotype_vector,
                "covariate_matrix": self.covariate_matrix,
                "gt_warnings": [],
                "mac_filtered": False,
            }
            result.update(self._empty_annotation_payload())
            return result

        geno, mafs, sample_mask, gt_warnings, keep_variants_mask = cast(
            tuple[np.ndarray, np.ndarray, list[bool], list[str], np.ndarray],
            build_genotype_matrix(
                self.gene_df,
                self.vcf_samples,
                self.gt_columns,
                is_binary=self.is_binary,
                missing_site_threshold=self.missing_site_threshold,
                missing_sample_threshold=self.missing_sample_threshold,
                phenotype_vector=self.phenotype_vector,
                return_keep_mask=True,
            ),
        )
        annotation_payload = self._aligned_annotation_payload(keep_variants_mask)
        # Apply sample mask to phenotype and covariates
        pv = self.phenotype_vector
        cm = self.covariate_matrix
        if not all(sample_mask):
            mask_arr = np.array(sample_mask, dtype=bool)
            pv = pv[mask_arr] if pv is not None else None
            cm = cm[mask_arr] if cm is not None else None
            geno = geno[mask_arr]

        # Per-gene MAC check: skip regression if < 5 minor allele copies
        mac_filtered = False
        total_mac = int(geno.sum()) if geno.size > 0 else 0
        if total_mac < 5:
            geno = np.zeros((geno.shape[0], 0), dtype=float)
            mafs = np.zeros(0, dtype=float)
            annotation_payload = self._empty_annotation_payload()
            mac_filtered = True

        result = {
            "genotype_matrix": geno,
            "variant_mafs": mafs,
            "phenotype_vector": pv,
            "covariate_matrix": cm,
            "gt_warnings": gt_warnings,
            "mac_filtered": mac_filtered,
        }
        result.update(annotation_payload)
        return result


def _find_cadd_column(df: pd.DataFrame) -> str | None:
    return next((c for c in df.columns if c.lower() in ("dbnsfp_cadd_phred", "cadd_phred")), None)


def _find_revel_column(df: pd.DataFrame) -> str | None:
    return next((c for c in df.columns if c.lower() in ("dbnsfp_revel_score", "revel_score")), None)


def _find_effect_column(df: pd.DataFrame) -> str | None:
    return next((c for c in df.columns if c.upper() in ("EFFECT", "ANN_0__EFFECT")), None)


# Standard column aliases: sanitized name -> canonical short name.
# When ANN[0].GENE is sanitized to ANN_0__GENE, downstream stages still expect "GENE".
_STANDARD_COLUMN_ALIASES = {
    "ANN_0__GENE": "GENE",
    "ANN____GENE": "GENE",  # ANN[*].GENE
    "EFF____GENE": "GENE",  # EFF[*].GENE
    "ANN_0__EFFECT": "EFFECT",
    "ANN_0__IMPACT": "IMPACT",
    "ANN_0__HGVS_C": "HGVS_C",
    "ANN_0__HGVS_P": "HGVS_P",
}


def _create_standard_column_aliases(df: pd.DataFrame) -> None:
    """Create standard short-name aliases for sanitized ANN columns.

    After column sanitization (ANN[0].GENE -> ANN_0__GENE), many downstream
    stages expect plain column names like GENE, EFFECT, IMPACT. This function
    creates those aliases in-place if the sanitized source column exists and
    the short name doesn't already exist.
    """
    for sanitized_name, short_name in _STANDARD_COLUMN_ALIASES.items():
        if short_name not in df.columns and sanitized_name in df.columns:
            df[short_name] = df[sanitized_name]
            logger.debug(f"Created column alias: {sanitized_name} -> {short_name}")
