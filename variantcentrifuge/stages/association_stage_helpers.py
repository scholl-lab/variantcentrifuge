"""Helpers for the association analysis stage: artifact, cohort, covariates, COAST, warnings."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


def _emit_score_artifact(config: dict, path: str, gene_burden_data: list[dict]) -> None:
    """Write the meta-analysis score artifact from this run's per-gene matrix builders."""
    import hashlib
    import json

    from ..association.meta.artifact import emit_stratum_artifact
    from ..links import normalize_reference

    reference = str(config.get("reference") or "")
    # Strata are pooled only when they selected variants the same way.
    selection = {k: config.get(k) for k in ("filters", "final_filter", "bcftools_prefilter")}
    mask_hash = hashlib.sha256(json.dumps(selection, sort_keys=True).encode()).hexdigest()
    emit_stratum_artifact(
        path,
        {
            d["GENE"]: d["_genotype_matrix_builder"]
            for d in gene_burden_data
            if "_genotype_matrix_builder" in d
        },
        stratum_id=Path(path).name.split(".")[0],
        genome_build=normalize_reference(reference) or reference,
        annotation_version=reference,
        mask_id=mask_hash[:12],
        mask_hash=mask_hash,
    )


def _cohort_size_ok(
    n_cases_total: int, n_controls_total: int, max_ratio: float | None = 20.0
) -> bool:
    """Log binary cohort size/balance warnings; False when too few cases to analyse.

    ``max_ratio=None`` skips the imbalance warning (diagnostics emit it instead).
    """
    if n_cases_total < 10:
        logger.error(
            f"Association analysis: only {n_cases_total} case(s) found — "
            "fewer than 10 cases produces invalid results. Aborting."
        )
        return False
    if n_cases_total < 50:
        logger.warning(
            f"Association analysis: only {n_cases_total} case(s) — "
            "fewer than 50 cases provides no practical power for regression tests"
        )
    elif n_cases_total < 200:
        logger.warning(
            f"Association analysis: {n_cases_total} case(s) — "
            "fewer than 200 cases is underpowered for SKAT; interpret results with caution"
        )
    from ..association.diagnostics import imbalance_message

    msg = imbalance_message(n_cases_total, n_controls_total, max_ratio) if max_ratio else None
    if msg:
        logger.warning(msg)
    return True


def _load_covariates_and_pcs(
    assoc_config: Any,
    assigned_samples: list[str],
    remove_sample_substring: str | None,
    raw_categorical_out: dict[str, Any] | None = None,
) -> tuple[Any, list[str]]:
    """Load covariates and merge PCA components; (None, []) when nothing is configured."""
    covariate_matrix = None
    covariate_col_names: list[str] = []
    if assoc_config.covariate_file and assigned_samples:
        from ..association.covariates import load_covariates

        covariate_matrix, covariate_col_names = load_covariates(
            assoc_config.covariate_file,
            assigned_samples,
            assoc_config.covariate_columns,
            assoc_config.categorical_covariates,
            raw_categorical_out,
        )
        logger.info(
            f"Association analysis: loaded {covariate_matrix.shape[1]} "
            f"covariate(s): {covariate_col_names}"
        )

    pca_file = assoc_config.pca_file
    if pca_file and assigned_samples:
        from ..association.pca import load_pca_file, merge_pca_covariates

        n_pcs = assoc_config.pca_components
        pca_matrix, pca_col_names = load_pca_file(
            pca_file,
            assigned_samples,
            n_components=n_pcs,
            remove_sample_substring=remove_sample_substring,
        )
        covariate_matrix, covariate_col_names = merge_pca_covariates(
            pca_matrix,
            pca_col_names,
            covariate_matrix,
            covariate_col_names,
        )
        logger.info(
            f"PCA: merged {len(pca_col_names)} PCs with existing covariates "
            f"-> {len(covariate_col_names)} total columns"
        )
    return covariate_matrix, covariate_col_names


def _coast_missing_columns(df: pd.DataFrame, assoc_config: Any) -> list[str]:
    """Return the required COAST annotation columns absent from *df*."""
    _coast_base_required = [
        ("ANN_0__EFFECT", "ANN[0].EFFECT"),
        ("ANN_0__IMPACT", "ANN[0].IMPACT"),
    ]
    _coast_model = getattr(assoc_config, "coast_classification", None) or ""
    _is_sift_model = (
        not _coast_model
        or _coast_model.endswith("sift_polyphen")
        or _coast_model == "sift_polyphen"
    )
    if _is_sift_model:
        _coast_base_required.extend(
            [
                ("dbNSFP_SIFT_pred", "dbNSFP_SIFT_pred"),
                ("dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_Polyphen2_HDIV_pred"),
            ]
        )
    elif _coast_model.endswith("cadd") or _coast_model == "cadd":
        _coast_base_required.extend(
            [
                ("dbNSFP_CADD_phred", "dbNSFP_CADD_phred"),
            ]
        )
    missing = [
        san for san, alt in _coast_base_required if san not in df.columns and alt not in df.columns
    ]
    return missing


def _per_sample_gt_frame(df: pd.DataFrame, unpacked_samples: list[str]) -> pd.DataFrame:
    """Build a per-sample GT column copy of *df* from the packed GT column."""
    from .sample_column_helpers import create_sample_columns_from_gt_intelligent

    out = create_sample_columns_from_gt_intelligent(df.copy(), vcf_samples=unpacked_samples)
    gt_series_str = out["GT"].fillna("").astype(str)
    for s_col in unpacked_samples:
        if s_col not in out.columns:
            out[s_col] = gt_series_str.str.extract(
                rf"(?:^|;)\s*{re.escape(s_col)}\(([^:)]+)", expand=False
            ).fillna("./.")
        listed = gt_series_str.str.contains(rf"(?:^|;)\s*{re.escape(s_col)}\(", regex=True)
        fill_ref = (out[s_col] == "./.") & ~listed
        out.loc[fill_ref, s_col] = "0/0"
    return out


def _merge_carrier_warnings(
    results_df: pd.DataFrame,
    gene_data_by_gene: dict[str, Any],
    assoc_config: Any,
    quantitative: bool,
) -> None:
    """Append per-gene case-carrier warnings to results_df['warnings'] in place."""
    from ..association.diagnostics import compute_per_gene_warnings

    warnings_by_gene: dict[str, str] = {}
    # Case-carrier warnings are a binary-trait concept.
    for gene in [] if quantitative else results_df["gene"].unique():
        gdata = gene_data_by_gene.get(gene, {})
        case_carriers = gdata.get("proband_carrier_count", 0)
        gene_warnings = compute_per_gene_warnings(gene, case_carriers, assoc_config)
        if gene_warnings:
            warnings_by_gene[gene] = ";".join(gene_warnings)
    # Keep test-reported codes (e.g. linear/logistic burden skips) and append carrier warnings.
    test_codes = results_df.get("warnings", pd.Series("", index=results_df.index))
    results_df["warnings"] = [
        ";".join(c for c in (*(v if isinstance(v, list) else [v] if v else []), w) if c)
        for v, w in zip(
            test_codes, results_df["gene"].map(warnings_by_gene).fillna(""), strict=True
        )
    ]
