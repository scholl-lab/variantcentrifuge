"""Analysis stages facade module for backwards compatibility."""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

from ..analyze_variants import analyze_variants, analyze_variants_df
from ..annotator import annotate_dataframe_with_features, load_custom_features
from ..gene_burden import perform_gene_burden_analysis
from ..gt_reconstruction import _find_per_sample_gt_columns, reconstruct_gt_column
from ..inheritance import analyze_inheritance
from ..scoring import apply_scoring
from ..stats_engine import StatsEngine
from .analysis_stage_helpers import (
    VALID_ASSOCIATION_KEYS,
    _build_assoc_config_from_context,
    _count_significant_association_results,
    _resolve_association_workers,
    _validate_association_config_dict,
)
from .association_analysis_stages import AssociationAnalysisStage
from .chunk_helpers import (
    apply_chunk_annotations,
    apply_chunk_inheritance_analysis,
    calculate_inheritance_safe_chunk_size,
    ensure_sorted_by_gene,
    read_tsv_in_gene_chunks,
)
from .chunked_analysis_stages import ChunkedAnalysisStage
from .data_loading_stages import (
    CustomAnnotationStage,
    DataFrameLoadingStage,
)
from .gene_burden_stages import GeneBurdenAnalysisStage
from .genotype_matrix_builder import (
    _STANDARD_COLUMN_ALIASES,
    _create_standard_column_aliases,
    _find_cadd_column,
    _find_effect_column,
    _find_revel_column,
    _GenotypeMatrixBuilder,
)
from .inheritance_and_scoring_stages import (
    InheritanceAnalysisStage,
    VariantScoringStage,
)
from .sample_column_helpers import (
    cleanup_sample_columns,
    create_sample_columns_from_gt,
    create_sample_columns_from_gt_intelligent,
    create_sample_columns_from_gt_vectorized,
)
from .statistics_and_clinvar_stages import (
    ClinVarPM5Stage,
    StatisticsGenerationStage,
)
from .variant_analysis_stages import VariantAnalysisStage

__all__ = [
    "VALID_ASSOCIATION_KEYS",
    "_STANDARD_COLUMN_ALIASES",
    "AssociationAnalysisStage",
    "ChunkedAnalysisStage",
    "ClinVarPM5Stage",
    "CustomAnnotationStage",
    "DataFrameLoadingStage",
    "GeneBurdenAnalysisStage",
    "InheritanceAnalysisStage",
    "StatisticsGenerationStage",
    "StatsEngine",
    "VariantAnalysisStage",
    "VariantScoringStage",
    "_GenotypeMatrixBuilder",
    "_build_assoc_config_from_context",
    "_count_significant_association_results",
    "_create_standard_column_aliases",
    "_find_cadd_column",
    "_find_effect_column",
    "_find_per_sample_gt_columns",
    "_find_revel_column",
    "_resolve_association_workers",
    "_validate_association_config_dict",
    "analyze_inheritance",
    "analyze_variants",
    "analyze_variants_df",
    "annotate_dataframe_with_features",
    "apply_chunk_annotations",
    "apply_chunk_inheritance_analysis",
    "apply_scoring",
    "calculate_inheritance_safe_chunk_size",
    "cleanup_sample_columns",
    "create_sample_columns_from_gt",
    "create_sample_columns_from_gt_intelligent",
    "create_sample_columns_from_gt_vectorized",
    "ensure_sorted_by_gene",
    "load_custom_features",
    "perform_gene_burden_analysis",
    "read_tsv_in_gene_chunks",
    "reconstruct_gt_column",
]
