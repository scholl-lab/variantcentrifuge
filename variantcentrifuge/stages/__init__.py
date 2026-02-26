"""
Pipeline stages for VariantCentrifuge.

This package contains all stage implementations organized by category:
- setup_stages: Configuration and data loading stages
- processing_stages: Variant extraction and data transformation stages
- analysis_stages: Analysis and computation stages
- output_stages: Report generation and output stages
"""

from .analysis_stages import (
    AssociationAnalysisStage,
    ChunkedAnalysisStage,
    ClinVarPM5Stage,
    CustomAnnotationStage,
    DataFrameLoadingStage,
    GeneBurdenAnalysisStage,
    InheritanceAnalysisStage,
    StatisticsGenerationStage,
    VariantAnalysisStage,
    VariantScoringStage,
)
from .output_stages import (
    ArchiveCreationStage,
    ExcelReportStage,
    FinalFilteringStage,
    HTMLReportStage,
    IGVReportStage,
    MetadataGenerationStage,
    PseudonymizationStage,
    TSVOutputStage,
    VariantIdentifierStage,
)
from .processing_stages import (
    DataSortingStage,
    ExtraColumnRemovalStage,
    FieldExtractionStage,
    GeneBedCreationStage,
    GenotypeReplacementStage,
    MultiAllelicSplitStage,
    ParallelCompleteProcessingStage,
    PCAComputationStage,
    PhenotypeIntegrationStage,
    SnpSiftFilterStage,
    VariantExtractionStage,
)

# Import all stages for convenient access
from .setup_stages import (
    AnnotationConfigLoadingStage,
    ConfigurationLoadingStage,
    PedigreeLoadingStage,
    PhenotypeCaseControlAssignmentStage,
    PhenotypeLoadingStage,
    SampleConfigLoadingStage,
    ScoringConfigLoadingStage,
)

__all__ = [
    "AnnotationConfigLoadingStage",
    "ArchiveCreationStage",
    "AssociationAnalysisStage",
    "ChunkedAnalysisStage",
    "ClinVarPM5Stage",
    "ConfigurationLoadingStage",
    "CustomAnnotationStage",
    "DataFrameLoadingStage",
    "DataSortingStage",
    "ExcelReportStage",
    "ExtraColumnRemovalStage",
    "FieldExtractionStage",
    "FinalFilteringStage",
    "GeneBedCreationStage",
    "GeneBurdenAnalysisStage",
    "GenotypeReplacementStage",
    "HTMLReportStage",
    "IGVReportStage",
    "InheritanceAnalysisStage",
    "MetadataGenerationStage",
    "MultiAllelicSplitStage",
    "PCAComputationStage",
    "ParallelCompleteProcessingStage",
    "PedigreeLoadingStage",
    "PhenotypeCaseControlAssignmentStage",
    "PhenotypeIntegrationStage",
    "PhenotypeLoadingStage",
    "PseudonymizationStage",
    "SampleConfigLoadingStage",
    "ScoringConfigLoadingStage",
    "SnpSiftFilterStage",
    "StatisticsGenerationStage",
    "TSVOutputStage",
    "VariantAnalysisStage",
    "VariantExtractionStage",
    "VariantIdentifierStage",
    "VariantScoringStage",
]
