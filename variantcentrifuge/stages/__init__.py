"""
Pipeline stages for VariantCentrifuge.

This package contains all stage implementations organized by category:
- setup_stages: Configuration and data loading stages
- processing_stages: Variant extraction and data transformation stages
- analysis_stages: Analysis and computation stages
- output_stages: Report generation and output stages
"""

from .analysis_stages import (
    ChunkedAnalysisStage,
    CustomAnnotationStage,
    DataFrameLoadingStage,
    GeneBurdenAnalysisStage,
    InheritanceAnalysisStage,
    StatisticsGenerationStage,
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
    ExtraColumnRemovalStage,
    FieldExtractionStage,
    GeneBedCreationStage,
    GenotypeReplacementStage,
    MultiAllelicSplitStage,
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
    "ChunkedAnalysisStage",
    # Setup stages
    "ConfigurationLoadingStage",
    "CustomAnnotationStage",
    # Analysis stages
    "DataFrameLoadingStage",
    "ExcelReportStage",
    "ExtraColumnRemovalStage",
    "FieldExtractionStage",
    "FinalFilteringStage",
    # Processing stages
    "GeneBedCreationStage",
    "GeneBurdenAnalysisStage",
    "GenotypeReplacementStage",
    "HTMLReportStage",
    "IGVReportStage",
    "InheritanceAnalysisStage",
    "MetadataGenerationStage",
    "MultiAllelicSplitStage",
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
    "VariantExtractionStage",
    # Output stages
    "VariantIdentifierStage",
    "VariantScoringStage",
]
