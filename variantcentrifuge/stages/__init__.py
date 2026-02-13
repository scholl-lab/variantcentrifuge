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
    ParallelAnalysisOrchestrator,
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
    ParallelReportGenerationStage,
    PseudonymizationStage,
    TSVOutputStage,
    VariantIdentifierStage,
)
from .processing_stages import (
    BCFToolsPrefilterStage,
    ExtraColumnRemovalStage,
    FieldExtractionStage,
    GeneBedCreationStage,
    GenotypeReplacementStage,
    MultiAllelicSplitStage,
    ParallelVariantExtractionStage,
    PhenotypeIntegrationStage,
    SnpSiftFilterStage,
    StreamingDataProcessingStage,
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
    "BCFToolsPrefilterStage",
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
    "ParallelAnalysisOrchestrator",
    "ParallelReportGenerationStage",
    "ParallelVariantExtractionStage",
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
