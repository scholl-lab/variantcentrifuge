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
    PhenotypeLoadingStage,
    SampleConfigLoadingStage,
    ScoringConfigLoadingStage,
)

__all__ = [
    # Setup stages
    "ConfigurationLoadingStage",
    "PhenotypeLoadingStage",
    "ScoringConfigLoadingStage",
    "PedigreeLoadingStage",
    "AnnotationConfigLoadingStage",
    "SampleConfigLoadingStage",
    # Processing stages
    "GeneBedCreationStage",
    "VariantExtractionStage",
    "ParallelVariantExtractionStage",
    "BCFToolsPrefilterStage",
    "MultiAllelicSplitStage",
    "SnpSiftFilterStage",
    "FieldExtractionStage",
    "GenotypeReplacementStage",
    "PhenotypeIntegrationStage",
    "ExtraColumnRemovalStage",
    # Analysis stages
    "DataFrameLoadingStage",
    "CustomAnnotationStage",
    "InheritanceAnalysisStage",
    "VariantScoringStage",
    "StatisticsGenerationStage",
    "GeneBurdenAnalysisStage",
    "ChunkedAnalysisStage",
    "ParallelAnalysisOrchestrator",
    # Output stages
    "VariantIdentifierStage",
    "FinalFilteringStage",
    "PseudonymizationStage",
    "TSVOutputStage",
    "ExcelReportStage",
    "HTMLReportStage",
    "IGVReportStage",
    "MetadataGenerationStage",
    "ArchiveCreationStage",
    "ParallelReportGenerationStage",
]
