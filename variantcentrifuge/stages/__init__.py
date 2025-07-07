"""
Pipeline stages for VariantCentrifuge.

This package contains all stage implementations organized by category:
- setup_stages: Configuration and data loading stages
- processing_stages: Variant extraction and data transformation stages
- analysis_stages: Analysis and computation stages
- output_stages: Report generation and output stages
"""

# Import all stages for convenient access
from .setup_stages import (
    ConfigurationLoadingStage,
    PhenotypeLoadingStage,
    ScoringConfigLoadingStage,
    PedigreeLoadingStage,
    AnnotationConfigLoadingStage,
    SampleConfigLoadingStage,
)

from .processing_stages import (
    GeneBedCreationStage,
    VariantExtractionStage,
    ParallelVariantExtractionStage,
    BCFToolsPrefilterStage,
    MultiAllelicSplitStage,
    SnpSiftFilterStage,
    FieldExtractionStage,
    GenotypeReplacementStage,
    PhenotypeIntegrationStage,
    ExtraColumnRemovalStage,
    StreamingDataProcessingStage,
)

from .analysis_stages import (
    DataFrameLoadingStage,
    CustomAnnotationStage,
    InheritanceAnalysisStage,
    VariantScoringStage,
    StatisticsGenerationStage,
    GeneBurdenAnalysisStage,
    ChunkedAnalysisStage,
    ParallelAnalysisOrchestrator,
)

from .output_stages import (
    VariantIdentifierStage,
    FinalFilteringStage,
    PseudonymizationStage,
    TSVOutputStage,
    ExcelReportStage,
    HTMLReportStage,
    IGVReportStage,
    MetadataGenerationStage,
    ArchiveCreationStage,
    ParallelReportGenerationStage,
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
