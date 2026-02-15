"""
Stage Registry - Centralized management of all pipeline stages.

This module provides a central registry for all pipeline stages, their metadata,
dependencies, and relationships. It enables stage discovery, validation, and
dependency analysis for the selective resume system.
"""

import logging
from collections import defaultdict
from dataclasses import dataclass

from ..pipeline_core.stage import Stage

logger = logging.getLogger(__name__)


@dataclass
class StageInfo:
    """Information about a pipeline stage."""

    name: str
    class_ref: type[Stage]
    category: str
    aliases: list[str]
    description: str
    estimated_runtime: float = 0.0

    def __post_init__(self):
        """Validate stage info after initialization."""
        if not self.name:
            raise ValueError("Stage name cannot be empty")
        if not self.class_ref:
            raise ValueError("Stage class reference cannot be None")
        if not issubclass(self.class_ref, Stage):
            raise ValueError(f"Class {self.class_ref} must inherit from Stage")


class StageRegistry:
    """Central registry for all pipeline stages."""

    def __init__(self):
        """Initialize the stage registry."""
        self._stages: dict[str, StageInfo] = {}
        self._categories: dict[str, list[str]] = defaultdict(list)
        self._aliases: dict[str, str] = {}
        self._initialized = False

    def register_stage(
        self,
        stage_class: type[Stage],
        category: str,
        aliases: list[str] | None = None,
        estimated_runtime: float = 0.0,
    ) -> None:
        """Register a stage with the registry.

        Parameters
        ----------
        stage_class : Type[Stage]
            The stage class to register
        category : str
            Category name (e.g., 'setup', 'processing', 'analysis', 'output')
        aliases : Optional[List[str]]
            Alternative names for the stage
        estimated_runtime : float
            Estimated runtime in seconds
        """
        if not issubclass(stage_class, Stage):
            raise ValueError(f"Class {stage_class} must inherit from Stage")

        # Create temporary instance to get name and description
        try:
            temp_instance = stage_class()
            stage_name = temp_instance.name
            description = temp_instance.description
        except Exception as e:
            logger.warning(f"Could not instantiate {stage_class} to get metadata: {e}")
            stage_name = stage_class.__name__.replace("Stage", "").lower()
            description = f"Stage: {stage_name}"

        # Validate stage name is unique
        if stage_name in self._stages:
            raise ValueError(f"Stage '{stage_name}' is already registered")

        # Process aliases
        aliases = aliases or []

        # Check for alias conflicts
        for alias in aliases:
            if alias in self._aliases:
                existing_stage = self._aliases[alias]
                raise ValueError(f"Alias '{alias}' already maps to stage '{existing_stage}'")
            if alias in self._stages:
                raise ValueError(f"Alias '{alias}' conflicts with stage name")

        # Create stage info
        stage_info = StageInfo(
            name=stage_name,
            class_ref=stage_class,
            category=category,
            aliases=aliases,
            description=description,
            estimated_runtime=estimated_runtime,
        )

        # Register the stage
        self._stages[stage_name] = stage_info
        self._categories[category].append(stage_name)

        # Register aliases
        for alias in aliases:
            self._aliases[alias] = stage_name

        logger.debug(
            f"Registered stage '{stage_name}' in category '{category}' with aliases {aliases}"
        )

    def get_stage_info(self, stage_name: str) -> StageInfo | None:
        """Get information about a stage.

        Parameters
        ----------
        stage_name : str
            Name or alias of the stage

        Returns
        -------
        Optional[StageInfo]
            Stage information or None if not found
        """
        # First try direct name lookup
        if stage_name in self._stages:
            return self._stages[stage_name]

        # Then try alias lookup
        if stage_name in self._aliases:
            canonical_name = self._aliases[stage_name]
            return self._stages[canonical_name]

        return None

    def get_stage_class(self, stage_name: str) -> type[Stage] | None:
        """Get the class for a stage.

        Parameters
        ----------
        stage_name : str
            Name or alias of the stage

        Returns
        -------
        Optional[Type[Stage]]
            Stage class or None if not found
        """
        stage_info = self.get_stage_info(stage_name)
        return stage_info.class_ref if stage_info else None

    def resolve_stage_name(self, stage_name: str) -> str | None:
        """Resolve an alias to the canonical stage name.

        Parameters
        ----------
        stage_name : str
            Name or alias of the stage

        Returns
        -------
        Optional[str]
            Canonical stage name or None if not found
        """
        if stage_name in self._stages:
            return stage_name
        if stage_name in self._aliases:
            return self._aliases[stage_name]
        return None

    def get_all_stages(self) -> dict[str, StageInfo]:
        """Get all registered stages.

        Returns
        -------
        Dict[str, StageInfo]
            Dictionary of stage name to stage info
        """
        return self._stages.copy()

    def get_stages_by_category(self, category: str) -> dict[str, StageInfo]:
        """Get all stages in a specific category.

        Parameters
        ----------
        category : str
            Category name

        Returns
        -------
        Dict[str, StageInfo]
            Dictionary of stage name to stage info for the category
        """
        stage_names = self._categories.get(category, [])
        return {name: self._stages[name] for name in stage_names}

    def get_categories(self) -> list[str]:
        """Get all available categories.

        Returns
        -------
        List[str]
            List of category names
        """
        return list(self._categories.keys())

    def stage_exists(self, stage_name: str) -> bool:
        """Check if a stage exists.

        Parameters
        ----------
        stage_name : str
            Name or alias of the stage

        Returns
        -------
        bool
            True if stage exists, False otherwise
        """
        return self.get_stage_info(stage_name) is not None

    def build_dependency_graph(self, active_stages: list[Stage]) -> dict[str, set[str]]:
        """Build dependency graph for active stages.

        Parameters
        ----------
        active_stages : List[Stage]
            List of stage instances

        Returns
        -------
        Dict[str, Set[str]]
            Dictionary mapping stage names to their dependencies
        """
        dependency_graph = {}
        active_stage_names = {stage.name for stage in active_stages}

        for stage in active_stages:
            # Get hard dependencies that exist in active stages
            hard_deps = stage.dependencies & active_stage_names

            # Get soft dependencies that exist in active stages
            soft_deps = stage.soft_dependencies & active_stage_names

            # Combine dependencies
            all_deps = hard_deps | soft_deps
            dependency_graph[stage.name] = all_deps

        return dependency_graph

    def validate_stage_chain(
        self, start_stage: str, active_stages: list[Stage]
    ) -> tuple[bool, list[str]]:
        """Validate that a stage chain is valid for execution.

        Parameters
        ----------
        start_stage : str
            Name of the starting stage
        active_stages : List[Stage]
            List of stage instances

        Returns
        -------
        Tuple[bool, List[str]]
            (is_valid, error_messages)
        """
        errors = []

        # Check if start stage exists
        if not self.stage_exists(start_stage):
            errors.append(f"Stage '{start_stage}' does not exist")
            return False, errors

        # Resolve canonical name
        canonical_start = self.resolve_stage_name(start_stage)

        # Build dependency graph
        dependency_graph = self.build_dependency_graph(active_stages)

        # Check if start stage is in active stages
        active_names = {stage.name for stage in active_stages}
        if canonical_start not in active_names:
            errors.append(f"Stage '{canonical_start}' is not in the active pipeline")
            return False, errors

        # Find all stages that need to execute from start_stage onwards
        stages_to_execute = self._find_stages_from(canonical_start, dependency_graph)

        # Validate dependencies for each stage that will execute
        for stage_name in stages_to_execute:
            stage_deps = dependency_graph.get(stage_name, set())

            # Check if all dependencies will be executed or are before start_stage
            for dep in stage_deps:
                if dep not in stages_to_execute:
                    # This dependency won't be executed, which might be a problem
                    errors.append(
                        f"Stage '{stage_name}' depends on '{dep}' which won't be executed "
                        f"when starting from '{canonical_start}'"
                    )

        return len(errors) == 0, errors

    def _find_stages_from(
        self, start_stage: str, dependency_graph: dict[str, set[str]]
    ) -> set[str]:
        """Find all stages that should execute from a given start stage.

        Parameters
        ----------
        start_stage : str
            Starting stage name
        dependency_graph : Dict[str, Set[str]]
            Stage dependency graph

        Returns
        -------
        Set[str]
            Set of stage names that should execute
        """
        stages_to_execute = {start_stage}

        # Find all stages that depend on stages in our execution set
        changed = True
        while changed:
            changed = False
            for stage_name, deps in dependency_graph.items():
                if stage_name not in stages_to_execute and deps & stages_to_execute:
                    # If any dependency is in our execution set, this stage should execute too
                    stages_to_execute.add(stage_name)
                    changed = True

        return stages_to_execute

    def get_stage_summary(self) -> str:
        """Get a summary of all registered stages.

        Returns
        -------
        str
            Human-readable summary
        """
        lines = ["Stage Registry Summary:"]
        lines.append(f"  Total stages: {len(self._stages)}")
        lines.append(f"  Total aliases: {len(self._aliases)}")
        lines.append(f"  Categories: {len(self._categories)}")

        for category, stage_names in self._categories.items():
            lines.append(f"    {category}: {len(stage_names)} stages")
            for stage_name in sorted(stage_names):
                stage_info = self._stages[stage_name]
                aliases_str = (
                    f" (aliases: {', '.join(stage_info.aliases)})" if stage_info.aliases else ""
                )
                lines.append(f"      - {stage_name}{aliases_str}")

        return "\n".join(lines)


# Global registry instance
_registry = StageRegistry()


def get_registry() -> StageRegistry:
    """Get the global stage registry instance.

    Returns
    -------
    StageRegistry
        The global registry instance
    """
    return _registry


def register_stage(
    stage_class: type[Stage],
    category: str,
    aliases: list[str] | None = None,
    estimated_runtime: float = 0.0,
) -> None:
    """Register a stage with the global registry.

    Parameters
    ----------
    stage_class : Type[Stage]
        The stage class to register
    category : str
        Category name
    aliases : Optional[List[str]]
        Alternative names for the stage
    estimated_runtime : float
        Estimated runtime in seconds
    """
    _registry.register_stage(stage_class, category, aliases, estimated_runtime)


def initialize_registry() -> None:
    """Initialize the registry with all available stages."""
    if _registry._initialized:
        return

    logger.info("Initializing stage registry...")

    # Import and register all stages
    _register_setup_stages()
    _register_processing_stages()
    _register_analysis_stages()
    _register_output_stages()

    _registry._initialized = True
    logger.info(f"Stage registry initialized with {len(_registry._stages)} stages")


def _register_setup_stages():
    """Register all setup stages."""
    from .setup_stages import (
        AnnotationConfigLoadingStage,
        ConfigurationLoadingStage,
        PedigreeLoadingStage,
        PhenotypeCaseControlAssignmentStage,
        PhenotypeLoadingStage,
        SampleConfigLoadingStage,
        ScoringConfigLoadingStage,
    )

    register_stage(ConfigurationLoadingStage, "setup", ["config", "configuration"], 1.0)
    register_stage(PhenotypeLoadingStage, "setup", ["phenotype"], 2.0)
    register_stage(ScoringConfigLoadingStage, "setup", ["scoring_config"], 1.0)
    register_stage(PedigreeLoadingStage, "setup", ["pedigree", "ped"], 1.0)
    register_stage(AnnotationConfigLoadingStage, "setup", ["annotation_config"], 1.0)
    register_stage(SampleConfigLoadingStage, "setup", ["sample_config", "samples"], 2.0)
    register_stage(PhenotypeCaseControlAssignmentStage, "setup", ["case_control_assignment"], 3.0)


def _register_processing_stages():
    """Register all processing stages."""
    from .processing_stages import (
        BCFToolsPrefilterStage,
        ExtraColumnRemovalStage,
        FieldExtractionStage,
        GeneBedCreationStage,
        GenotypeReplacementStage,
        MultiAllelicSplitStage,
        ParallelCompleteProcessingStage,
        ParallelVariantExtractionStage,
        PhenotypeIntegrationStage,
        SnpSiftFilterStage,
        StreamingDataProcessingStage,
        TranscriptFilterStage,
        VariantExtractionStage,
    )

    register_stage(GeneBedCreationStage, "processing", ["gene_bed", "bed_creation"], 10.0)
    register_stage(VariantExtractionStage, "processing", ["extraction", "variant_extraction"], 30.0)
    register_stage(ParallelVariantExtractionStage, "processing", ["parallel_extraction"], 20.0)
    register_stage(BCFToolsPrefilterStage, "processing", ["bcftools_prefilter", "prefilter"], 15.0)
    register_stage(
        MultiAllelicSplitStage, "processing", ["split_multiallelic", "multiallelic"], 20.0
    )
    register_stage(TranscriptFilterStage, "processing", ["transcript_filter"], 10.0)
    register_stage(SnpSiftFilterStage, "processing", ["snpsift_filter", "filter"], 25.0)
    register_stage(FieldExtractionStage, "processing", ["field_extraction", "extract_fields"], 15.0)
    # Deprecated (Phase 11): Stage no-ops immediately, GT reconstruction deferred to output time.
    # Kept for dependency graph compatibility (other stages reference "genotype_replacement").
    register_stage(
        GenotypeReplacementStage, "processing", ["genotype_replacement", "replace_genotypes"], 20.0
    )
    register_stage(PhenotypeIntegrationStage, "processing", ["phenotype_integration"], 5.0)
    register_stage(ExtraColumnRemovalStage, "processing", ["column_cleanup", "cleanup"], 2.0)
    register_stage(StreamingDataProcessingStage, "processing", ["streaming_processing"], 30.0)
    register_stage(ParallelCompleteProcessingStage, "processing", ["parallel_processing"], 60.0)


def _register_analysis_stages():
    """Register all analysis stages."""
    from .analysis_stages import (
        ChunkedAnalysisStage,
        CustomAnnotationStage,
        DataFrameLoadingStage,
        GeneBurdenAnalysisStage,
        GenotypeFilterStage,
        InheritanceAnalysisStage,
        ParallelAnalysisOrchestrator,
        StatisticsGenerationStage,
        VariantAnalysisStage,
        VariantScoringStage,
    )

    register_stage(DataFrameLoadingStage, "analysis", ["dataframe_loading", "load_data"], 10.0)
    register_stage(CustomAnnotationStage, "analysis", ["custom_annotation", "annotation"], 15.0)
    register_stage(
        InheritanceAnalysisStage, "analysis", ["inheritance_analysis", "inheritance"], 20.0
    )
    register_stage(VariantScoringStage, "analysis", ["variant_scoring", "scoring"], 15.0)
    register_stage(GenotypeFilterStage, "analysis", ["genotype_filter"], 5.0)
    register_stage(StatisticsGenerationStage, "analysis", ["statistics_generation", "stats"], 10.0)
    register_stage(VariantAnalysisStage, "analysis", ["variant_analysis", "analysis"], 25.0)
    register_stage(GeneBurdenAnalysisStage, "analysis", ["gene_burden", "burden_analysis"], 30.0)
    register_stage(ChunkedAnalysisStage, "analysis", ["chunked_analysis"], 60.0)
    register_stage(ParallelAnalysisOrchestrator, "analysis", ["analysis_orchestrator"], 45.0)


def _register_output_stages():
    """Register all output stages."""
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

    register_stage(VariantIdentifierStage, "output", ["variant_id"], 2.0)
    register_stage(FinalFilteringStage, "output", ["final_filtering"], 5.0)
    register_stage(PseudonymizationStage, "output", ["pseudonymization"], 3.0)
    register_stage(TSVOutputStage, "output", ["tsv_output", "output"], 5.0)
    register_stage(ExcelReportStage, "output", ["excel_report", "excel"], 10.0)
    register_stage(HTMLReportStage, "output", ["html_report", "html"], 15.0)
    register_stage(IGVReportStage, "output", ["igv_report", "igv"], 20.0)
    register_stage(MetadataGenerationStage, "output", ["metadata_generation", "metadata"], 2.0)
    register_stage(ArchiveCreationStage, "output", ["archive_creation", "archive"], 10.0)
    register_stage(ParallelReportGenerationStage, "output", ["parallel_reports"], 25.0)
