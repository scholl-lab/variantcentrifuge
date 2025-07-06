# Pipeline.py Refactoring Plan - Final Architecture

## Executive Summary

This document presents the final refactoring plan for `pipeline.py`, incorporating senior developer feedback to create a unified, simplified architecture using a single Stage abstraction and PipelineContext pattern.

## Core Architecture

### 1. PipelineContext - Single Source of Truth

```python
# pipeline/context.py
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, List, Set
from pathlib import Path
import pandas as pd
from datetime import datetime
import argparse

@dataclass
class PipelineContext:
    """Container for all pipeline state and data - the single source of truth."""
    
    # --- Immutable Configuration ---
    args: argparse.Namespace
    config: Dict[str, Any]
    workspace: 'Workspace'  # Manages all file paths
    start_time: datetime = field(default_factory=datetime.now)
    
    # --- Mutable State & Data Artifacts ---
    # Core data that flows through pipeline
    data: Any = None  # Can be file path or DataFrame
    
    # Stage tracking
    completed_stages: Set[str] = field(default_factory=set)
    stage_results: Dict[str, Any] = field(default_factory=dict)
    
    # Loaded configurations
    vcf_samples: List[str] = field(default_factory=list)
    pedigree_data: Optional[Dict] = None
    phenotype_data: Optional[Dict] = None
    scoring_config: Optional[Dict] = None
    annotation_configs: Dict[str, Any] = field(default_factory=dict)
    
    # Processing artifacts
    gene_bed_file: Optional[Path] = None
    extracted_vcf: Optional[Path] = None
    filtered_vcf: Optional[Path] = None
    extracted_tsv: Optional[Path] = None
    
    # Analysis results
    current_dataframe: Optional[pd.DataFrame] = None
    statistics: Dict[str, Any] = field(default_factory=dict)
    gene_burden_results: Optional[pd.DataFrame] = None
    
    # Output paths
    final_output_path: Optional[Path] = None
    report_paths: Dict[str, Path] = field(default_factory=dict)
    
    # Checkpoint state (if enabled)
    checkpoint_state: Optional[Any] = None
    
    def mark_complete(self, stage_name: str, result: Any = None) -> None:
        """Mark a stage as complete with optional result storage."""
        self.completed_stages.add(stage_name)
        if result is not None:
            self.stage_results[stage_name] = result
    
    def is_complete(self, stage_name: str) -> bool:
        """Check if a stage has been completed."""
        return stage_name in self.completed_stages
```

### 2. Unified Stage Abstraction

```python
# pipeline/stage.py
from abc import ABC, abstractmethod
from typing import Set, List, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class Stage(ABC):
    """Abstract base class for all pipeline stages - unified abstraction."""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for the stage."""
        pass
    
    @property
    def dependencies(self) -> Set[str]:
        """Stage names that must complete before this stage."""
        return set()
    
    @property
    def description(self) -> str:
        """Human-readable description for logging."""
        return f"Stage: {self.name}"
    
    def __call__(self, context: PipelineContext) -> PipelineContext:
        """Execute the stage with pre/post processing."""
        # Validate dependencies
        for dep in self.dependencies:
            if not context.is_complete(dep):
                raise RuntimeError(
                    f"Stage '{self.name}' requires '{dep}' to be completed first"
                )
        
        # Log execution
        logger.info(f"Executing {self.description}")
        
        # Process
        try:
            updated_context = self._process(context)
            context.mark_complete(self.name)
            return updated_context
        except Exception as e:
            logger.error(f"Stage '{self.name}' failed: {e}")
            raise
    
    @abstractmethod
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Core processing logic - must be implemented by subclasses."""
        pass
    
    def get_input_files(self, context: PipelineContext) -> List[Path]:
        """Return input files for checkpoint tracking."""
        return []
    
    def get_output_files(self, context: PipelineContext) -> List[Path]:
        """Return output files for checkpoint tracking."""
        return []
```

### 3. Workspace Management

```python
# pipeline/workspace.py
from pathlib import Path
from datetime import datetime
import tempfile

class Workspace:
    """Manages all file paths for a pipeline run."""
    
    def __init__(self, output_dir: Path, base_name: str):
        self.output_dir = Path(output_dir)
        self.base_name = base_name
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.intermediate_dir = self.output_dir / "intermediate"
        self.intermediate_dir.mkdir(exist_ok=True)
        self.temp_dir = Path(tempfile.mkdtemp(prefix="vc_", dir=self.intermediate_dir))
    
    def get_output_path(self, suffix: str, extension: str = ".tsv") -> Path:
        """Generate output file path."""
        return self.output_dir / f"{self.base_name}{suffix}{extension}"
    
    def get_intermediate_path(self, name: str) -> Path:
        """Generate intermediate file path."""
        return self.intermediate_dir / name
    
    def get_temp_path(self, name: str) -> Path:
        """Generate temporary file path."""
        return self.temp_dir / name
    
    def cleanup(self):
        """Clean up temporary files."""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
```

### 4. Pipeline Runner

```python
# pipeline/runner.py
from typing import List, Dict, Set
import logging
from .stage import Stage
from .context import PipelineContext

logger = logging.getLogger(__name__)

class PipelineRunner:
    """Executes stages in dependency order."""
    
    def __init__(self, enable_checkpoints: bool = False):
        self.enable_checkpoints = enable_checkpoints
    
    def run(self, stages: List[Stage], context: PipelineContext) -> PipelineContext:
        """Execute all stages in dependency order."""
        # Build dependency graph
        stage_map = {stage.name: stage for stage in stages}
        
        # Topological sort
        sorted_stages = self._topological_sort(stages)
        
        # Execute in order
        for stage in sorted_stages:
            if context.checkpoint_state and context.is_complete(stage.name):
                logger.info(f"Skipping completed stage: {stage.name}")
                continue
                
            context = stage(context)
            
            # Save checkpoint if enabled
            if self.enable_checkpoints and context.checkpoint_state:
                context.checkpoint_state.save()
        
        return context
    
    def _topological_sort(self, stages: List[Stage]) -> List[Stage]:
        """Sort stages by dependencies."""
        # Build adjacency list
        graph = {stage.name: list(stage.dependencies) for stage in stages}
        stage_map = {stage.name: stage for stage in stages}
        
        # Kahn's algorithm
        in_degree = {name: 0 for name in graph}
        for deps in graph.values():
            for dep in deps:
                if dep in in_degree:
                    in_degree[dep] += 1
        
        queue = [name for name, degree in in_degree.items() if degree == 0]
        sorted_names = []
        
        while queue:
            current = queue.pop(0)
            sorted_names.append(current)
            
            for name, deps in graph.items():
                if current in deps:
                    in_degree[name] -= 1
                    if in_degree[name] == 0:
                        queue.append(name)
        
        if len(sorted_names) != len(stages):
            raise ValueError("Circular dependency detected in stages")
        
        return [stage_map[name] for name in sorted_names]
```

## Stage Implementations

### Example: Configuration Loading Stage

```python
# stages/setup_stages.py
from pathlib import Path
from ..pipeline.stage import Stage
from ..pipeline.context import PipelineContext

class ConfigurationLoadingStage(Stage):
    """Load and merge configuration from file and CLI arguments."""
    
    @property
    def name(self) -> str:
        return "configuration_loading"
    
    @property
    def description(self) -> str:
        return "Load configuration and merge with CLI arguments"
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # Configuration loading logic moved from pipeline.py
        config = load_config_file(context.args.config)
        
        # Merge CLI args (CLI takes precedence)
        for key, value in vars(context.args).items():
            if value is not None:
                config[key] = value
        
        # Update context
        context.config = config
        return context
```

### Example: Parallel Variant Extraction Stage

```python
# stages/extraction_stages.py
from concurrent.futures import ProcessPoolExecutor
from ..pipeline.stage import Stage

class ParallelVariantExtractionStage(Stage):
    """Extract variants in parallel by genomic region."""
    
    @property
    def name(self) -> str:
        return "variant_extraction"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"gene_bed_creation"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # All parallel complexity hidden within this stage
        bed_regions = self._split_bed_file(context.gene_bed_file)
        
        with ProcessPoolExecutor(max_workers=context.config.get('threads', 1)) as executor:
            # Submit extraction jobs
            futures = []
            for region_bed in bed_regions:
                future = executor.submit(
                    self._extract_region,
                    context.args.vcf_file,
                    region_bed,
                    context.workspace.get_temp_path(f"region_{len(futures)}.vcf.gz")
                )
                futures.append(future)
            
            # Collect results
            extracted_files = [f.result() for f in futures]
        
        # Merge results
        merged_vcf = self._merge_vcfs(
            extracted_files,
            context.workspace.get_intermediate_path("variants.vcf.gz")
        )
        
        context.extracted_vcf = merged_vcf
        context.data = merged_vcf  # Update primary data artifact
        return context
```

## Simplified Main Pipeline

```python
# pipeline.py (new, lean orchestrator)
import logging
from pathlib import Path
from .pipeline import PipelineContext, PipelineRunner, Workspace
from .stages import (
    # Setup stages
    ConfigurationLoadingStage,
    PhenotypeLoadingStage,
    ScoringConfigLoadingStage,
    
    # Processing stages  
    GeneBedCreationStage,
    ParallelVariantExtractionStage,
    FieldExtractionStage,
    GenotypeReplacementStage,
    
    # Analysis stages
    DataFrameLoadingStage,
    CustomAnnotationStage,
    InheritanceAnalysisStage,
    VariantScoringStage,
    
    # Output stages
    TSVOutputStage,
    ExcelReportStage,
    HTMLReportStage
)

logger = logging.getLogger(__name__)

def run_pipeline(args):
    """Main pipeline orchestrator - clean and simple."""
    # Initialize context
    workspace = Workspace(args.output_dir, args.base_name)
    context = PipelineContext(args=args, config={}, workspace=workspace)
    
    # Assemble stages based on configuration
    stages = [
        # Always run
        ConfigurationLoadingStage(),
        
        # Conditional stages
        PhenotypeLoadingStage() if args.phenotype_file else None,
        ScoringConfigLoadingStage() if args.scoring_config else None,
        
        # Core processing
        GeneBedCreationStage(),
        ParallelVariantExtractionStage(),
        FieldExtractionStage(),
        GenotypeReplacementStage(),
        
        # Analysis
        DataFrameLoadingStage(),
        CustomAnnotationStage(),
        InheritanceAnalysisStage() if args.ped_file else None,
        VariantScoringStage() if args.scoring_config else None,
        
        # Output
        TSVOutputStage(),
        ExcelReportStage() if args.xlsx else None,
        HTMLReportStage() if args.html_report else None,
    ]
    
    # Filter out None stages
    stages = [s for s in stages if s is not None]
    
    # Run pipeline
    runner = PipelineRunner(enable_checkpoints=args.enable_checkpoints)
    try:
        final_context = runner.run(stages, context)
        logger.info("Pipeline completed successfully")
        return final_context
    finally:
        workspace.cleanup()
```

## Benefits of the Refined Architecture

1. **Unified Abstraction**: Single `Stage` pattern for all components
2. **Explicit Data Flow**: `PipelineContext` passed through each stage
3. **Hidden Complexity**: Parallelization and chunking encapsulated within stages
4. **Simple Testing**: Each stage easily testable in isolation
5. **Clean Orchestrator**: Main pipeline file reduced to ~100 lines

## Implementation Strategy

### Phase 1: Core Infrastructure (Days 1-3)
- Implement PipelineContext, Stage, Workspace, and PipelineRunner
- Create comprehensive unit tests
- Validate with simple test stages

### Phase 2: Stage Migration (Weeks 2-3)
- Migrate functionality from pipeline.py to stages
- Implement stages in dependency order
- Test each stage thoroughly

### Phase 3: Integration (Week 4)
- Connect all stages
- Run full regression tests
- Optimize performance

### Phase 4: Polish (Week 5)
- Complete documentation
- Final testing
- Prepare for merge

## Success Criteria

✅ All regression tests pass (byte-for-byte output)
✅ Pipeline.py reduced from 2,831 to ~200 lines
✅ Each stage < 200 lines
✅ Test coverage > 85%
✅ Performance equal or better

This refined architecture incorporates senior developer feedback to create a simpler, more maintainable solution while preserving all the benefits of the original plan.