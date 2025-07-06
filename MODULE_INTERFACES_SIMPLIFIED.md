# Module Interfaces - Simplified Architecture

## Overview

This document defines the simplified module interfaces for the refactored VariantCentrifuge pipeline, incorporating senior developer feedback to use a unified Stage abstraction and PipelineContext pattern.

## Core Components

### 1. PipelineContext (pipeline/context.py)

The single source of truth for all pipeline state and data.

```python
@dataclass
class PipelineContext:
    """Container for all pipeline state and data."""
    # Immutable configuration
    args: argparse.Namespace
    config: Dict[str, Any]
    workspace: Workspace
    
    # Primary data artifact (flows through pipeline)
    data: Any = None
    
    # Stage tracking
    completed_stages: Set[str] = field(default_factory=set)
    stage_results: Dict[str, Any] = field(default_factory=dict)
    
    # All other state attributes...
    
    def mark_complete(self, stage_name: str, result: Any = None) -> None
    def is_complete(self, stage_name: str) -> bool
```

**Key Design Decisions:**
- Single object that flows through all stages
- Explicit data flow (no global state)
- Simple testing - just create context and pass to stage

### 2. Stage Abstract Base Class (pipeline/stage.py)

Unified abstraction for all pipeline components.

```python
class Stage(ABC):
    """Abstract base class for all pipeline stages."""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for the stage."""
        pass
    
    @property
    def dependencies(self) -> Set[str]:
        """Stage names that must complete before this stage."""
        return set()
    
    def __call__(self, context: PipelineContext) -> PipelineContext:
        """Execute the stage."""
        pass
    
    @abstractmethod
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Core processing logic."""
        pass
```

**Key Design Decisions:**
- Single pattern for all components (processors, analyzers, reporters)
- Dependencies declared as simple set of stage names
- Complexity hidden within each stage

### 3. Workspace (pipeline/workspace.py)

Centralized file path management.

```python
class Workspace:
    """Manages all file paths for a pipeline run."""
    
    def __init__(self, output_dir: Path, base_name: str)
    def get_output_path(self, suffix: str, extension: str = ".tsv") -> Path
    def get_intermediate_path(self, name: str) -> Path
    def get_temp_path(self, name: str) -> Path
    def cleanup(self)
```

### 4. PipelineRunner (pipeline/runner.py)

Simple execution engine with dependency resolution.

```python
class PipelineRunner:
    """Executes stages in dependency order."""
    
    def run(self, stages: List[Stage], context: PipelineContext) -> PipelineContext
```

## Stage Categories

### Setup Stages (stages/setup_stages.py)

```python
class ConfigurationLoadingStage(Stage):
    """Load and merge configuration."""
    name = "configuration_loading"
    dependencies = set()  # No dependencies

class PhenotypeLoadingStage(Stage):
    """Load phenotype data."""
    name = "phenotype_loading"
    dependencies = set()  # Can run in parallel with config

class ScoringConfigLoadingStage(Stage):
    """Load scoring configuration."""
    name = "scoring_config_loading"
    dependencies = set()
```

### Processing Stages (stages/processing_stages.py)

```python
class GeneBedCreationStage(Stage):
    """Convert genes to BED intervals."""
    name = "gene_bed_creation"
    dependencies = {"configuration_loading"}

class ParallelVariantExtractionStage(Stage):
    """Extract variants in parallel."""
    name = "variant_extraction"
    dependencies = {"gene_bed_creation"}
    
    # All parallel complexity hidden inside _process()

class FieldExtractionStage(Stage):
    """Extract fields from VCF to TSV."""
    name = "field_extraction"
    dependencies = {"variant_extraction", "snpsift_filtering"}
```

### Analysis Stages (stages/analysis_stages.py)

```python
class DataFrameLoadingStage(Stage):
    """Load TSV into DataFrame or setup chunked processing."""
    name = "dataframe_loading"
    dependencies = {"phenotype_integration"}

class InheritanceAnalysisStage(Stage):
    """Calculate inheritance patterns."""
    name = "inheritance_analysis"
    dependencies = {"dataframe_loading", "custom_annotation"}

class ChunkedAnalysisStage(Stage):
    """Process large files in chunks."""
    name = "chunked_analysis"
    dependencies = {"dataframe_loading"}
    
    # Encapsulates all chunking complexity
```

### Output Stages (stages/output_stages.py)

```python
class TSVOutputStage(Stage):
    """Write final TSV output."""
    name = "tsv_output"
    dependencies = {"final_filtering", "pseudonymization"}

class ParallelReportGenerationStage(Stage):
    """Generate all reports in parallel."""
    name = "parallel_reports"
    dependencies = {"tsv_output"}
    
    # Internally manages parallel report generation
```

## Usage Example

```python
# Simple, clean pipeline definition
def run_pipeline(args):
    # Setup
    workspace = Workspace(args.output_dir, args.base_name)
    context = PipelineContext(args=args, config={}, workspace=workspace)
    
    # Define stages
    stages = [
        ConfigurationLoadingStage(),
        GeneBedCreationStage(),
        ParallelVariantExtractionStage(),
        FieldExtractionStage(),
        InheritanceAnalysisStage(),
        TSVOutputStage(),
    ]
    
    # Run
    runner = PipelineRunner()
    final_context = runner.run(stages, context)
    
    return final_context
```

## Testing Pattern

```python
def test_inheritance_analysis_stage():
    # Setup
    context = PipelineContext(
        args=Mock(),
        config={'inheritance_mode': 'simple'},
        workspace=Mock()
    )
    context.current_dataframe = test_dataframe
    context.mark_complete('dataframe_loading')
    context.mark_complete('custom_annotation')
    
    # Execute
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    # Assert
    assert 'Inheritance_Pattern' in result.current_dataframe.columns
```

## Benefits of Simplified Architecture

1. **Single Pattern**: One Stage class for everything
2. **Explicit Flow**: Context passed through each stage
3. **Easy Testing**: Create context, run stage, check result
4. **Hidden Complexity**: Parallelism encapsulated in stages
5. **Clean Main**: Pipeline.py becomes ~100 lines

This simplified architecture makes the codebase more maintainable and testable while preserving all functionality.