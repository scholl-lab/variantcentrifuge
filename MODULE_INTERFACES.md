# VariantCentrifuge Module Interface Specifications

## Overview

This document defines the detailed interfaces and contracts for the proposed modular architecture of VariantCentrifuge. Each module has clearly defined inputs, outputs, and responsibilities.

## Core Interfaces

### 1. Pipeline Orchestrator Interface

```python
# pipeline/orchestrator.py

from typing import Dict, Any, Optional
from dataclasses import dataclass
import argparse
from datetime import datetime

@dataclass
class PipelineResult:
    """Result container for pipeline execution."""
    success: bool
    final_output_path: Optional[str]
    metadata_path: Optional[str]
    report_paths: Dict[str, str]  # {report_type: path}
    errors: List[str]
    warnings: List[str]
    execution_time: float

class IPipelineOrchestrator:
    """Main pipeline orchestration interface."""
    
    def initialize(self, args: argparse.Namespace, config: Dict[str, Any], 
                   start_time: datetime) -> None:
        """Initialize the pipeline with configuration."""
        pass
    
    def validate(self) -> bool:
        """Validate environment and configuration."""
        pass
    
    def execute(self) -> PipelineResult:
        """Execute the complete pipeline."""
        pass
    
    def cleanup(self, keep_intermediates: bool = False) -> None:
        """Clean up temporary files and resources."""
        pass
```

### 2. Configuration Manager Interface

```python
# pipeline/config_manager.py

from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass

@dataclass
class PhenotypeConfig:
    """Phenotype configuration container."""
    phenotypes: Dict[str, str]
    case_terms: List[str]
    control_terms: List[str]

@dataclass
class SampleConfig:
    """Sample configuration container."""
    case_samples: List[str]
    control_samples: List[str]
    all_samples: List[str]
    sample_substring_removal: Optional[str]

@dataclass
class AnnotationConfig:
    """Annotation configuration container."""
    bed_files: List[str]
    gene_lists: List[str]
    json_gene_data: Optional[Dict[str, Any]]
    json_mapping: Optional[Dict[str, Any]]

class IConfigurationManager:
    """Configuration management interface."""
    
    def __init__(self, args: argparse.Namespace, base_config: Dict[str, Any]):
        """Initialize with command line args and base configuration."""
        pass
    
    def load_phenotypes(self) -> PhenotypeConfig:
        """Load and validate phenotype configuration."""
        pass
    
    def load_samples(self, vcf_file: str) -> SampleConfig:
        """Load sample information from VCF and configuration."""
        pass
    
    def load_scoring_config(self) -> Optional[Dict[str, Any]]:
        """Load variant scoring configuration if specified."""
        pass
    
    def load_pedigree_data(self) -> Optional[Dict[str, Any]]:
        """Load pedigree data for inheritance analysis."""
        pass
    
    def get_annotation_config(self) -> AnnotationConfig:
        """Get unified annotation configuration."""
        pass
    
    def get_output_config(self) -> Dict[str, Any]:
        """Get output format and path configuration."""
        pass
```

### 3. Data Processor Interfaces

```python
# processors/base_processor.py

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

class IProcessor(ABC):
    """Base interface for all data processors."""
    
    @abstractmethod
    def process(self, input_data: Any, config: Dict[str, Any]) -> Any:
        """Process input data according to configuration."""
        pass
    
    @abstractmethod
    def validate_input(self, input_data: Any) -> bool:
        """Validate input data before processing."""
        pass
```

```python
# processors/gene_processor.py

from typing import List, Optional
from pathlib import Path

class IGeneProcessor(IProcessor):
    """Gene processing interface."""
    
    def normalize_genes(self, gene_names: Optional[str], 
                       gene_file: Optional[Path]) -> str:
        """Normalize gene names from various input sources."""
        pass
    
    def create_bed_file(self, genes: str, reference: str, 
                       output_dir: Path, interval_expansion: int = 0) -> Path:
        """Create BED file from gene list."""
        pass
    
    def split_bed_file(self, bed_file: Path, chunk_size: int) -> List[Path]:
        """Split BED file into chunks for parallel processing."""
        pass
```

```python
# processors/variant_processor.py

from typing import List, Dict, Optional, Iterator
from pathlib import Path

class IVariantProcessor(IProcessor):
    """Variant extraction and filtering interface."""
    
    def extract_variants_by_region(self, vcf_file: Path, bed_file: Path,
                                  output_file: Path) -> Path:
        """Extract variants within specified regions."""
        pass
    
    def apply_filters(self, vcf_file: Path, filters: str, 
                     output_file: Path) -> Path:
        """Apply SnpSift filters to VCF."""
        pass
    
    def extract_fields(self, vcf_file: Path, fields: List[str],
                      output_file: Path) -> Path:
        """Extract specified fields to TSV format."""
        pass
    
    def apply_bcftools_filter(self, vcf_file: Path, filter_expr: str,
                             output_file: Path) -> Path:
        """Apply bcftools filtering for performance."""
        pass
```

```python
# processors/genotype_processor.py

from typing import Iterator, Tuple, List
from pathlib import Path

class IGenotypeProcessor(IProcessor):
    """Genotype replacement and manipulation interface."""
    
    def replace_genotypes(self, tsv_file: Path, sample_map: Dict[str, str],
                         output_file: Path) -> Path:
        """Replace genotype encodings with sample IDs."""
        pass
    
    def filter_by_genotype(self, tsv_file: Path, 
                          genotype_filters: Dict[str, List[str]],
                          output_file: Path) -> Path:
        """Filter variants by genotype patterns."""
        pass
    
    def extract_sample_genotypes(self, tsv_line: str) -> Dict[str, str]:
        """Extract sample genotypes from a TSV line."""
        pass
```

### 4. Analyzer Interfaces

```python
# analyzers/base_analyzer.py

from abc import ABC, abstractmethod
from typing import Dict, Any, Iterator
import pandas as pd

class IAnalyzer(ABC):
    """Base interface for all analyzers."""
    
    @abstractmethod
    def analyze(self, data: pd.DataFrame, config: Dict[str, Any]) -> pd.DataFrame:
        """Perform analysis on input data."""
        pass
    
    @abstractmethod
    def get_summary(self) -> Dict[str, Any]:
        """Get analysis summary statistics."""
        pass
```

```python
# analyzers/variant_analyzer.py

from typing import Dict, List, Optional
import pandas as pd

class IVariantAnalyzer(IAnalyzer):
    """Variant-level analysis interface."""
    
    def analyze(self, variants_df: pd.DataFrame, 
                config: Dict[str, Any]) -> pd.DataFrame:
        """Perform variant-level analysis."""
        pass
    
    def calculate_allele_frequencies(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Calculate allele frequencies across samples."""
        pass
    
    def annotate_clinical_significance(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Add clinical significance annotations."""
        pass
```

```python
# analyzers/inheritance_analyzer.py

from typing import Dict, List, Optional
import pandas as pd

@dataclass
class InheritanceResult:
    """Container for inheritance analysis results."""
    pattern: str
    confidence: float
    affected_samples: List[str]
    description: str
    supporting_variants: List[str]

class IInheritanceAnalyzer(IAnalyzer):
    """Inheritance pattern analysis interface."""
    
    def analyze_inheritance(self, variants_df: pd.DataFrame,
                           pedigree_data: Dict[str, Any],
                           config: Dict[str, Any]) -> pd.DataFrame:
        """Analyze inheritance patterns in variants."""
        pass
    
    def detect_compound_heterozygous(self, gene_variants: pd.DataFrame,
                                   pedigree_data: Dict[str, Any]) -> List[InheritanceResult]:
        """Detect compound heterozygous variants."""
        pass
    
    def calculate_segregation(self, variant: pd.Series,
                            pedigree_data: Dict[str, Any]) -> float:
        """Calculate segregation score for a variant."""
        pass
```

```python
# analyzers/gene_burden_analyzer.py

from typing import Dict, List, Tuple
import pandas as pd

@dataclass
class GeneBurdenResult:
    """Container for gene burden test results."""
    gene: str
    case_variants: int
    control_variants: int
    odds_ratio: float
    p_value: float
    q_value: float

class IGeneBurdenAnalyzer(IAnalyzer):
    """Gene burden analysis interface."""
    
    def calculate_gene_burden(self, variants_df: pd.DataFrame,
                            case_samples: List[str],
                            control_samples: List[str]) -> List[GeneBurdenResult]:
        """Perform gene burden analysis."""
        pass
    
    def apply_multiple_testing_correction(self, 
                                        p_values: List[float]) -> List[float]:
        """Apply multiple testing correction."""
        pass
```

### 5. Reporter Interfaces

```python
# reporters/base_reporter.py

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional
from pathlib import Path

class IReporter(ABC):
    """Base interface for all output reporters."""
    
    @abstractmethod
    def generate(self, data: Any, output_path: Path, 
                config: Dict[str, Any]) -> Path:
        """Generate report output."""
        pass
    
    @abstractmethod
    def validate_output(self, output_path: Path) -> bool:
        """Validate generated output."""
        pass
```

```python
# reporters/tsv_reporter.py

from typing import List, Iterator
import pandas as pd
from pathlib import Path

class ITSVReporter(IReporter):
    """TSV output generation interface."""
    
    def write_tsv(self, data: pd.DataFrame, output_path: Path,
                 compression: Optional[str] = None) -> Path:
        """Write data to TSV format."""
        pass
    
    def write_streaming(self, data_iterator: Iterator[str], 
                       output_path: Path) -> Path:
        """Write data in streaming fashion for large files."""
        pass
    
    def add_variant_ids(self, data: pd.DataFrame) -> pd.DataFrame:
        """Add unique variant identifiers."""
        pass
```

```python
# reporters/excel_reporter.py

from typing import Dict, List
import pandas as pd
from pathlib import Path

class IExcelReporter(IReporter):
    """Excel output generation interface."""
    
    def create_workbook(self, output_path: Path) -> None:
        """Create new Excel workbook."""
        pass
    
    def add_sheet(self, data: pd.DataFrame, sheet_name: str,
                 formatting: Optional[Dict[str, Any]] = None) -> None:
        """Add sheet to workbook with optional formatting."""
        pass
    
    def add_hyperlinks(self, sheet_name: str, column: str,
                      link_template: str) -> None:
        """Add hyperlinks to specified column."""
        pass
    
    def finalize(self) -> Path:
        """Finalize and save the workbook."""
        pass
```

```python
# reporters/html_reporter.py

from typing import Dict, Any, List
import pandas as pd
from pathlib import Path

class IHTMLReporter(IReporter):
    """HTML report generation interface."""
    
    def generate_report(self, variants_data: pd.DataFrame,
                       summary_data: Dict[str, Any],
                       output_dir: Path,
                       config: Dict[str, Any]) -> Path:
        """Generate interactive HTML report."""
        pass
    
    def create_variant_table(self, data: pd.DataFrame) -> str:
        """Create HTML table for variants."""
        pass
    
    def create_summary_charts(self, summary_data: Dict[str, Any]) -> List[str]:
        """Create summary visualization charts."""
        pass
    
    def integrate_igv_links(self, igv_map: Dict[str, str]) -> None:
        """Integrate IGV report links."""
        pass
```

### 6. Utility Interfaces

```python
# utils/file_utils.py

from typing import Union, Iterator, Optional
from pathlib import Path
import io

class IFileHandler:
    """File handling utilities interface."""
    
    def smart_open(self, file_path: Union[str, Path], mode: str = 'r',
                  encoding: str = 'utf-8') -> io.IOBase:
        """Open file with automatic compression detection."""
        pass
    
    def read_lines(self, file_path: Path) -> Iterator[str]:
        """Read file lines with memory efficiency."""
        pass
    
    def get_file_info(self, file_path: Path) -> Dict[str, Any]:
        """Get file metadata and statistics."""
        pass
```

```python
# utils/validation_utils.py

from typing import List, Dict, Any, Tuple
from pathlib import Path

class IValidator:
    """Input validation utilities interface."""
    
    def validate_vcf(self, vcf_path: Path) -> Tuple[bool, List[str]]:
        """Validate VCF file format and content."""
        pass
    
    def validate_bed(self, bed_path: Path) -> Tuple[bool, List[str]]:
        """Validate BED file format."""
        pass
    
    def validate_tools(self, required_tools: List[str]) -> Dict[str, bool]:
        """Check availability of required external tools."""
        pass
    
    def validate_config(self, config: Dict[str, Any], 
                       schema: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate configuration against schema."""
        pass
```

## Data Flow Contracts

### 1. Pipeline Execution Flow

```python
# Example of data flow through the pipeline

# 1. Configuration Phase
config_manager = ConfigurationManager(args, base_config)
phenotype_config = config_manager.load_phenotypes()
sample_config = config_manager.load_samples(vcf_file)
annotation_config = config_manager.get_annotation_config()

# 2. Processing Phase
gene_processor = GeneProcessor()
normalized_genes = gene_processor.normalize_genes(args.gene_name, args.gene_file)
bed_file = gene_processor.create_bed_file(normalized_genes, reference, output_dir)

variant_processor = VariantProcessor()
variants_vcf = variant_processor.extract_variants_by_region(vcf_file, bed_file, output_path)
filtered_vcf = variant_processor.apply_filters(variants_vcf, filters, output_path)
variants_tsv = variant_processor.extract_fields(filtered_vcf, fields, output_path)

# 3. Analysis Phase
variants_df = pd.read_csv(variants_tsv, sep='\t')

variant_analyzer = VariantAnalyzer()
analyzed_df = variant_analyzer.analyze(variants_df, config)

if pedigree_data:
    inheritance_analyzer = InheritanceAnalyzer()
    analyzed_df = inheritance_analyzer.analyze_inheritance(analyzed_df, pedigree_data, config)

# 4. Reporting Phase
tsv_reporter = TSVReporter()
final_tsv = tsv_reporter.write_tsv(analyzed_df, output_path)

if config.get('xlsx'):
    excel_reporter = ExcelReporter()
    excel_reporter.create_workbook(excel_path)
    excel_reporter.add_sheet(analyzed_df, "Variants")
    excel_path = excel_reporter.finalize()

if config.get('html_report'):
    html_reporter = HTMLReporter()
    html_path = html_reporter.generate_report(analyzed_df, summary_data, report_dir, config)
```

### 2. Error Handling Contract

```python
class PipelineError(Exception):
    """Base exception for pipeline errors."""
    pass

class ConfigurationError(PipelineError):
    """Configuration-related errors."""
    pass

class ProcessingError(PipelineError):
    """Data processing errors."""
    pass

class AnalysisError(PipelineError):
    """Analysis-related errors."""
    pass

class ReportingError(PipelineError):
    """Report generation errors."""
    pass

# All modules should follow this error handling pattern
try:
    result = processor.process(input_data, config)
except ValidationError as e:
    logger.error(f"Validation failed: {e}")
    raise ProcessingError(f"Failed to process: {e}") from e
except Exception as e:
    logger.error(f"Unexpected error: {e}")
    raise ProcessingError(f"Processing failed: {e}") from e
```

### 3. Checkpoint Contract

```python
# pipeline/checkpoint_manager.py

from typing import Dict, Any, Optional
from pathlib import Path
import json

@dataclass
class CheckpointState:
    """State information for pipeline checkpointing."""
    step_name: str
    status: str  # 'pending', 'in_progress', 'completed', 'failed'
    start_time: Optional[datetime]
    end_time: Optional[datetime]
    input_files: List[Path]
    output_files: List[Path]
    metadata: Dict[str, Any]

class ICheckpointManager:
    """Checkpoint and resume functionality interface."""
    
    def __init__(self, output_dir: Path, enable_checksum: bool = False):
        """Initialize checkpoint manager."""
        pass
    
    def initialize(self, config: Dict[str, Any], version: str) -> None:
        """Initialize new checkpoint session."""
        pass
    
    def load(self) -> bool:
        """Load existing checkpoint state."""
        pass
    
    def save_checkpoint(self, state: CheckpointState) -> None:
        """Save checkpoint state."""
        pass
    
    def should_skip_step(self, step_name: str) -> bool:
        """Check if step should be skipped."""
        pass
    
    def mark_step_started(self, step_name: str) -> None:
        """Mark step as started."""
        pass
    
    def mark_step_completed(self, step_name: str, 
                           output_files: List[Path]) -> None:
        """Mark step as completed."""
        pass
    
    def mark_step_failed(self, step_name: str, error: str) -> None:
        """Mark step as failed."""
        pass
    
    def can_resume(self, config: Dict[str, Any], version: str) -> bool:
        """Check if pipeline can be resumed."""
        pass
```

## Module Communication Patterns

### 1. Event-Based Communication

```python
# pipeline/events.py

from enum import Enum
from typing import Any, Callable, Dict

class PipelineEvent(Enum):
    """Pipeline event types."""
    STEP_STARTED = "step_started"
    STEP_COMPLETED = "step_completed"
    STEP_FAILED = "step_failed"
    PROGRESS_UPDATE = "progress_update"
    WARNING = "warning"

class IEventBus:
    """Event communication interface."""
    
    def subscribe(self, event_type: PipelineEvent, 
                 handler: Callable[[Dict[str, Any]], None]) -> None:
        """Subscribe to pipeline events."""
        pass
    
    def publish(self, event_type: PipelineEvent, 
                data: Dict[str, Any]) -> None:
        """Publish pipeline event."""
        pass
```

### 2. Shared Data Store

```python
# pipeline/data_store.py

from typing import Any, Dict, Optional
import threading

class IPipelineDataStore:
    """Thread-safe shared data store for pipeline components."""
    
    def __init__(self):
        self._store: Dict[str, Any] = {}
        self._lock = threading.Lock()
    
    def set(self, key: str, value: Any) -> None:
        """Store a value."""
        pass
    
    def get(self, key: str, default: Any = None) -> Any:
        """Retrieve a value."""
        pass
    
    def exists(self, key: str) -> bool:
        """Check if key exists."""
        pass
    
    def update(self, key: str, updater: Callable[[Any], Any]) -> Any:
        """Atomically update a value."""
        pass
```

## Testing Interfaces

### 1. Mock Factories

```python
# tests/mocks/mock_factories.py

from typing import Dict, Any
from unittest.mock import Mock

class MockFactory:
    """Factory for creating mock objects for testing."""
    
    @staticmethod
    def create_mock_processor() -> IProcessor:
        """Create mock processor for testing."""
        pass
    
    @staticmethod
    def create_mock_analyzer() -> IAnalyzer:
        """Create mock analyzer for testing."""
        pass
    
    @staticmethod
    def create_mock_reporter() -> IReporter:
        """Create mock reporter for testing."""
        pass
    
    @staticmethod
    def create_test_dataframe(num_variants: int = 100) -> pd.DataFrame:
        """Create test variant dataframe."""
        pass
```

### 2. Integration Test Helpers

```python
# tests/integration/helpers.py

from typing import Dict, Any, List
from pathlib import Path

class IntegrationTestHelper:
    """Helper utilities for integration testing."""
    
    @staticmethod
    def setup_test_environment(test_dir: Path) -> Dict[str, Path]:
        """Set up test environment with sample files."""
        pass
    
    @staticmethod
    def create_test_vcf(output_path: Path, num_variants: int = 100) -> Path:
        """Create test VCF file."""
        pass
    
    @staticmethod
    def create_test_config(overrides: Dict[str, Any] = None) -> Dict[str, Any]:
        """Create test configuration."""
        pass
    
    @staticmethod
    def compare_outputs(expected: Path, actual: Path, 
                       ignore_columns: List[str] = None) -> bool:
        """Compare output files for testing."""
        pass
```

## Performance Contracts

### 1. Memory Usage

- Streaming processors should not load entire files into memory
- Use generators for large file processing
- Maximum memory per component: 2GB (configurable)

### 2. Processing Speed

- Variant extraction: > 10,000 variants/second
- TSV parsing: > 50,000 lines/second
- Report generation: < 10 seconds for 100,000 variants

### 3. Parallelization

- All processors must be thread-safe
- Support for multiprocessing with shared memory
- Efficient work distribution for parallel execution

## Conclusion

These interfaces provide clear contracts between components, enabling:

1. **Independent Development**: Teams can work on different modules
2. **Easy Testing**: Clear interfaces enable comprehensive mocking
3. **Flexibility**: Components can be swapped or extended
4. **Maintainability**: Clear boundaries and responsibilities

The modular design with well-defined interfaces will make the codebase more maintainable, testable, and extensible while preserving all existing functionality.