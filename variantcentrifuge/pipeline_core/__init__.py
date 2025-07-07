"""
Pipeline infrastructure for VariantCentrifuge.

This package provides the core abstractions for the modular pipeline architecture:
- PipelineContext: Container for all pipeline state and data
- Stage: Abstract base class for all pipeline components
- Workspace: Centralized file path management
- PipelineRunner: Orchestrates stage execution with parallelization
"""

from .context import PipelineContext
from .runner import PipelineRunner
from .stage import Stage
from .workspace import Workspace

__all__ = [
    "PipelineContext",
    "Stage",
    "Workspace",
    "PipelineRunner",
]

# Version of the pipeline infrastructure
__version__ = "1.0.0"
