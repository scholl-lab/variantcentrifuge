# File: variantcentrifuge/association/__init__.py
# Location: variantcentrifuge/variantcentrifuge/association/__init__.py
"""
variantcentrifuge.association — Modular rare variant association framework.

v0.15.0 introduces a plugin-style association testing system where each
statistical test implements the AssociationTest ABC. Phase 18 delivers the
core abstractions and FisherExactTest (bit-identical reimplementation of
gene_burden.py's Fisher test). Subsequent phases add burden regression
(Phase 19), SKAT (Phases 20-21), and ACAT-O (Phase 22).

Public API
----------
AssociationTest   : Abstract base class for all association tests
TestResult        : Dataclass holding per-gene test results
AssociationConfig : Configuration dataclass with defaults mirroring gene_burden.py
AssociationEngine : Orchestrator: dispatches tests, applies correction, returns DataFrame
apply_correction  : Standalone FDR/Bonferroni correction function
"""

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.correction import apply_correction
from variantcentrifuge.association.engine import AssociationEngine

__all__ = [
    "AssociationConfig",
    "AssociationEngine",
    "AssociationTest",
    "TestResult",
    "apply_correction",
]
