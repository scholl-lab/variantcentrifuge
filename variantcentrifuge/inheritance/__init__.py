"""
Inheritance pattern analysis module for VariantCentrifuge.

This module provides functionality to deduce Mendelian inheritance patterns
from variant data and pedigree information.
"""

from .analyzer import analyze_inheritance

__all__ = ["analyze_inheritance"]
