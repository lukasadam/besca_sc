"""BESCA: Utilities for single-cell transcriptomics analysis.

Provides modules for preprocessing, quality control, downstream analysis,
visualization, and single-cell dataset management.
"""

from besca_sc.pp._filter import filter as filter_cells  

__all__: list[str] = ["filter_cells"]
