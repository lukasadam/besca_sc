"""BESCA: Utilities for single-cell transcriptomics analysis.

Provides modules for preprocessing, quality control, downstream analysis,
visualization, and single-cell dataset management.
"""

from besca_sc.pp._filter import filter_cells_and_genes as filter_cells_and_genes
from besca_sc.pp._subset import subset_adata
from besca_sc.pp._qc import qc
from besca_sc.pl._thresh import plot_thresholds
from besca_sc.pp._outlier import calculate_thresholds, get_threshold
from besca_sc.pp._lognorm import lognorm

__all__: list[str] = [
    "filter_cells_and_genes",
    "subset_adata",
    "qc",
    "plot_thresholds",
    "get_threshold",
    "calculate_thresholds",
    "lognorm",
]
