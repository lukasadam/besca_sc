"""Filtering helpers for single-cell data."""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from loguru import logger

from besca_sc.pp._subset import subset_adata


def _validate_thresholds(thresholds: dict[str, int | float | None]) -> None:
    """Validate supplied filtering thresholds."""
    for name, value in thresholds.items():
        if value is not None and value < 0:
            raise ValueError(f"{name} must be non-negative")

    if (
        thresholds["min_genes"] is not None
        and thresholds["max_genes"] is not None
        and thresholds["min_genes"] > thresholds["max_genes"]
    ):
        raise ValueError("min_genes cannot be greater than max_genes")
    if (
        thresholds["min_counts"] is not None
        and thresholds["max_counts"] is not None
        and thresholds["min_counts"] > thresholds["max_counts"]
    ):
        raise ValueError("min_counts cannot be greater than max_counts")


def _metric(adata: AnnData, name: str, context: str):
    """Return a required Scanpy metric, or raise a useful error."""
    if name not in adata.obs:
        raise ValueError(f"{context} requires adata.obs['{name}']")
    return adata.obs[name].to_numpy()


def filter_cells_and_genes(
    adata: AnnData,
    max_genes: int | None = None,
    min_genes: int | None = None,
    max_counts: int | None = None,
    min_counts: int | None = None,
    max_mito: float | None = None,
    max_ribo: float | None = None,
    exclude_doublets: bool = False,
    min_cells: int | None = None,
) -> AnnData:
    """Filter cells and genes using quality-control metrics.

    Parameters
    ----------
    adata
        Annotated data matrix containing the expression data and metrics.
    max_genes
        Maximum number of expressed genes allowed per cell.
    min_genes
        Minimum number of expressed genes required per cell.
    max_counts
        Maximum total count value allowed per cell.
    min_counts
        Minimum total count value required per cell.
    max_mito
        Maximum mitochondrial fraction allowed per cell. Cells equal to or
        above this value are removed using ``adata.obs['percent_mito']``.
    max_ribo
        Maximum ribosomal fraction allowed per cell. Cells equal to or above
        this value are removed using ``adata.obs['percent_ribo']``.
    exclude_doublets
        Remove cells marked ``True`` in the Scanpy-compatible
        ``adata.obs['predicted_doublet']`` column.
    min_cells
            Minimum number of cells in which a gene must be expressed. The metric
            must be stored in ``adata.var['n_cells']``.

    Returns
    -------
    AnnData
        A copied AnnData object containing the filtered cells and genes.

    Raises
    ------
    TypeError
        If ``adata`` is not an AnnData object.
    ValueError
        If a threshold is invalid or a required quality-control metric is
        missing.

    Examples
    --------
    >>> filtered = filter(
    ...     adata, min_genes=500, max_mito=0.1, max_ribo=0.4, min_cells=3
    ... )
    """
    if not isinstance(adata, AnnData):
        raise TypeError("adata must be an AnnData object")

    thresholds = {
        "max_genes": max_genes,
        "min_genes": min_genes,
        "max_counts": max_counts,
        "min_counts": min_counts,
        "min_cells": min_cells,
        "max_mito": max_mito,
        "max_ribo": max_ribo,
    }
    _validate_thresholds(thresholds)

    logger.info("started with {} cells and {} genes", adata.n_obs, adata.n_vars)

    cell_mask = np.ones(adata.n_obs, dtype=bool)
    if max_genes is not None:
        cell_mask &= _metric(adata, "n_genes_by_counts", "max_genes") <= max_genes
    if min_genes is not None:
        cell_mask &= _metric(adata, "n_genes_by_counts", "min_genes") >= min_genes
    if max_counts is not None:
        cell_mask &= _metric(adata, "total_counts", "max_counts") <= max_counts
    if min_counts is not None:
        cell_mask &= _metric(adata, "total_counts", "min_counts") >= min_counts
    if max_mito is not None:
        cell_mask &= _metric(adata, "pct_counts_mt", "max_mito") < max_mito
    if max_ribo is not None:
        cell_mask &= _metric(adata, "pct_counts_ribo", "max_ribo") < max_ribo
    if exclude_doublets:
        doublets = _metric(adata, "predicted_doublet", "exclude_doublets")
        if not np.issubdtype(doublets.dtype, np.bool_):
            raise TypeError("adata.obs['predicted_doublet'] must contain booleans")
        cell_mask &= ~doublets

    if not np.all(cell_mask):
        logger.info("removed {} cells", int((~cell_mask).sum()))
        adata = subset_adata(adata, cell_mask, axis=0)

    if min_cells is not None:
        if "total_counts" not in adata.var:
            raise ValueError("min_cells requires adata.var['total_counts']")
        gene_mask = adata.var["total_counts"] >= min_cells
        logger.info("removed {} genes", int((~gene_mask).sum()))
        adata = subset_adata(adata, gene_mask, axis=1)

    logger.info("finished with {} cells and {} genes", adata.n_obs, adata.n_vars)
    return adata
