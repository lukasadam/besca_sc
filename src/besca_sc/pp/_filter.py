"""Filtering helpers for single-cell data."""

from __future__ import annotations

from collections.abc import Iterable

from anndata import AnnData
from loguru import logger


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


def _metric(adata: AnnData, names: Iterable[str], label: str):
    """Return the first available metric, or raise a useful error."""
    for name in names:
        if name in adata.obs:
            return adata.obs[name]
    available = ", ".join(names)
    raise ValueError(f"{label} requires one of these adata.obs columns: {available}")


def _filter_cells(
    adata: AnnData, metric, condition, description: str, threshold
) -> AnnData:
    """Apply a cell mask and log the number of removed cells."""
    before = adata.n_obs
    filtered = adata[condition(metric), :].copy()
    logger.info(
        "removed {} cells {} {}", before - filtered.n_obs, description, threshold
    )
    return filtered


def filter(
    adata: AnnData,
    max_genes: int | None = None,
    min_genes: int | None = None,
    max_counts: int | None = None,
    min_counts: int | None = None,
    min_cells: int | None = None,
    max_mito: float | None = None,
    annotation_type: str | None = None,
    species: str = "human",
) -> AnnData:
    """Filter cells and genes using quality-control metrics.

    Filters are applied in this order: maximum genes, minimum genes, maximum
    counts, minimum counts, minimum cells per gene, and maximum mitochondrial
    percentage. The input object is not modified.

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
    min_cells
        Minimum number of cells in which a gene must be expressed. The metric
        must be stored in ``adata.var['n_cells']``.
    max_mito
        Maximum mitochondrial fraction allowed per cell. Cells equal to or
        above this value are removed using ``adata.obs['percent_mito']``.
    annotation_type
        Retained for API compatibility and not used by this function.
    species
        Retained for API compatibility and not used by this function.

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
    >>> filtered = filter(adata, min_genes=500, max_mito=0.1, min_cells=3)
    """
    del annotation_type, species
    if not isinstance(adata, AnnData):
        raise TypeError("adata must be an AnnData object")

    thresholds = {
        "max_genes": max_genes,
        "min_genes": min_genes,
        "max_counts": max_counts,
        "min_counts": min_counts,
        "min_cells": min_cells,
        "max_mito": max_mito,
    }
    _validate_thresholds(thresholds)

    logger.info("started with {} cells and {} genes", adata.n_obs, adata.n_vars)

    if max_genes is not None:
        adata = _filter_cells(
            adata,
            _metric(adata, ("n_genes", "n_genes_by_counts"), "max_genes"),
            lambda values: values <= max_genes,
            "with more than",
            max_genes,
        )
    if min_genes is not None:
        adata = _filter_cells(
            adata,
            _metric(adata, ("n_genes", "n_genes_by_counts"), "min_genes"),
            lambda values: values >= min_genes,
            "with fewer than",
            min_genes,
        )
    if max_counts is not None:
        adata = _filter_cells(
            adata,
            _metric(adata, ("n_counts", "total_counts"), "max_counts"),
            lambda values: values <= max_counts,
            "with more than",
            max_counts,
        )
    if min_counts is not None:
        adata = _filter_cells(
            adata,
            _metric(adata, ("n_counts", "total_counts"), "min_counts"),
            lambda values: values >= min_counts,
            "with fewer than",
            min_counts,
        )
    if min_cells is not None:
        if "n_cells" not in adata.var:
            raise ValueError("min_cells requires adata.var['n_cells']")
        before = adata.n_vars
        adata = adata[:, adata.var["n_cells"] >= min_cells].copy()
        logger.info(
            "removed {} genes expressed in fewer than {} cells",
            before - adata.n_vars,
            min_cells,
        )
    if max_mito is not None:
        adata = _filter_cells(
            adata,
            _metric(adata, ("percent_mito",), "max_mito"),
            lambda values: values < max_mito,
            "with mitochondrial percentage at least",
            max_mito,
        )

    logger.info("finished with {} cells and {} genes", adata.n_obs, adata.n_vars)
    return adata
