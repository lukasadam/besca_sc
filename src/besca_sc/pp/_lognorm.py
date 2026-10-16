"""Log-normalization helper for single-cell data."""

from __future__ import annotations

import scanpy as sc
from anndata import AnnData

def lognorm(adata: AnnData, target_sum: float | None = None) -> AnnData:
    """
    Log-normalize the data in an AnnData object.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix to be log-normalized.
    target_sum : float, optional
        The target sum for normalization. Default is None.

    Returns
    -------
    AnnData
        The log-normalized annotated data matrix.
    """
    # First make copy of adata to avoid modifying the original object
    adata = adata.copy()
    # Check whether adata.raw is not None
    if adata.raw is None:
        adata.raw = adata.copy()
        adata.layers["counts"] = adata.X.copy()
    # Normalize by library size and log-transform the data
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    return adata