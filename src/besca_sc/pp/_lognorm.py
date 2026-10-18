"""Log-normalization helper for single-cell data."""

from __future__ import annotations

import scanpy as sc
from anndata import AnnData
from loguru import logger
from besca_sc._helper import is_count_matrix


def lognorm(adata: AnnData, target_sum: float | None = 1e4) -> AnnData:
    """
    Log-normalize the data in an AnnData object.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix to be log-normalized.
    target_sum : float, optional
        The target sum for normalization. Default is 1e4.

    Returns
    -------
    AnnData
        The log-normalized annotated data matrix.
    """
    logger.info(
        "started log-normalization for {} cells and {} genes",
        adata.n_obs,
        adata.n_vars,
    )

    if not is_count_matrix(adata.X):
        raise ValueError(
            "adata.X appears to be already normalized or contains non-integer values."
        )

    # First make copy of adata to avoid modifying the original object
    adata = adata.copy()
    # Check whether adata.raw is not None
    if adata.raw is None:
        adata.raw = adata.copy()
    # Store the original counts in a separate layer
    adata.layers["counts"] = adata.X.copy()
    # Normalize by library size and log-transform the data
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    logger.info(
        "finished log-normalization for {} cells and {} genes",
        adata.n_obs,
        adata.n_vars,
    )
    return adata
