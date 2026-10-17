"""Preprocessing pipeline for single-cell data."""

from __future__ import annotations

from collections.abc import Sequence

import anndata as ad
import scanpy as sc
from loguru import logger

from besca_sc.pp._lognorm import lognorm as _lognorm


def preprocess(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    batch_key: str | None = None,
    covariates: Sequence[str] | None = None,
    layer: str | None = None,
    *,
    lognorm: bool = True,
    hvg: bool = True,
    regression: bool | None = None,
) -> ad.AnnData:
    """Run selected preprocessing steps on an AnnData object.

    Steps run in order: log-normalization,
    highly variable gene identification, and regression.
    Regression is enabled automatically when ``covariates``
    is provided, preserving the original behavior.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix to preprocess.
    n_top_genes : int, optional
        Number of highly variable genes to select. Default is 2000.
    flavor : str, optional
        The method to use for selecting highly variable genes. Default is "seurat_v3".
    batch_key : str, optional
        The key in adata.obs to use for batch-aware HVG selection.
        Default is None (no batch-aware selection).
    covariates : Sequence[str], optional
        List of covariate keys in adata.obs to regress out.
        Default is None (no regression).
    layer : str, optional
        The layer in adata to use for HVG selection. Default is None (use adata.X).
    lognorm : bool, optional
        Whether to perform log-normalization. Default is True.
    hvg : bool, optional
        Whether to select highly variable genes. Default is True.
    regression : bool, optional
        Whether to regress out covariates. If None, automatically enabled when
        covariates are provided. Default is None.

    Returns
    -------
    AnnData
        The preprocessed annotated data matrix.

    """
    if regression is None:
        regression = covariates is not None
    if regression and not covariates:
        raise ValueError("regression requires at least one covariate")

    logger.info(
        "started preprocessing for {} cells and {} genes",
        adata.n_obs,
        adata.n_vars,
    )
    result = _lognorm(adata) if lognorm else adata.copy()
    if lognorm:
        logger.info("preprocessing step: log-normalization")
    else:
        logger.info("preprocessing step: log-normalization skipped")

    if hvg:
        logger.info(
            "preprocessing step: identifying {} highly variable genes",
            n_top_genes,
        )
        sc.pp.highly_variable_genes(
            result,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
            layer=layer,
        )
    else:
        logger.info("preprocessing step: highly variable gene selection skipped")

    if regression:
        logger.info("preprocessing step: regressing out {}", list(covariates))
        sc.pp.regress_out(result, keys=list(covariates))
    else:
        logger.info("preprocessing step: regression skipped")

    logger.info(
        "finished preprocessing for {} cells and {} genes",
        result.n_obs,
        result.n_vars,
    )
    return result
