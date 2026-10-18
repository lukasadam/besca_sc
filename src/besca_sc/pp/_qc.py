"""Quality-control metric calculation for AnnData objects."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import anndata as ad
import scanpy as sc
from loguru import logger
from scipy.special import logit


def qc(
    adata: ad.AnnData,
    *,
    species: Literal["human", "mouse"] = "human",
    mitochondrial_prefixes: Sequence[str] | None = None,
    ribosomal_prefixes: Sequence[str] | None = None,
    percent_top: Sequence[int] | None = (20, 50, 100),
    log1p: bool = True,
) -> ad.AnnData:
    """Annotate genes and calculate cell-level QC metrics.

    Mitochondrial and ribosomal genes are identified using species-specific
    or caller-provided prefixes. Scanpy adds the resulting metrics to
    ``adata.obs`` and the gene masks to ``adata.var`` in place.

    Parameters
    ----------
    adata : AnnData
        Annotated count matrix to update.
    species : {"human", "mouse"}, optional
        Species convention used for gene prefixes. Defaults to ``"human"``.
        Human uses ``MT-``, ``RPS`` and ``RPL``; mouse uses ``mt-``, ``Rps``
        and ``Rpl``.
    mitochondrial_prefixes : sequence of str, optional
        Override the species-specific mitochondrial prefixes.
    ribosomal_prefixes : sequence of str, optional
        Override the species-specific ribosomal prefixes.
    percent_top : sequence of int or None, optional
        Ranks for top-gene percentage metrics. Defaults to ``(20, 50, 100)``.
        Set to ``None`` to skip these metrics.
    log1p : bool, optional
        Whether to add log-transformed QC metrics. Defaults to ``True``.

    Returns
    -------
    AnnData
        The input object with QC annotations added.
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")

    prefixes = {
        "human": (("MT-",), ("RPS", "RPL")),
        "mouse": (("mt-",), ("Rps", "Rpl")),
    }
    if species not in prefixes:
        raise ValueError("species must be 'human' or 'mouse'")
    default_mt, default_ribo = prefixes[species]
    mt_prefixes = tuple(
        default_mt if mitochondrial_prefixes is None else mitochondrial_prefixes
    )
    ribo_prefixes = tuple(
        default_ribo if ribosomal_prefixes is None else ribosomal_prefixes
    )
    if not mt_prefixes or not ribo_prefixes:
        raise ValueError("gene prefix collections must not be empty")

    adata.var["ribo"] = adata.var_names.str.startswith(ribo_prefixes)
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefixes)

    mt_count = int(adata.var["mt"].sum())
    ribo_count = int(adata.var["ribo"].sum())
    if mt_count == 0:
        logger.warning("no mitochondrial genes detected using prefixes {}", mt_prefixes)
    if ribo_count == 0:
        logger.warning("no ribosomal genes detected using prefixes {}", ribo_prefixes)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=percent_top,
        log1p=log1p,
        inplace=True,
    )
    # Use a smoothed mitochondrial fraction so it remains finite for cells
    # with zero mitochondrial counts or zero total counts.
    adata.obs["logit_mt_fraction"] = logit(
        (adata.obs["total_counts_mt"] + 1)
        / (adata.obs["total_counts"] + 2)
    )
    logger.info(
        "calculated QC metrics for {} cells and {} genes "
        "({} mitochondrial, {} ribosomal); added logit_mt_fraction",
        adata.n_obs,
        adata.n_vars,
        mt_count,
        ribo_count,
    )
    return adata
