"""Outlier detection helpers for AnnData objects."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
from loguru import logger
from scipy.special import logit
from scipy.stats import median_abs_deviation


def mad_thresholds(values: pd.Series, nmads: float) -> tuple[float, float]:
    """Return lower and upper thresholds based on the median absolute deviation."""
    if nmads < 0:
        raise ValueError("nmads must be non-negative")
    if not pd.api.types.is_numeric_dtype(values):
        raise TypeError("metric must contain numeric values")

    median = float(np.nanmedian(values))
    mad = float(median_abs_deviation(values, nan_policy="omit"))
    return median - nmads * mad, median + nmads * mad


def is_outlier(
    adata: ad.AnnData,
    metric: str,
    nmads: float,
    side: Literal["both", "upper"] = "both",
    groupby: str | None = None,
) -> pd.Series:
    """Return MAD-based outlier flags, globally or independently per group."""
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")
    if metric not in adata.obs:
        raise ValueError(f"metric must be present in adata.obs: {metric!r}")
    if side not in ("both", "upper"):
        raise ValueError("side must be 'both' or 'upper'")
    if groupby is not None and groupby not in adata.obs:
        raise ValueError(f"groupby must be present in adata.obs: {groupby!r}")

    values = adata.obs[metric]
    if groupby is None:
        groups = [(None, values.index)]
    else:
        groups = adata.obs.groupby(groupby, dropna=False).groups.items()

    flags = pd.Series(False, index=values.index, dtype=bool)
    for _, indices in groups:
        group_values = values.loc[indices]
        lower, upper = mad_thresholds(group_values, nmads)
        if side == "upper":
            group_flags = group_values > upper
        else:
            group_flags = (group_values < lower) | (group_values > upper)
        flags.loc[indices] = group_flags
    return flags.fillna(False).astype(bool)


def flag_outliers(
    adata: ad.AnnData,
    qc_filters: Mapping[str, tuple[float, Literal["both", "upper"]]] | None = None,
    mt_threshold: float = 8.0,
    groupby: str | None = None,
) -> ad.AnnData:
    """Add combined QC outlier flags to an AnnData object.

    The mitochondrial flag is restricted to cells whose ``pct_counts_mt`` is
    above ``mt_threshold``. By default, the function uses common Scanpy QC
    metrics; callers can provide a custom mapping of metric names to
    ``(nmads, side)`` pairs.

    Parameters
    ----------
    adata
        AnnData object containing the QC metrics in ``adata.obs``.
    qc_filters
        Mapping of metric names to MAD settings. The special metric
        ``logit_mt_fraction`` is used to create ``mt_outlier``.
    mt_threshold
        Minimum ``pct_counts_mt`` value required for a mitochondrial outlier.
    groupby
        Optional observation column used to calculate MAD thresholds
        independently for each sample or biological group.

    Returns
    -------
    AnnData
        The input ``adata`` with ``outlier`` and ``mt_outlier`` columns added.
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")
    if mt_threshold < 0:
        raise ValueError("mt_threshold must be non-negative")

    logger.info(
        "started outlier detection for {} cells and {} genes",
        adata.n_obs,
        adata.n_vars,
    )

    filters = dict(
        qc_filters
        or {
            "log1p_total_counts": (5, "both"),
            "log1p_n_genes_by_counts": (5, "both"),
            "pct_counts_in_top_20_genes": (5, "both"),
            "logit_mt_fraction": (3, "upper"),
        }
    )
    filters.setdefault("logit_mt_fraction", (3, "upper"))

    required = {"total_counts_mt", "total_counts", "pct_counts_mt"}
    missing = required - set(adata.obs.columns)
    if missing:
        raise ValueError(f"missing mitochondrial QC metrics: {sorted(missing)}")

    mt_fraction = (adata.obs["total_counts_mt"] + 1) / (adata.obs["total_counts"] + 2)
    adata.obs["logit_mt_fraction"] = logit(mt_fraction)

    flags = {
        metric: is_outlier(adata, metric, nmads, side, groupby=groupby)
        for metric, (nmads, side) in filters.items()
    }
    for metric, metric_flags in flags.items():
        scope = f" across groups in '{groupby}'" if groupby else ""
        logger.info(
            "metric '{}' flagged {} outliers{}",
            metric,
            int(metric_flags.sum()),
            scope,
        )
    non_mito_flags = [
        flags[metric] for metric in filters if metric != "logit_mt_fraction"
    ]
    adata.obs["outlier"] = (
        np.logical_or.reduce(non_mito_flags)
        if non_mito_flags
        else np.zeros(adata.n_obs, dtype=bool)
    )
    adata.obs["mt_outlier"] = flags["logit_mt_fraction"] & (
        adata.obs["pct_counts_mt"] > mt_threshold
    )
    logger.info(
        "finished outlier detection: {} outliers and {} mitochondrial outliers",
        int(adata.obs["outlier"].sum()),
        int(adata.obs["mt_outlier"].sum()),
    )
    return adata
