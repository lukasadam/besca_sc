"""Outlier detection helpers for AnnData objects."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

Side = Literal["both", "lower", "upper"]
FilterSpec = tuple[float, Side]

# Default MAD-based QC filters based on single-cell best-practice guidelines:
# See https://www.sc-best-practices.org/preprocessing-visualization/quality-control.
default_nmads: dict[str, FilterSpec] = {
    "log1p_total_counts": (5, "both"),
    "log1p_n_genes_by_counts": (5, "both"),
    "pct_counts_in_top_20_genes": (5, "both"),
    "logit_mt_fraction": (3, "upper"),
}

def get_threshold(
    adata: ad.AnnData,
    cell: bool,
    metric: str,
    nmads: float,
    direction: str,
    lower_bound: float | None = None,
    upper_bound: float | None = None,
) -> float:
    """Calculate a bounded MAD threshold for an observation or gene metric."""
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")
    if direction not in {"low", "high"}:
        raise ValueError("direction must be 'low' or 'high'")
    if not np.isfinite(nmads) or nmads < 0:
        raise ValueError("nmads must be finite and non-negative")
    table = adata.obs if cell else adata.var
    if metric not in table:
        scope = "adata.obs" if cell else "adata.var"
        raise ValueError(f"unknown metric {metric!r} in {scope}")
    values = pd.to_numeric(table[metric], errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError(f"metric {metric!r} contains no finite values")
    median = float(np.median(values))
    mad = float(median_abs_deviation(values))
    estimate = median + nmads * mad if direction == "high" else median - nmads * mad
    if lower_bound is not None:
        estimate = max(estimate, lower_bound)
    if upper_bound is not None:
        estimate = min(estimate, upper_bound)
    return float(estimate)


def calculate_thresholds(
    adata: ad.AnnData,
    metric_bounds: Mapping[str, tuple[float | None, float | None]] | None = None,
    *,
    nmads: float | Mapping[str, float] | None = None,
    cell: bool = True,
) -> dict[str, tuple[float, float]]:
    """Calculate low/high MAD thresholds for QC metrics.

    If ``metric_bounds`` is omitted, the metrics and MAD multipliers from
    ``DEFAULT_FILTERS`` are used. A scalar ``nmads`` overrides all defaults;
    a mapping overrides selected metrics.
    """
    if metric_bounds is None:
        metric_bounds = {metric: (None, None) for metric in default_nmads}

    def metric_nmads(metric: str) -> float:
        if isinstance(nmads, Mapping):
            return float(nmads.get(metric, default_nmads.get(metric, (5, "both"))[0]))
        if nmads is not None:
            return float(nmads)
        return float(default_nmads.get(metric, (5, "both"))[0])

    return {
        metric: (
            get_threshold(
                adata, cell, metric, metric_nmads(metric), "low", lower, upper
            ),
            get_threshold(
                adata, cell, metric, metric_nmads(metric), "high", lower, upper
            ),
        )
        for metric, (lower, upper) in metric_bounds.items()
    }
