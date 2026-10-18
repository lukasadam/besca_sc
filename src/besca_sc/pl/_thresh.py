"""Plot QC metrics with MAD-based threshold lines."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc

Threshold = tuple[float, float | None]


def plot_thresholds(
    adata: ad.AnnData,
    thresholds: Mapping[str, Threshold] | None = None,
    *,
    metrics: Sequence[str] | None = None,
    figsize: tuple[float, float] | None = None,
    axes: Sequence[plt.Axes] | None = None,
) -> tuple[plt.Figure, list[plt.Axes]]:
    """Plot QC distributions and optional MAD thresholds.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing QC metrics in ``adata.obs``.
    thresholds : mapping, optional
        Maps each metric to ``(lower, upper)``. Use ``None`` for a missing
        bound, e.g. ``{"pct_counts_mt": (None, 8.0)}``.
    metrics : sequence of str, optional
        Columns from ``adata.obs`` to plot.
    figsize : tuple of float, optional
        Figure size. Defaults to ``(5 * len(metrics), 4)``.
    axes : sequence of Axes, optional
        Existing axes to draw on.

    Returns
    -------
    figure, axes
        The Matplotlib figure and axes used for the plots.
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")
    if metrics is None:
        metrics = tuple(thresholds) if thresholds is not None else (
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        )
    else:
        metrics = tuple(metrics)
    missing = [metric for metric in metrics if metric not in adata.obs]
    if missing:
        raise ValueError(f"missing QC metrics in adata.obs: {missing}")

    thresholds = {} if thresholds is None else thresholds
    if axes is not None:
        raise ValueError(
            "axes cannot be supplied when using Scanpy multi_panel violins"
        )

    sc.pl.violin(
        adata,
        keys=list(metrics),
        rotation=45,
        jitter=0.15,
        stripplot=True,
        multi_panel=True,
        show=False,
    )
    fig = plt.gcf()
    if figsize is not None:
        fig.set_size_inches(figsize)
    plot_axes = list(fig.axes)
    for axis, metric in zip(plot_axes, metrics):
        lower, upper = thresholds.get(metric, (None, None))
        if lower is not None:
            axis.axhline(
                lower, color="tab:orange", linestyle="--", label="lower threshold"
            )
        if upper is not None:
            axis.axhline(
                upper, color="tab:red", linestyle="--", label="upper threshold"
            )
        if lower is not None or upper is not None:
            axis.legend(loc="best")

    return fig, plot_axes
