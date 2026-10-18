from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd


def subset_adata(
    adata: ad.AnnData,
    mask: pd.Series | Sequence[bool] | np.ndarray,
    *,
    axis: Literal[0, 1] = 0,
    use_raw: bool = False,
) -> ad.AnnData:
    """Return an independent subset of an AnnData object.

    Parameters
    ----------
    adata
        AnnData object to subset.
    mask
        Boolean mask, or a pandas Series indexed by observations or variables.
    use_raw
        Use ``adata.raw`` for the returned expression matrix and variables.
    axis
        ``0`` filters observations; ``1`` filters variables.

    Returns
    -------
    AnnData
        A copied AnnData subset.

    Examples
    --------
    >>> import besca_sc as bsc
    >>> adata = bsc.datasets.simulated_pbmc3k_processed()
    >>> mask = adata.obs["dblabel"].eq("naive B cell")
    >>> adata_BCELL = subset_adata(adata, mask)

    """
    if axis not in (0, 1):
        raise ValueError("axis must be 0 (observations) or 1 (variables)")

    axis_names = adata.obs_names if axis == 0 else adata.var_names
    if isinstance(mask, pd.Series):
        series = mask
        if not series.index.is_unique or set(series.index) != set(axis_names):
            raise ValueError(
                "filter criteria index must match the selected AnnData names"
            )
        mask = series.reindex(axis_names).to_numpy()
    else:
        mask = np.asarray(mask)

    if mask.ndim != 1 or len(mask) != len(axis_names):
        raise ValueError(
            "filter criteria must be a one-dimensional mask "
            f"of length {len(axis_names)}"
        )
    if not np.issubdtype(mask.dtype, np.bool_):
        raise TypeError("filter criteria must contain boolean values")

    subset_full = adata[mask, :].copy() if axis == 0 else adata[:, mask].copy()
    if not use_raw:
        return subset_full
    if subset_full.raw is None:
        raise ValueError("use_raw=True requires adata.raw to be set")
    return subset_full.raw.to_adata()
