"""Tests for MAD threshold helpers."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from besca_sc.pp._outlier import calculate_thresholds, get_threshold


@pytest.fixture
def adata() -> ad.AnnData:
    """Build a small AnnData object with QC metrics."""
    return ad.AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {"score": [1.0, 1.0, 1.0, 1.0, 100.0]},
            index=[f"c{i}" for i in range(5)],
        ),
    )


def test_get_threshold_flags_upper_tail(adata: ad.AnnData) -> None:
    """Calculate an upper median/MAD threshold."""
    threshold = get_threshold(adata, True, "score", 2, "high")
    assert threshold == pytest.approx(1.0)


def test_get_threshold_applies_bounds(adata: ad.AnnData) -> None:
    """Clip calculated thresholds to supplied bounds."""
    threshold = get_threshold(adata, True, "score", 2, "high", upper_bound=0.5)
    assert threshold == pytest.approx(0.5)


def test_calculate_thresholds_uses_default_mads() -> None:
    """Use default per-metric MAD multipliers."""
    adata = ad.AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "log1p_total_counts": [1.0, 1.0, 1.0, 1.0, 2.0],
                "log1p_n_genes_by_counts": [1.0] * 5,
                "pct_counts_in_top_20_genes": [1.0] * 5,
                "logit_mt_fraction": [1.0] * 5,
            }
        ),
    )
    thresholds = calculate_thresholds(adata)
    assert set(thresholds) == {
        "log1p_total_counts",
        "log1p_n_genes_by_counts",
        "pct_counts_in_top_20_genes",
        "logit_mt_fraction",
    }


def test_get_threshold_validation(adata: ad.AnnData) -> None:
    """Reject invalid directions and missing metrics."""
    with pytest.raises(ValueError, match="direction"):
        get_threshold(adata, True, "score", 2, "invalid")
    with pytest.raises(ValueError, match="unknown metric"):
        get_threshold(adata, True, "missing", 2, "low")
