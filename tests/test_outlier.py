"""Tests for outlier detection helpers."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from besca_sc.pp._outlier import flag_outliers, is_outlier, mad_thresholds


@pytest.fixture
def adata() -> ad.AnnData:
    """Build an AnnData object with QC metrics and known outliers."""
    return ad.AnnData(
        X=np.ones((6, 1)),
        obs=pd.DataFrame(
            {
                "score": [1.0, 1.0, 1.0, 1.0, 1.0, 100.0],
                "sample": ["a", "a", "a", "b", "b", "b"],
                "total_counts_mt": [1, 1, 1, 1, 1, 100],
                "total_counts": [100, 100, 100, 100, 100, 100],
                "pct_counts_mt": [1, 1, 1, 1, 1, 20],
            },
            index=[f"cell{i}" for i in range(6)],
        ),
    )


def test_mad_thresholds() -> None:
    """Calculate symmetric thresholds around the median."""
    lower, upper = mad_thresholds(pd.Series([1.0, 2.0, 3.0]), nmads=2)

    assert lower == pytest.approx(0.0)
    assert upper == pytest.approx(4.0)


def test_is_outlier_flags_upper_values(adata: ad.AnnData) -> None:
    """Flag values above the upper MAD threshold."""
    flags = is_outlier(adata, "score", nmads=2, side="upper")

    assert flags.index.equals(adata.obs.index)
    assert flags.tolist() == [False, False, False, False, False, True]


def test_is_outlier_calculates_thresholds_per_group(adata: ad.AnnData) -> None:
    """Calculate independent thresholds for each sample."""
    flags = is_outlier(adata, "score", nmads=2, groupby="sample")

    assert flags.tolist() == [False, False, False, False, False, True]


def test_flag_outliers_adds_combined_flags(adata: ad.AnnData) -> None:
    """Add combined general and mitochondrial flags to obs."""
    result = flag_outliers(
        adata,
        qc_filters={"score": (2, "upper"), "logit_mt_fraction": (2, "upper")},
        mt_threshold=8,
        groupby="sample",
    )

    assert result is adata
    assert result.obs["outlier"].tolist() == [False, False, False, False, False, True]
    assert result.obs["mt_outlier"].tolist() == [
        False,
        False,
        False,
        False,
        False,
        True,
    ]
    assert "logit_mt_fraction" in result.obs


def test_flag_outliers_requires_mitochondrial_metrics(adata: ad.AnnData) -> None:
    """Reject high-level flagging when mitochondrial metrics are missing."""
    del adata.obs["pct_counts_mt"]

    with pytest.raises(ValueError, match="mitochondrial"):
        flag_outliers(adata, qc_filters={"score": (2, "upper")})


def test_outlier_validation(adata: ad.AnnData) -> None:
    """Reject invalid metrics, sides, groups, and MAD counts."""
    with pytest.raises(ValueError, match="metric"):
        is_outlier(adata, "missing", nmads=2)
    with pytest.raises(ValueError, match="side"):
        is_outlier(adata, "score", nmads=2, side="lower")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="groupby"):
        is_outlier(adata, "score", nmads=2, groupby="missing")
    with pytest.raises(ValueError, match="nmads"):
        mad_thresholds(adata.obs["score"], -1)
