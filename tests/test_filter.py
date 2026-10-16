"""Tests for quality-control filtering."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from besca_sc.pp._filter import filter


@pytest.fixture
def adata() -> AnnData:
    """Build a small AnnData object with deterministic QC metrics."""
    return AnnData(
        X=np.ones((3, 3)),
        obs=pd.DataFrame(
            {
                "n_genes_by_counts": [2, 4, 6],
                "total_counts": [5, 10, 15],
                "percent_mito": [0.01, 0.02, 0.20],
            },
            index=["cell1", "cell2", "cell3"],
        ),
        var=pd.DataFrame({"n_cells": [1, 2, 3]}, index=["g1", "g2", "g3"]),
    )


def test_filter_applies_cell_and_gene_thresholds(adata: AnnData) -> None:
    """Apply multiple thresholds in the documented order."""
    result = filter(adata, min_genes=3, max_counts=12, min_cells=2, max_mito=0.1)

    assert list(result.obs_names) == ["cell2"]
    assert list(result.var_names) == ["g2", "g3"]
    assert adata.n_obs == 3


def test_max_genes_filters_cells(adata: AnnData) -> None:
    """Remove cells above the maximum gene threshold."""
    result = filter(adata, max_genes=4)

    assert list(result.obs_names) == ["cell1", "cell2"]


def test_min_genes_is_inclusive(adata: AnnData) -> None:
    """Keep cells exactly at the minimum gene threshold."""
    result = filter(adata, min_genes=4)

    assert list(result.obs_names) == ["cell2", "cell3"]


def test_count_thresholds_use_aliases(adata: AnnData) -> None:
    """Support total_counts as the count metric alias."""
    result = filter(adata, min_counts=10, max_counts=10)

    assert list(result.obs_names) == ["cell2"]


def test_max_mito_is_exclusive(adata: AnnData) -> None:
    """Remove cells at or above the mitochondrial threshold."""
    result = filter(adata, max_mito=0.20)

    assert list(result.obs_names) == ["cell1", "cell2"]


def test_min_cells_filters_genes(adata: AnnData) -> None:
    """Remove genes expressed in fewer than the requested cells."""
    result = filter(adata, min_cells=2)

    assert list(result.var_names) == ["g2", "g3"]


def test_filter_does_not_modify_input(adata: AnnData) -> None:
    """Leave the input AnnData unchanged."""
    original_obs_names = list(adata.obs_names)
    original_var_names = list(adata.var_names)

    filter(adata, min_genes=4, min_cells=2)

    assert list(adata.obs_names) == original_obs_names
    assert list(adata.var_names) == original_var_names


def test_filter_rejects_conflicting_thresholds(adata: AnnData) -> None:
    """Reject minimum thresholds greater than their maximums."""
    with pytest.raises(ValueError, match="min_genes"):
        filter(adata, min_genes=5, max_genes=2)


def test_filter_rejects_negative_thresholds(adata: AnnData) -> None:
    """Reject negative quality-control thresholds."""
    with pytest.raises(ValueError, match="min_cells"):
        filter(adata, min_cells=-1)


def test_filter_reports_missing_requested_metric(adata: AnnData) -> None:
    """Raise a useful error when a requested metric is absent."""
    del adata.obs["n_genes_by_counts"]

    with pytest.raises(ValueError, match="max_genes"):
        filter(adata, max_genes=10)
