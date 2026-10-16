"""Tests for log-normalization helpers."""

import numpy as np
from anndata import AnnData

from besca_sc.pp._lognorm import lognorm


def test_lognorm_normalizes_and_preserves_counts() -> None:
    """Normalize counts while retaining the original matrix and layers."""
    counts = np.array([[1.0, 3.0], [2.0, 2.0]])
    adata = AnnData(counts.copy())

    result = lognorm(adata, target_sum=10.0)

    expected = np.log1p(np.array([[2.5, 7.5], [5.0, 5.0]]))
    np.testing.assert_allclose(result.X, expected)
    np.testing.assert_allclose(result.layers["lognorm"], expected)
    np.testing.assert_array_equal(result.layers["counts"], counts)
    np.testing.assert_array_equal(result.raw.X, counts)
    np.testing.assert_array_equal(adata.X, counts)
    assert result is not adata

