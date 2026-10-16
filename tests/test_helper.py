"""Tests for AnnData helpers."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from besca_sc._helper import convert_ensembl_to_symbol, subset_adata


@pytest.fixture
def adata() -> AnnData:
    """Build a small AnnData fixture."""
    result = AnnData(
        X=np.arange(12).reshape(3, 4),
        obs=pd.DataFrame({"group": ["a", "b", "a"]}, index=["c", "a", "b"]),
        var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
    )
    result.raw = result.copy()
    return result


def test_subset_aligns_series_and_uses_raw(adata: AnnData) -> None:
    """Align a labeled mask and return the raw matrix."""
    mask = pd.Series([True, False, True], index=["b", "c", "a"])
    result = subset_adata(adata, mask)

    assert list(result.obs_names) == ["a", "b"]
    np.testing.assert_array_equal(result.X, [[4, 5, 6, 7], [8, 9, 10, 11]])


def test_subset_variables_and_non_raw(adata: AnnData) -> None:
    """Subset variables from the regular matrix."""
    result = subset_adata(adata, [True, False, True, False], use_raw=False, axis=1)
    assert list(result.var_names) == ["g1", "g3"]
    np.testing.assert_array_equal(result.X, [[0, 2], [4, 6], [8, 10]])


@pytest.mark.parametrize("criteria", [[True], [1, 0, 1], [True, None, False]])
def test_rejects_invalid_masks(adata: AnnData, criteria: list) -> None:
    """Reject malformed or non-boolean masks."""
    with pytest.raises((TypeError, ValueError)):
        subset_adata(adata, criteria, use_raw=False)


def test_raw_requires_raw_layer() -> None:
    """Require a raw layer when requested."""
    adata = AnnData(np.ones((2, 2)))
    with pytest.raises(ValueError, match="adata.raw"):
        subset_adata(adata, [True, False], use_raw=True)


def test_convert_ensembl_to_symbol(monkeypatch: pytest.MonkeyPatch) -> None:
    """Convert Ensembl IDs using the requested MyGene query parameters."""

    class FakeMyGeneInfo:
        def querymany(
            self, gene_list: list[str], **kwargs: object
        ) -> list[dict[str, str]]:
            return [{"symbol": "TP53"}, {"symbol": "BRCA1"}]

    fake_mygene = type("FakeMyGene", (), {"MyGeneInfo": FakeMyGeneInfo})
    monkeypatch.setitem(__import__("sys").modules, "mygene", fake_mygene)

    result = convert_ensembl_to_symbol(["ENSG00000141510", "ENSG00000012048"])

    assert result == ["TP53", "BRCA1"]
