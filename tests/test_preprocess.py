"""Tests for the preprocessing pipeline."""

import anndata as ad
import numpy as np
import pytest

from besca_sc.pp import _preprocess


@pytest.fixture
def adata() -> ad.AnnData:
    """Build a small expression matrix for pipeline tests."""
    return ad.AnnData(np.ones((2, 3)))


def test_preprocess_runs_enabled_steps_in_order(
    adata: ad.AnnData, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run normalization, HVG selection, and regression in sequence."""
    calls: list[tuple[str, object]] = []

    def fake_lognorm(value: ad.AnnData) -> ad.AnnData:
        calls.append(("lognorm", None))
        return value.copy()

    def fake_hvg(value: ad.AnnData, **kwargs: object) -> None:
        calls.append(("hvg", kwargs))

    def fake_regress(value: ad.AnnData, **kwargs: object) -> None:
        calls.append(("regression", kwargs))

    monkeypatch.setattr(_preprocess, "_lognorm", fake_lognorm)
    monkeypatch.setattr(_preprocess.sc.pp, "highly_variable_genes", fake_hvg)
    monkeypatch.setattr(_preprocess.sc.pp, "regress_out", fake_regress)

    result = _preprocess.preprocess(
        adata,
        n_top_genes=25,
        flavor="cell_ranger",
        batch_key="sample",
        covariates=["total_counts"],
        layer="counts",
    )

    assert result is not adata
    assert [name for name, _ in calls] == ["lognorm", "hvg", "regression"]
    assert calls[1][1] == {
        "n_top_genes": 25,
        "flavor": "cell_ranger",
        "batch_key": "sample",
        "layer": "counts",
    }
    assert calls[2][1] == {"keys": ["total_counts"]}


def test_preprocess_can_disable_each_step(
    adata: ad.AnnData, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Skip all processing steps when their switches are disabled."""
    calls: list[str] = []
    monkeypatch.setattr(
        _preprocess, "_lognorm", lambda value: calls.append("lognorm") or value
    )
    monkeypatch.setattr(
        _preprocess.sc.pp,
        "highly_variable_genes",
        lambda *args, **kwargs: calls.append("hvg"),
    )
    monkeypatch.setattr(
        _preprocess.sc.pp,
        "regress_out",
        lambda *args, **kwargs: calls.append("regression"),
    )

    result = _preprocess.preprocess(
        adata, lognorm=False, hvg=False, regression=False, covariates=["batch"]
    )

    assert result is not adata
    assert calls == []


def test_preprocess_requires_covariates_for_regression(adata: ad.AnnData) -> None:
    """Reject an explicitly enabled regression step without covariates."""
    with pytest.raises(ValueError, match="regression requires"):
        _preprocess.preprocess(adata, lognorm=False, hvg=False, regression=True)
