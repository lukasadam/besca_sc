"""Smoke tests for the package scaffold."""


def test_package_imports() -> None:
    """The package can be imported before optional modules are added."""
    import besca_sc

    assert "filter_cells_and_genes" in besca_sc.__all__
    assert "calculate_thresholds" in besca_sc.__all__
