"""Smoke tests for the package scaffold."""


def test_package_imports() -> None:
    """The package can be imported before optional modules are added."""
    import besca_sc

    assert besca_sc.__all__ == ["filter_cells", "flag_outliers"]
