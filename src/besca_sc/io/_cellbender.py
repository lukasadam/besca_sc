"""Load CellBender-filtered 10x data."""

from pathlib import Path

import anndata as ad
import scanpy as sc

from besca_sc.paths import RAW_DIR


def load_cellbender(
    library: str,
    *,
    raw_dir: Path = RAW_DIR,
    filename: str = "cellbender_filtered.h5",
) -> ad.AnnData:
    """Load a CellBender-filtered 10x library from ``data/raw``.

    Parameters
    ----------
    library
        Library subdirectory, for example ``"Infected"`` or ``"Uninfected"``.
    raw_dir
        Root raw-data directory. Defaults to the project ``data/raw`` folder.
    filename
        CellBender output filename.

    Returns
    -------
    AnnData
        Gene Expression counts with non-Gene Expression features removed.
    """
    path = raw_dir / library / filename
    if not path.exists():
        raise FileNotFoundError(f"CellBender file not found: {path}")

    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()

    if "feature_types" in adata.var:
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

    return adata
