from __future__ import annotations

import numpy as np
from scipy import sparse


def is_count_matrix(X, atol: float = 1e-8) -> bool:
    """Check whether a matrix is consistent with non-negative counts."""
    values = X.data if sparse.issparse(X) else np.asarray(X).ravel()

    return bool(
        np.all(np.isfinite(values))
        and np.all(values >= 0)
        and np.all(
            np.isclose(values, np.rint(values), atol=atol, rtol=0)
        )
    )

def convert_ensembl_to_symbol(gene_list: list, species: str = "human") -> list:
    """Convert Ensembl gene IDs to symbols.

    Uses the ``mygene`` package to
    look up the supplied list of ENSEMBLE Ids and return
    the equivalent list of Symbols. Species needs to be supplied.

    Parameters
    ----------
    gene_list: `list`
        list of ensemble_Ids that need to be converted
    species: `str` | default = 'human'
        string identifying the species of the supplied Ensemble IDs

    Returns
    -------
    list
        List containing the converted Symbols

    """
    from mygene import MyGeneInfo

    mg = MyGeneInfo()
    gene_symbols = mg.querymany(
        gene_list, scopes="ensembl.gene", fields="symbol", species=species
    )

    symbols = []
    for x in gene_symbols:
        symbols.append(x.get("symbol"))

    return symbols
