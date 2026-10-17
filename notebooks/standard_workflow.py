"""Interactive standard workflow for single-cell data."""

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys
    from pathlib import Path

    # Keep the notebook runnable when Marimo is launched outside Pixi.
    project_src = Path(__file__).resolve().parents[1] / "src"
    if str(project_src) not in sys.path:
        sys.path.insert(0, str(project_src))

    import os
    import marimo as mo
    import scanpy as sc
    import anndata as ad
    import pandas as pd
    import numpy as np

    sc.settings.verbosity = 2  # show INFO messages from scanpy
    sc.settings.set_figure_params(dpi=100, facecolor="white")

    print(f"scanpy  version: {sc.__version__}")
    print(f"anndata version: {ad.__version__}")
    print(f"pandas  version: {pd.__version__}")
    print(f"numpy   version: {np.__version__}")

    from besca_sc import filter_cells

    return mo, os


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Step 1 — Define file paths
    """)
    return


@app.cell
def _(os):
    BASE = "/Users/adaml9/besca_sc/data"

    # CellBender ambient-corrected, filtered count matrices
    CB_INF = os.path.join(BASE, "raw/Infected/cellbender_filtered.h5")
    CB_UNINF = os.path.join(BASE, "raw/Uninfected/cellbender_filtered.h5")

    # HADGE consensus demultiplexing results (7-tool majority vote)
    # Columns: Barcode, [7 tool assignments], Consensus_Assignment, Confidence,
    #          Consensus_DropletType, Donor, Condition
    DEMUX_INF = os.path.join(BASE, "raw/Infected/Infected.consensus.csv")
    DEMUX_UNINF = os.path.join(BASE, "raw/Uninfected/Uninfected.consensus.csv")

    # Output directory
    OUT_DIR = os.path.join(BASE, "scanpy_analysis")
    os.makedirs(OUT_DIR, exist_ok=True)

    # Verify all input files exist before proceeding
    for f in [CB_INF, CB_UNINF, DEMUX_INF, DEMUX_UNINF]:
        status = "OK" if os.path.exists(f) else "MISSING!"
        print(f"  [{status}]  {f}")
    return


if __name__ == "__main__":
    app.run()
