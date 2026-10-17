# Marimo Cheat Sheet

Quick reference for running and developing reactive Python notebooks in this project with [Marimo](https://marimo.io/).

## Start the workflow notebook

Run the notebook from the `besca` Pixi environment:

```bash
pixi run -e besca marimo edit notebooks/standard_workflow.py
```

This opens the notebook editor in your browser. The notebook is stored as a regular Python file, so it can be reviewed, formatted, and committed like normal source code.

Run it without opening the editor:

```bash
pixi run -e besca marimo run notebooks/standard_workflow.py
```

## Notebook structure

A Marimo notebook is a Python application made of cells:

```python
import marimo

app = marimo.App()


@app.cell
def _():
    import scanpy as sc

    return sc


if __name__ == "__main__":
    app.run()
```

Cells are reactive: when a value changes, Marimo reruns cells that depend on it and keeps the dependency graph consistent.

## Use the local BESCA package

The project is installed in editable mode through Pixi:

```toml
[tool.pixi.pypi-dependencies]
besca_sc = { path = ".", editable = true }
```

This allows notebook code to import the package directly:

```python
from besca_sc import filter_cells, flag_outliers
```

Changes made under `src/besca_sc/` are available after rerunning the affected cell; reinstalling the package is normally unnecessary.

## Load an AnnData file

The workflow notebook accepts an `.h5ad` path through a text input. A minimal equivalent is:

```python
adata = sc.read_h5ad("data/input.h5ad")
adata
```

Keep large datasets outside Git and add their directory to `.gitignore` where appropriate.

## Add interactive controls

Marimo UI elements return values that can be used by later cells:

```python
min_genes = mo.ui.number(value=500, start=0, label="Minimum genes")
min_genes
```

Use the value in another cell:

```python
filtered = filter_cells(adata, min_genes=min_genes.value)
filtered
```

## Validate and format

Compile-check the notebook:

```bash
pixi run -e besca python -m py_compile notebooks/standard_workflow.py
```

Run the project checks:

```bash
pixi run -e dev lint
pixi run -e dev format-check
pixi run -e dev test
```

Format the notebook source:

```bash
pixi run -e dev format
```

## Export and convert

Export a notebook as HTML:

```bash
pixi run -e besca marimo export html notebooks/standard_workflow.py -o notebooks/standard_workflow.html
```

Convert an existing Jupyter notebook to Marimo format:

```bash
pixi run -e besca marimo convert notebooks/example.ipynb -o notebooks/example.py
```

## Troubleshooting

- **`marimo: command not found`:** use `pixi run -e besca marimo ...` and ensure the `besca` environment includes the `notebook` feature.
- **`No module named besca_sc`:** install the environment and confirm the editable dependency is present in `pyproject.toml` or `pixi.toml`.
- **Input file not found:** paths are resolved from the directory where Marimo is launched, usually the repository root.
- **Pixi certificate errors:** use the system certificate store with `PIXI_TLS_ROOT_CERTS=system`.
- **Notebook changes are not reflected:** rerun the affected cell or restart the Marimo session after changing package imports.

## Useful references

- [Marimo documentation](https://docs.marimo.io/)
- [Marimo editor guide](https://docs.marimo.io/guides/editor_features/)
- [Pixi cheatsheet](PIXI_CHEATSHEET.md)
