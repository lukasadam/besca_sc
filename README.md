# besca-sc

Utilities and notebooks for single-cell RNA-seq quality control, preprocessing,
doublet detection, and clustering.

## Requirements

- [Pixi](https://pixi.sh/latest/)
- macOS Apple Silicon or Linux x86-64

Install Pixi:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Restart your shell, then verify the installation:

```bash
pixi --version
```

## Install the project

From the repository root:

```bash
pixi install -e besca
```

The `besca` environment includes the package, Scanpy, scikit-misc for Seurat-v3
highly variable gene selection, and Leiden clustering (`leidenalg` and
`python-igraph`).

To inspect the environment:

```bash
pixi list -e besca
```

## Register the Jupyter kernel

Install the project kernel so it is selectable in Jupyter or VS Code:

```bash
pixi run -e besca python -m ipykernel install --user \
  --name besca-sc \
  --display-name "Python (besca-sc)"
```

Start Jupyter from the same environment:

```bash
pixi run -e besca jupyter lab
```

In VS Code, select the kernel named **Python (besca-sc)** before running a
notebook. Restart the kernel after changing Pixi dependencies.

## Run the notebooks

Open the notebooks in this order:

1. `notebooks/01_load_and_merge.ipynb`
2. `notebooks/02_initial_qc.ipynb`
3. `notebooks/03_doublet_detection.ipynb`
4. `notebooks/04_final_qc.ipynb`
5. `notebooks/05_pre_integration_clustering.ipynb`

The notebooks expect project data under `data/`. In particular, raw Matrix
Market input is expected under `data/processed/raw/`, together with the
corresponding gene, barcode, and metadata files.

## Development checks

Install the development environment and run the checks with:

```bash
pixi install -e dev
pixi run -e dev ruff check src tests
pixi run -e dev pytest -q
```

If a local pytest plugin causes collection errors, disable third-party plugin
autoloading:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 pixi run -e dev pytest -q
```

## Useful commands

```bash
pixi shell -e besca       # interactive project shell
pixi run -e besca python  # Python in the project environment
pixi update               # update the lockfile
pixi install              # synchronize environments with pixi.lock
```
