# Pixi Cheat Sheet

Quick reference for managing Python projects, dependencies, environments, and development workflows with [Pixi](https://pixi.sh/).

## Project initialization

```bash
pixi init                                # Create a project with pixi.toml
pixi init --format pyproject myproject   # Initialize a Python project
pixi info                                # Show workspace and environment information
pixi --version                           # Check Pixi version
```

## Dependency management

```bash
pixi add python=3.12                     # Add conda dependency
pixi add scanpy anndata                  # Add multiple conda dependencies
pixi add --pypi scdrs                    # Add package from PyPI
pixi add --feature dev pytest            # Add to dev feature
pixi add --feature dev --pypi pre-commit # PyPI dependency in dev feature

pixi remove scanpy                       # Remove conda dependency
pixi remove --pypi scdrs                 # Remove PyPI dependency
pixi upgrade scanpy                      # Upgrade dependency
pixi upgrade                             # Upgrade all dependencies
pixi list -e dev                         # List installed packages
```

Prefer conda-forge for scientific Python packages when available. Use PyPI when necessary. Avoid declaring conflicting versions of the same package in conda and PyPI.

## Features and environments

A **feature** groups dependencies and tasks. An **environment** combines features into something installable.

```toml
[feature.dev.dependencies]
pytest = "*"
ruff = "*"

[feature.dev.pypi-dependencies]
pre-commit = "*"

[environments]
dev = ["dev"]
```

Equivalent CLI:

```bash
pixi workspace environment add dev --feature dev
```

Run commands in an environment:

```bash
pixi install -e dev
pixi run -e dev python
pixi run -e dev pytest
pixi shell -e dev
```

Features do not automatically create environments. Top-level dependencies are included in each environment by default.

## Running commands and tasks

```bash
pixi run python script.py          # Run in default environment
pixi run -e dev pytest             # Run in dev environment
pixi shell                        # Enter default environment
pixi shell -e dev                 # Enter dev environment
pixi task list                    # List registered tasks
```

Define reusable tasks in `pyproject.toml`:

```toml
[tool.pixi.tasks]
test = "pytest tests/"
lint = "ruff check src/ tests/"
format = "ruff format src/ tests/"
```

Or in `pixi.toml`, use `[tasks]` instead.

```bash
pixi task add test "pytest tests/"
pixi run test
```

## Python package development

Recommended structure:

```text
myproject/
├── pyproject.toml
├── pixi.lock
├── src/
│   └── myproject/
│       ├── __init__.py
│       └── qc.py
└── tests/
    └── test_qc.py
```

For projects using Hatchling, a typical local editable dependency is:

```toml
[tool.pixi.pypi-dependencies]
myproject = { path = ".", editable = true }
```

This allows Python modules to be imported directly from `src/` and edited without reinstalling after every change.

```bash
pixi run python -c "import myproject"
pixi run pytest
```

## Reproducibility

```bash
pixi install                 # Install and update lockfile if needed
pixi install --locked        # Fail if lockfile is outdated
pixi install --frozen        # Use existing lockfile without updating
pixi lock                    # Resolve and generate lockfile
pixi clean                   # Remove workspace environments
```

Commit `pixi.lock` to Git for reproducibility. Usually exclude `.pixi/` from version control.

## Troubleshooting

```bash
pixi info                    # Inspect workspace configuration
pixi list -e dev             # Inspect resolved packages
pixi run -e dev python -V    # Check active Python version
pixi run -e dev which python # Locate environment Python (macOS/Linux)
pixi --help                  # Main command reference
pixi add --help              # Command-specific help
```

Common problems:

- **Unused feature:** Add it to an environment under `[environments]`.
- **Unknown environment:** Define the environment before running `pixi install -e`.
- **Conda/PyPI conflict:** Check whether a package is already pinned by the conda solver.
- **Import fails:** Check whether the local Python package is installed in editable mode.
- **Unexpected dependency versions:** Inspect `pixi.lock` and the resolved environment.

## Useful references

- [Pixi documentation](https://pixi.sh/)
- [CLI reference](https://pixi.sh/latest/reference/cli/pixi/)
- [Multi-environment tutorial](https://pixi.sh/latest/tutorials/multi_environment/)
- [Manifest reference](https://pixi.sh/latest/reference/pixi_manifest/)