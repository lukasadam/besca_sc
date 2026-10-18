"""Standard project paths for data, intermediate files, and results."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
EXTERNAL_DIR = DATA_DIR / "external"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

for directory in (
    RAW_DIR,
    EXTERNAL_DIR,
    INTERIM_DIR,
    PROCESSED_DIR,
    RESULTS_DIR,
):
    directory.mkdir(parents=True, exist_ok=True)
