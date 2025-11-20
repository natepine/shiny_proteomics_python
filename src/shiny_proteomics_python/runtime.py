from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import os
import shutil
import logging

from shiny_proteomics_python.cache import Cache
from shiny_proteomics_python.config import AppConfig, load_config


logger = logging.getLogger(__name__)

APP_DIR = Path(os.environ.get("APP_DIR", Path(__file__).resolve().parents[1]))
DATA_DIR = APP_DIR / "data"
CACHE_DIR = DATA_DIR / "cache"
DB_PATH = DATA_DIR / "database.db"
CONFIG_ENV_VAR = "CONFIG_PATH"


@lru_cache(maxsize=1)
def get_config() -> AppConfig:
    configured_path = os.environ.get(CONFIG_ENV_VAR)
    config_path = Path(configured_path) if configured_path else DATA_DIR / "conf.yml"

    if not config_path.is_absolute():
        config_path = (APP_DIR / config_path).resolve()

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not config_path.exists():
        config_path = _bootstrap_config(config_path)

    return load_config(config_path)


def _bootstrap_config(target_path: Path) -> Path:
    template_path = APP_DIR.parent / "conf.tmp"

    if template_path.exists():
        shutil.copy(template_path, target_path)
        logger.warning(
            "Configuration file was missing; copied template from %s to %s. "
            "Please edit the new file with your MassPike credentials and species list before using the app.",
            template_path,
            target_path,
        )
        return target_path

    raise FileNotFoundError(
        f"Missing configuration file: {target_path}. Provide one or copy conf.tmp to this path."
    )


@lru_cache(maxsize=1)
def get_cache() -> Cache:
    return Cache(CACHE_DIR)


def get_database_path() -> Path:
    return DB_PATH
