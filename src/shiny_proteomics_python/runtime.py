from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import os

from shiny_proteomics_python.cache import Cache
from shiny_proteomics_python.config import AppConfig, load_config


APP_DIR = Path(os.environ.get("APP_DIR", Path(__file__).resolve().parents[1]))
DATA_DIR = APP_DIR / "data"
CACHE_DIR = DATA_DIR / "cache"
DB_PATH = DATA_DIR / "database.db"


@lru_cache(maxsize=1)
def get_config() -> AppConfig:
    config_path = DATA_DIR / "conf.yml"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing configuration file: {config_path}")
    return load_config(config_path)


@lru_cache(maxsize=1)
def get_cache() -> Cache:
    return Cache(CACHE_DIR)


def get_database_path() -> Path:
    return DB_PATH
