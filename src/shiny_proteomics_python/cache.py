from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence
import pickle


class Cache:
    """Filesystem-backed cache that mirrors the original R helper."""

    def __init__(self, base_dir: Path, ext: str = ".pkl") -> None:
        self.base_dir = base_dir
        self.ext = ext
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def _resolve(self, key: Sequence[str]) -> Optional[Path]:
        if not key:
            return None
        path = self.base_dir
        *dirs, filename = key
        for directory in dirs:
            candidate = (path / directory).resolve()
            if not str(candidate).startswith(str(self.base_dir.resolve())):
                return None
            candidate.mkdir(exist_ok=True)
            path = candidate
        filename = f"{filename}{self.ext}" if self.ext else filename
        resolved = (path / filename).resolve()
        if not str(resolved).startswith(str(self.base_dir.resolve())):
            return None
        return resolved

    def read(self, key: Sequence[str]):
        path = self._resolve(key)
        if not path or not path.exists():
            return None
        with path.open("rb") as fh:
            return pickle.load(fh)

    def write(self, data, key: Sequence[str], overwrite: bool = False) -> bool:
        if data is None:
            return False
        path = self._resolve(key)
        if path is None:
            return False
        if path.exists() and not overwrite:
            return False
        with path.open("wb") as fh:
            pickle.dump(data, fh)
        return True

    def delete(self, key: Sequence[str]) -> bool:
        path = self._resolve(key)
        if path and path.exists():
            path.unlink()
            return True
        return False


def mosaic_cache_key(identifier: int | str, is_site_quant: bool, payload: str, origin: str) -> Iterable[str]:
    prefix = "sq" if is_site_quant else "pq"
    filename = f"{prefix}_{identifier}"
    return origin, payload, filename
