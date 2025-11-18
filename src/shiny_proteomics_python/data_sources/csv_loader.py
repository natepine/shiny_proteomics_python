from __future__ import annotations

from pathlib import Path
from typing import Literal
import pandas as pd


def load_table(path: Path, kind: Literal["csv", "tsv"]) -> pd.DataFrame:
    if kind == "csv":
        return pd.read_csv(path)
    if kind == "tsv":
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported file type: {kind}")
