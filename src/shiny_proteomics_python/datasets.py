from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
import pandas as pd


INFO_CANDIDATES = [
    "UniprotID",
    "Protein Id",
    "ProteinId",
    "GeneSymbol",
    "Gene Symbol",
    "Description",
]


@dataclass
class ProteinData:
    info: pd.DataFrame
    quant: pd.DataFrame
    id_label: str


@dataclass
class ColumnMetadata:
    table: pd.DataFrame
    group_names: List[str]


def split_protein_tables(data: pd.DataFrame, is_site_quant: bool) -> ProteinData:
    numeric_cols = [col for col in data.columns if pd.api.types.is_numeric_dtype(data[col])]
    quant = data[numeric_cols].copy()
    info_cols = [col for col in data.columns if col not in numeric_cols]
    info = data[info_cols].copy()
    if "GeneSymbol" not in info.columns and "Gene Symbol" in info.columns:
        info = info.rename(columns={"Gene Symbol": "GeneSymbol"})
    if "Protein Id" in info.columns and "UniprotID" not in info.columns:
        info = info.rename(columns={"Protein Id": "UniprotID"})
    id_label = "Site" if is_site_quant else "Protein"
    if "GeneSymbol" not in info.columns and "UniprotID" in info.columns:
        info["GeneSymbol"] = info["UniprotID"]
    return ProteinData(info=info, quant=quant, id_label=id_label)


def _parse_column_name(column: str) -> Tuple[str, str]:
    parts = column.split("~")
    if len(parts) >= 2:
        return parts[0], "~".join(parts[1:])
    if "-" in column:
        head, tail = column.split("-", 1)
        return head, tail
    return "Class", column


def derive_column_metadata(quant: pd.DataFrame) -> ColumnMetadata:
    rows = []
    for order, column in enumerate(quant.columns, start=1):
        class_name, sample_name = _parse_column_name(column)
        rows.append({
            "Order": order,
            "Source": column,
            "Class": class_name,
            "Name": sample_name,
            "Group": order,
            "Use": True,
        })
    table = pd.DataFrame(rows)
    groups = sorted(table["Group"].unique())
    group_names = [f"Group {idx}" for idx in groups]
    table["Group"] = _condense_groups(table["Group"].tolist())
    return ColumnMetadata(table=table, group_names=group_names)


def _condense_groups(values: List[int]) -> List[int]:
    mapping = {original: i + 1 for i, original in enumerate(sorted(set(values)))}
    return [mapping[v] for v in values]


def guess_initial_id(info: pd.DataFrame) -> str:
    for column in ("MosaicID", "GeneSymbol", "UniprotID"):
        if column in info.columns:
            value = info[column].iloc[0]
            if pd.notna(value):
                return str(value)
    return ""


def compute_plex_structure(table: pd.DataFrame) -> Tuple[int, int]:
    classes = table["Class"].unique()
    num_plexes = len(classes)
    if num_plexes == 0:
        return 0, 0
    plex_level = int(len(table) / num_plexes)
    return num_plexes, max(1, plex_level)


def to_records(df: pd.DataFrame) -> List[dict]:
    return df.to_dict(orient="records")
