from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import secrets
import sqlite3
import string


JSON_COLUMNS = (
    "GroupNames",
    "GroupColors",
    "GroupIDs",
    "ColumnClasses",
    "ColumnNames",
    "ColumnIDs",
)


@dataclass
class ViewerMetadata:
    qid: int
    username: str
    dataset: str
    notes: str
    num_samples: int
    num_groups: int
    plex_level: int
    num_plexes: int
    server: str
    species: Optional[str]
    initial_id: str
    are_reps: bool
    is_site_quant: bool
    group_names: List[str]
    group_colors: List[str]
    group_ids: List[int]
    column_classes: List[str]
    column_names: List[str]
    column_ids: List[int]
    key: Optional[str] = None


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS metadata (
    QID           INTEGER NOT NULL,
    Key           TEXT PRIMARY KEY,
    Username      TEXT,
    Dataset       TEXT NOT NULL,
    Notes         TEXT NOT NULL,
    NumSamples    INTEGER NOT NULL,
    NumGroups     INTEGER NOT NULL,
    PlexLevel     INTEGER NOT NULL,
    NumPlexes     INTEGER NOT NULL,
    Server        TEXT,
    Species       TEXT,
    InitialID     TEXT NOT NULL,
    AreReps       INTEGER NOT NULL,
    IsSiteQuant   INTEGER NOT NULL,
    Date          TEXT NOT NULL,
    ColumnClasses TEXT NOT NULL,
    GroupNames    TEXT NOT NULL,
    GroupColors   TEXT NOT NULL,
    GroupIDs      TEXT NOT NULL,
    ColumnNames   TEXT NOT NULL,
    ColumnIDs     TEXT NOT NULL
);
"""

UPLOAD_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS uploadID (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    key TEXT NOT NULL
);
"""


def _conn(path: Path, read_only: bool = False) -> sqlite3.Connection:
    if read_only:
        uri = f"file:{path}?mode=ro"
        return sqlite3.connect(uri, uri=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    return sqlite3.connect(path)


def ensure_schema(path: Path) -> None:
    with _conn(path) as conn:
        conn.execute(SCHEMA_SQL)
        conn.execute(UPLOAD_TABLE_SQL)
        conn.commit()


def _jsonify(record: Dict[str, Any]) -> Dict[str, Any]:
    payload = dict(record)
    for column in JSON_COLUMNS:
        value = payload[column]
        payload[column] = json.dumps(value)
    payload["AreReps"] = int(bool(payload["AreReps"]))
    payload["IsSiteQuant"] = int(bool(payload["IsSiteQuant"]))
    return payload


def _restore_row(row: sqlite3.Row) -> ViewerMetadata:
    data = dict(row)
    for column in JSON_COLUMNS:
        data[column] = json.loads(data[column]) if data[column] else []
    return ViewerMetadata(
        qid=data["QID"],
        username=data.get("Username", ""),
        dataset=data["Dataset"],
        notes=data["Notes"],
        num_samples=data["NumSamples"],
        num_groups=data["NumGroups"],
        plex_level=data["PlexLevel"],
        num_plexes=data["NumPlexes"],
        server=data["Server"],
        species=data.get("Species"),
        initial_id=data["InitialID"],
        are_reps=bool(data["AreReps"]),
        is_site_quant=bool(data["IsSiteQuant"]),
        group_names=data["GroupNames"],
        group_colors=data["GroupColors"],
        group_ids=[int(v) for v in data["GroupIDs"]],
        column_classes=data["ColumnClasses"],
        column_names=data["ColumnNames"],
        column_ids=[int(v) for v in data["ColumnIDs"]],
        key=data["Key"],
    )


def _generate_state_key() -> str:
    alphabet = string.ascii_uppercase + string.digits
    return "".join(secrets.choice(alphabet) for _ in range(32))


def save_metadata(path: Path, metadata: ViewerMetadata) -> str:
    ensure_schema(path)
    payload = asdict(metadata)
    key = metadata.key or _generate_state_key()
    payload["key"] = key
    record = {
        "QID": payload["qid"],
        "Key": key,
        "Username": payload["username"],
        "Dataset": payload["dataset"],
        "Notes": payload["notes"],
        "NumSamples": payload["num_samples"],
        "NumGroups": payload["num_groups"],
        "PlexLevel": payload["plex_level"],
        "NumPlexes": payload["num_plexes"],
        "Server": payload["server"],
        "Species": payload["species"],
        "InitialID": payload["initial_id"],
        "AreReps": payload["are_reps"],
        "IsSiteQuant": payload["is_site_quant"],
        "ColumnClasses": payload["column_classes"],
        "GroupNames": payload["group_names"],
        "GroupColors": payload["group_colors"],
        "GroupIDs": payload["group_ids"],
        "ColumnNames": payload["column_names"],
        "ColumnIDs": payload["column_ids"],
    }
    row = _jsonify(record)
    with _conn(path) as conn:
        conn.execute(SCHEMA_SQL)
        conn.execute(
            "DELETE FROM metadata WHERE Key = ?", (key,)
        )
        conn.execute(
            """
            INSERT INTO metadata
            (QID, Key, Username, Dataset, Notes, NumSamples, NumGroups, PlexLevel, NumPlexes,
             Server, Species, InitialID, AreReps, IsSiteQuant, Date,
             ColumnClasses, GroupNames, GroupColors, GroupIDs, ColumnNames, ColumnIDs)
            VALUES
            (:QID, :Key, :Username, :Dataset, :Notes, :NumSamples, :NumGroups, :PlexLevel, :NumPlexes,
             :Server, :Species, :InitialID, :AreReps, :IsSiteQuant, date('now'),
             :ColumnClasses, :GroupNames, :GroupColors, :GroupIDs, :ColumnNames, :ColumnIDs)
            """,
            row,
        )
        conn.commit()
    return key


def fetch_metadata_by_key(path: Path, key: str) -> Optional[ViewerMetadata]:
    if not path.exists():
        return None
    with _conn(path, read_only=True) as conn:
        conn.row_factory = sqlite3.Row
        result = conn.execute(
            "SELECT * FROM metadata WHERE Key = ?", (key,)
        ).fetchone()
    return _restore_row(result) if result else None


def fetch_metadata_for_dataset(path: Path, qid: int, is_site_quant: bool, server: str) -> Optional[ViewerMetadata]:
    if not path.exists():
        return None
    with _conn(path, read_only=True) as conn:
        conn.row_factory = sqlite3.Row
        result = conn.execute(
            "SELECT * FROM metadata WHERE QID = ? AND IsSiteQuant = ? AND Server = ?",
            (qid, int(is_site_quant), server),
        ).fetchone()
    return _restore_row(result) if result else None


def delete_viewer(path: Path, key: str) -> None:
    if not path.exists():
        return
    with _conn(path) as conn:
        conn.execute("DELETE FROM metadata WHERE Key = ?", (key,))
        conn.commit()


def admin_table(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    with _conn(path, read_only=True) as conn:
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT Date, QID, Username, Dataset, Key, Notes, Server, Species, NumSamples, NumGroups, PlexLevel, AreReps, IsSiteQuant, InitialID FROM metadata"
        ).fetchall()
    return [dict(row) for row in rows]


def next_upload_id(path: Path) -> int:
    ensure_schema(path)
    token = secrets.token_hex(8)
    with _conn(path) as conn:
        conn.execute("INSERT INTO uploadID (key) VALUES (?)", (token,))
        rowid = conn.execute("SELECT id FROM uploadID WHERE key = ?", (token,)).fetchone()[0]
        conn.execute("DELETE FROM uploadID WHERE key = ?", (token,))
        conn.commit()
    return int(rowid)
