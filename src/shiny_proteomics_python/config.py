from __future__ import annotations
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional
import os
import yaml


@dataclass(frozen=True)
class ServerDefinition:
    host: str
    label: str
    user: str
    key: str
    verify_ssl: bool = True

    @property
    def api_root(self) -> str:
        return f"https://{self.host}/gfy/www/modules/api/v1"


@dataclass(frozen=True)
class AppConfig:
    servers: Dict[str, ServerDefinition]
    viewer_subdir: str
    super_user_password: Optional[str]
    admin_email: Optional[str]
    annotated_species: Mapping[str, Optional[int]]

    @property
    def server_choices(self) -> Mapping[str, str]:
        return {srv.label: srv.host for srv in self.servers.values()}


def _coerce_species(raw: Optional[Mapping[str, object]]) -> Mapping[str, Optional[int]]:
    if not raw:
        return {}
    normalized: Dict[str, Optional[int]] = {}
    for label, taxid in raw.items():
        if taxid in (None, ".na", "NA"):
            normalized[label] = None
        else:
            try:
                normalized[label] = int(taxid)
            except (TypeError, ValueError):
                normalized[label] = None
    return normalized


def _coerce_servers(raw: Mapping[str, Mapping[str, object]]) -> Dict[str, ServerDefinition]:
    servers: Dict[str, ServerDefinition] = {}
    for host, cfg in raw.items():
        user = cfg.get("User")
        key = cfg.get("Key")
        if not user or not key:
            raise ValueError(f"Server '{host}' is missing User/Key values")
        label = cfg.get("Label") or host
        verify = bool(cfg.get("VerifySSL", True))
        servers[host] = ServerDefinition(host=host, user=user, key=key, label=label, verify_ssl=verify)
    if not servers:
        raise ValueError("At least one MassPike server must be configured")
    return servers


def _default_viewer_subdir(raw: Optional[str]) -> str:
    if raw:
        return raw.strip("/")
    app_name = os.getenv("APP_NAME") or "shinyproteomics"
    return f"{app_name}/IsoViewer"


def load_config(path: Path) -> AppConfig:
    raw = yaml.safe_load(path.read_text())
    servers = _coerce_servers(raw.get("SERVERS", {}))
    viewer_subdir = _default_viewer_subdir(raw.get("VIEWER_SUBDIR"))
    annotated_species = _coerce_species(raw.get("ANNOTATED_SPECIES"))
    super_user_password = raw.get("SUPER_USER_PASSWORD")
    admin_email = raw.get("ADMIN_EMAIL")
    return AppConfig(
        servers=servers,
        viewer_subdir=viewer_subdir,
        super_user_password=super_user_password,
        admin_email=admin_email,
        annotated_species=annotated_species,
    )
