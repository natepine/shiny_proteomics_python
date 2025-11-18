from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional
import requests

from shiny_proteomics_python.config import ServerDefinition


@dataclass
class MassPikeClient:
    server: ServerDefinition
    username: str
    api_key: str

    def _request(self, method: str, path: str, **kwargs) -> Dict[str, Any]:
        url = f"{self.server.api_root}/{path.lstrip('/')}"
        response = requests.request(
            method,
            url,
            auth=(self.username, self.api_key),
            verify=self.server.verify_ssl,
            timeout=60,
            **kwargs,
        )
        response.raise_for_status()
        data = response.json()
        if isinstance(data, dict) and data.get("status") == "error":
            raise RuntimeError(data.get("error", "Unknown MassPike error"))
        return data

    def get(self, path: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        return self._request("GET", path, params=params)


LOGIN_SUFFIX = "login"


def authenticate(server: ServerDefinition, username: str, password: str) -> Optional[MassPikeClient]:
    url = f"https://{server.host}/gfy/www/modules/api/v1/{LOGIN_SUFFIX}"
    response = requests.post(
        url,
        data={"username": username, "password": password},
        verify=server.verify_ssl,
        timeout=30,
    )
    response.raise_for_status()
    payload = response.json()
    if payload.get("status") != "success":
        return None
    return MassPikeClient(server=server, username=username, api_key=payload["key"])


def fetch_proteins(client: MassPikeClient, qid: int, is_site_quant: bool) -> Dict[str, Any]:
    if is_site_quant:
        module = "combined_site_quant"
        params = {"all_users": 1, "peptide_parsimony": "UR"}
    else:
        module = "protein_quant"
        params = {"all_users": 1}
    return client.get(f"{module}/{qid}", params=params)


def fetch_peptides(client: MassPikeClient, qid: int, is_site_quant: bool) -> Dict[str, Any]:
    if is_site_quant:
        module = "site_quant_peptide"
        params = {"all_users": 1, "peptide_parsimony": "UR"}
    else:
        module = "protein_quant_peptide"
        params = {"all_users": 1}
    return client.get(f"{module}/{qid}", params=params)


def fetch_normalization(client: MassPikeClient, qid: int, is_site_quant: bool) -> Optional[Dict[str, Any]]:
    module = "site_quant_info" if is_site_quant else "protein_quant_info"
    try:
        return client.get(f"{module}/{qid}", params={"all_users": 1})
    except Exception:
        return None
