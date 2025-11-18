from __future__ import annotations

from urllib.parse import parse_qs
from typing import Optional
import pandas as pd
import plotly.express as px
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui

from shiny_proteomics_python import runtime
from shiny_proteomics_python.cache import mosaic_cache_key
from shiny_proteomics_python.config import AppConfig
from shiny_proteomics_python.data_sources.masspike import MassPikeClient, fetch_proteins
from shiny_proteomics_python.database import fetch_metadata_by_key
from shiny_proteomics_python.datasets import split_protein_tables


def _viewer_ui(config: AppConfig):
    return ui.page_fluid(
        ui.h2("TMT Mosaic Viewer"),
        ui.navset_tab(
            ui.nav_panel(
                "Overview",
                ui.input_text("session_key", "Viewer Key", placeholder="Paste sessid if not provided in URL"),
                ui.input_action_button("load_viewer", "Load Viewer"),
                ui.output_text("overview"),
            ),
            ui.nav_panel(
                "Explore",
                ui.input_selectize("protein_select", "Protein", choices=[]),
                ui.output_plot("protein_plot"),
                ui.output_data_frame("protein_table"),
            ),
            ui.nav_panel(
                "Contact",
                ui.output_text("contact"),
            ),
        ),
    )


def _service_client(config: AppConfig, server_host: Optional[str]) -> Optional[MassPikeClient]:
    if not server_host or server_host not in config.servers:
        return None
    server = config.servers[server_host]
    return MassPikeClient(server=server, username=server.user, api_key=server.key)


def _viewer_server(config: AppConfig):
    cache = runtime.get_cache()
    database_path = runtime.get_database_path()

    def server(input: Inputs, output: Outputs, session: Session):
        state = reactive.Values(metadata=None, proteins=None, quant=None)

        def _query_key() -> Optional[str]:
            search = None
            try:
                search = input[".clientdata_url_search"]()  # type: ignore[index]
            except Exception:
                pass
            if search:
                params = parse_qs(search.lstrip("?"))
                key = params.get("sessid", [None])[0]
                if key:
                    return key
            manual = input.session_key()
            return manual or None

        def _load_metadata():
            key = _query_key()
            if not key:
                return None
            metadata = fetch_metadata_by_key(database_path, key)
            if metadata is None:
                session.notification_show("Viewer not found", type="error")
            return metadata

        def _load_dataset(metadata) -> Optional[pd.DataFrame]:
            origin = metadata.server or "Upload"
            cache_key = mosaic_cache_key(metadata.qid, metadata.is_site_quant, "proteins", origin)
            df = cache.read(cache_key)
            if df is None and metadata.server in config.servers:
                client = _service_client(config, metadata.server)
                if client:
                    raw = fetch_proteins(client, metadata.qid, metadata.is_site_quant)
                    df = pd.DataFrame(raw)
                    cache.write(df, cache_key, overwrite=True)
            return df

        def _hydrate():
            metadata = _load_metadata()
            if metadata is None:
                return
            df = _load_dataset(metadata)
            if df is None:
                session.notification_show("Unable to load dataset", type="error")
                return
            proteins = split_protein_tables(df, metadata.is_site_quant)
            state.metadata = metadata
            state.proteins = proteins
            try:
                state.quant = proteins.quant[metadata.column_names]
            except KeyError:
                state.quant = proteins.quant
            choices = [metadata.initial_id]
            if "GeneSymbol" in proteins.info.columns:
                choices = proteins.info["GeneSymbol"].fillna(proteins.info.get("UniprotID", "")).astype(str).tolist()
            session.send_input_message("protein_select", {"choices": choices, "selected": metadata.initial_id})

        @reactive.effect
        @reactive.event(input.load_viewer)
        def _hydrate_on_click():
            _hydrate()

        @reactive.effect
        def _hydrate_from_query():
            key = _query_key()
            if key and state.metadata is None:
                session.send_input_message("session_key", {"value": key})
                _hydrate()

        @output
        @render.text
        def overview():
            if state.metadata is None:
                return "Provide a viewer key to begin."
            meta = state.metadata
            parts = [
                f"Dataset: {meta.dataset}",
                f"Notes: {meta.notes}",
                f"Server: {meta.server}",
                f"Samples: {meta.num_samples}",
            ]
            return " | ".join(parts)

        def _selected_row():
            if state.proteins is None or state.quant is None:
                return None
            protein_id = input.protein_select()
            if not protein_id and state.metadata:
                protein_id = state.metadata.initial_id
            info = state.proteins.info
            if "GeneSymbol" in info.columns:
                match = info[info["GeneSymbol"] == protein_id]
                if not match.empty:
                    return state.quant.loc[match.index[0]]
            return state.quant.iloc[0]

        @output
        @render.plot
        def protein_plot():
            values = _selected_row()
            req(values is not None)
            fig = px.bar(x=values.index, y=values.values)
            fig.update_layout(xaxis_title="Sample", yaxis_title="Intensity")
            return fig

        @output
        @render.data_frame
        def protein_table():
            values = _selected_row()
            req(values is not None)
            return pd.DataFrame({"Sample": values.index, "Intensity": values.values})

        @output
        @render.text
        def contact():
            if config.admin_email:
                return f"Contact administrator: {config.admin_email}"
            return ""

    return server


def create_viewer_app() -> App:
    config = runtime.get_config()
    return App(_viewer_ui(config), _viewer_server(config))
