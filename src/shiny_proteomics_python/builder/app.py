from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import List, Optional
import pandas as pd
import plotly.express as px
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui

from shiny_proteomics_python import runtime
from shiny_proteomics_python.cache import mosaic_cache_key
from shiny_proteomics_python.colors import COLOR_PRESETS, preset_colors
from shiny_proteomics_python.config import AppConfig
from shiny_proteomics_python.data_sources import csv_loader
from shiny_proteomics_python.data_sources.masspike import authenticate, fetch_proteins
from shiny_proteomics_python.database import (
    ViewerMetadata,
    admin_table,
    next_upload_id,
    save_metadata,
)
from shiny_proteomics_python.datasets import (
    ColumnMetadata,
    ProteinData,
    compute_plex_structure,
    derive_column_metadata,
    guess_initial_id,
    split_protein_tables,
)


NAME_LIMIT = 30


def _species_choices(config: AppConfig):
    if not config.annotated_species:
        return None
    return list(config.annotated_species.keys())


def _column_table_to_csv(df: pd.DataFrame) -> str:
    buffer = StringIO()
    df.to_csv(buffer, index=False)
    return buffer.getvalue()


def _csv_to_dataframe(payload: str) -> pd.DataFrame:
    return pd.read_csv(StringIO(payload))


def _group_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    groups = sorted(df["Group"].unique())
    group_names = [f"Group {idx}" for idx in groups]
    return pd.DataFrame({"Group": groups, "Group Name": group_names})


def _class_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    classes = sorted(df["Class"].unique())
    return pd.DataFrame({"Class": classes, "Class Name": classes})


def _default_group_colors(n: int) -> List[str]:
    return preset_colors("Basic", n)


def _apply_group_colors(n: int, palette: str, existing: Optional[List[str]]) -> List[str]:
    if existing and len(existing) == n:
        return existing
    return preset_colors(palette, n)


def _sync_tables_to_inputs(session: Session, state: reactive.Values) -> None:
    if state.column_table is not None:
        session.send_input_message("column_table_text", {"value": _column_table_to_csv(state.column_table)})
    if state.group_table is not None:
        session.send_input_message("group_table_text", {"value": _column_table_to_csv(state.group_table)})
    if state.class_table is not None:
        session.send_input_message("class_table_text", {"value": _column_table_to_csv(state.class_table)})


def _sync_color_inputs(session: Session, state: reactive.Values) -> None:
    if state.group_table is None:
        return
    for idx in range(len(state.group_table)):
        color = state.group_colors[idx] if idx < len(state.group_colors) else "#CCCCCC"
        session.send_input_message(f"group_color_{idx + 1}", {"value": color})


def _read_color_inputs(input: Inputs, state: reactive.Values) -> List[str]:
    if state.group_table is None:
        return []
    colors: List[str] = []
    for idx in range(len(state.group_table)):
        field = f"group_color_{idx + 1}"
        value = None
        try:
            value = input[field]()  # type: ignore[index]
        except Exception:
            pass
        if not value:
            value = state.group_colors[idx] if idx < len(state.group_colors) else "#CCCCCC"
        colors.append(value)
    return colors


def _plot_preview(df: pd.DataFrame, column_table: pd.DataFrame) -> px.bar:
    used = column_table[column_table["Use"]]
    if used.empty:
        return px.bar(title="No samples selected")
    columns = used["Source"].tolist()
    preview = df[columns].mean().reset_index()
    preview.columns = ["Source", "Intensity"]
    preview["Sample"] = used["Name"].tolist()
    fig = px.bar(preview, x="Sample", y="Intensity")
    fig.update_layout(margin=dict(l=40, r=10, t=40, b=80))
    return fig


def _validate_tables(column_df: pd.DataFrame, group_df: pd.DataFrame) -> List[str]:
    errors: List[str] = []
    used = column_df[column_df["Use"]]
    if used.empty:
        errors.append("At least one sample must be marked as 'Use'.")
        return errors
    for series in (used["Name"], group_df["Group Name"]):
        duplicates = series[series.duplicated()].unique()
        if len(duplicates) > 0:
            errors.append(
                f"Duplicate names detected: {', '.join(map(str, duplicates))}."
            )
    long_names = used["Name"][used["Name"].str.len() > NAME_LIMIT]
    if len(long_names) > 0:
        errors.append(
            f"Sample names must be ≤ {NAME_LIMIT} characters (violations: {', '.join(long_names)})."
        )
    return errors


def _app_ui(config: AppConfig):
    species_choices = _species_choices(config)
    sidebar = ui.sidebar(
        ui.input_select("workflow", "Section", [
            ("source", "Data Source"),
            ("name", "Display Names"),
            ("color", "Colors"),
            ("notes", "Notes"),
            ("review", "Finalize"),
            ("viewers", "My Viewers"),
        ], selected="source"),
        open=True,
    )
    data_tab = ui.layout_column_wrap(
        ui.card(
            ui.card_header("Authentication"),
            ui.input_select("server", "Server", config.server_choices, selected=next(iter(config.server_choices.values()))),
            ui.input_text("username", "Username"),
            ui.input_password("password", "Password"),
            ui.input_action_button("login", "Login", class_="btn-primary"),
            ui.output_text("login_status"),
        ),
        ui.card(
            ui.card_header("Data Source"),
            ui.input_radio_buttons(
                "data_type",
                "Source",
                {"server": "Server", "upload": "CSV/TSV"},
                selected="server",
                inline=True,
            ),
            ui.input_numeric("qid", "Quant ID", value=None),
            ui.input_checkbox("is_site_quant", "Site Quant", value=False),
            ui.input_checkbox("refresh_cache", "Refresh Cache", value=False),
            ui.input_file("quant_file", "Protein Quant", multiple=False, accept=[".csv", ".tsv"]),
            ui.input_file("peptide_file", "Peptide Quant (optional)", multiple=False, accept=[".csv", ".tsv"]),
            ui.input_radio_buttons("species", "Species", species_choices, inline=True) if species_choices else ui.div(),
            ui.input_action_button("submit_source", "Check Data", class_="btn-success"),
            ui.output_text("source_feedback"),
        ),
    )
    table_editors = ui.layout_column_wrap(
        ui.card(
            ui.card_header("Columns"),
            ui.input_text_area("column_table_text", "Editable CSV", rows=12),
            ui.input_action_button("apply_tables", "Apply Edits"),
        ),
        ui.card(
            ui.card_header("Groups"),
            ui.input_text_area("group_table_text", "Editable CSV", rows=6),
            ui.input_action_button("group_make_unique", "Make Names Unique"),
        ),
        ui.card(
            ui.card_header("Classes"),
            ui.input_text_area("class_table_text", "Editable CSV", rows=6),
        ),
        ui.card(
            ui.card_header("Preview"),
            ui.output_plot("preview_plot"),
        ),
    )
    color_tab = ui.layout_column_wrap(
        ui.card(
            ui.input_select("palette", "Preset", list(COLOR_PRESETS.keys()), selected="Basic"),
            ui.output_ui("group_color_inputs"),
        ),
        ui.card(
            ui.card_header("Protein Preview"),
            ui.output_plot("color_preview"),
        ),
    )
    notes_tab = ui.layout_column_wrap(
        ui.card(
            ui.input_text("dataset_name", "Dataset Name"),
            ui.input_text_area("dataset_notes", "Notes", rows=5),
            ui.input_checkbox("are_reps", "Perfect Replicates", value=False),
        )
    )
    review_tab = ui.layout_column_wrap(
        ui.card(
            ui.card_header("Summary"),
            ui.output_text_verbatim("review_text"),
            ui.input_action_button("create_viewer", "Create Viewer", class_="btn-primary"),
            ui.output_ui("viewer_link"),
        )
    )
    viewers_tab = ui.layout_column_wrap(
        ui.card(
            ui.input_password("admin_password", "Admin Password"),
            ui.input_action_button("refresh_viewers", "Refresh"),
            ui.output_data_frame("viewer_table"),
        )
    )
    return ui.page_sidebar(
        sidebar,
        ui.navset_tab(
            ui.nav_panel("Data Source", data_tab, value="source"),
            ui.nav_panel("Names", table_editors, value="name"),
            ui.nav_panel("Colors", color_tab, value="color"),
            ui.nav_panel("Notes", notes_tab, value="notes"),
            ui.nav_panel("Finalize", review_tab, value="review"),
            ui.nav_panel("My Viewers", viewers_tab, value="viewers"),
            id="workflow_tabs",
        ),
        title="TMT Editor",
    )


def _builder_server(config: AppConfig):
    cache = runtime.get_cache()
    database_path = runtime.get_database_path()

    def server(input: Inputs, output: Outputs, session: Session):
        state = reactive.Values(
            client=None,
            proteins=None,
            column_table=None,
            group_table=None,
            class_table=None,
            group_colors=[],
            metadata={},
            viewer_key=None,
        )

        @output
        @render.text
        def login_status():
            if state.client is None:
                return "Not authenticated"
            return f"Authenticated as {state.client.username}"

        @output
        @render.text
        def source_feedback():
            meta = state.metadata
            if not meta:
                return ""
            return f"Loaded {meta.get('NumSamples', 0)} samples."

        @output
        @render.plot
        def preview_plot():
            req(state.proteins is not None and state.column_table is not None)
            return _plot_preview(state.proteins.quant, state.column_table)

        @output
        @render.plot
        def color_preview():
            req(state.proteins is not None and state.column_table is not None)
            used = state.column_table[state.column_table["Use"]]
            cols = used["Name"].tolist()
            if not cols:
                return _plot_preview(state.proteins.quant, state.column_table)
            data = state.proteins.quant[cols].iloc[:5]
            fig = px.line(data.transpose())
            fig.update_layout(showlegend=False)
            return fig

        @output
        @render.ui
        def group_color_inputs():
            if state.group_table is None:
                return ui.div("Load data to edit colors.")
            controls = []
            for idx, (_, row) in enumerate(state.group_table.iterrows(), start=1):
                label = row["Group Name"]
                color = state.group_colors[idx - 1] if idx - 1 < len(state.group_colors) else "#CCCCCC"
                controls.append(ui.input_text(f"group_color_{idx}", label, value=color))
            return ui.TagList(*controls)

        @reactive.effect
        @reactive.event(input.login)
        def _login():
            server_host = req(input.server())
            username = req(input.username()).strip()
            password = req(input.password())
            server_def = config.servers[server_host]
            client = authenticate(server_def, username, password)
            if client is None:
                session.notification_show("Login failed", type="warning")
            else:
                state.client = client
                session.notification_show("Login successful", type="message")

        def _load_from_server(qid: int, is_site_quant: bool, refresh: bool) -> pd.DataFrame:
            req(state.client is not None)
            server_name = state.client.server.host
            cache_key = mosaic_cache_key(qid, is_site_quant, "raw", server_name)
            if not refresh:
                cached = cache.read(cache_key)
                if isinstance(cached, pd.DataFrame):
                    return cached
            raw = fetch_proteins(state.client, qid, is_site_quant)
            df = pd.DataFrame(raw)
            cache.write(df, cache_key, overwrite=True)
            return df

        def _load_from_upload(kind: str, path: Path) -> pd.DataFrame:
            return csv_loader.load_table(path, kind)

        def _initialize_tables(proteins: ProteinData):
            col_meta = derive_column_metadata(proteins.quant)
            state.column_table = col_meta.table.copy()
            state.group_table = _group_dataframe(state.column_table)
            state.class_table = _class_dataframe(state.column_table)
            state.group_colors = _default_group_colors(len(state.group_table))
            _sync_tables_to_inputs(session, state)
            _sync_color_inputs(session, state)

        @reactive.effect
        @reactive.event(input.submit_source)
        def _load_data():
            data_type = input.data_type()
            is_site = bool(input.is_site_quant())
            refresh = bool(input.refresh_cache())
            if data_type == "server":
                qid = req(input.qid())
                df = _load_from_server(int(qid), is_site, refresh)
                server_name = state.client.server.host if state.client else ""
            else:
                uploaded = req(input.quant_file())
                path = Path(uploaded["datapath"])
                ext = path.suffix.lower().lstrip(".")
                df = _load_from_upload(ext, path)
                server_name = "Upload"
            proteins = split_protein_tables(df, is_site)
            state.proteins = proteins
            _initialize_tables(proteins)
            num_plexes, plex_level = compute_plex_structure(state.column_table)
            metadata = {
                "Server": server_name,
                "IsSiteQuant": is_site,
                "NumSamples": len(state.column_table),
                "NumGroups": len(state.group_table),
                "NumPlexes": num_plexes,
                "PlexLevel": plex_level,
                "InitialID": guess_initial_id(proteins.info),
                "QID": int(input.qid() or next_upload_id(database_path)),
            }
            state.metadata = metadata
            state.viewer_key = None

        @reactive.effect
        @reactive.event(input.apply_tables)
        def _apply_tables_event():
            req(state.column_table is not None)
            try:
                column_df = _csv_to_dataframe(req(input.column_table_text()))
                group_df = _csv_to_dataframe(req(input.group_table_text()))
                class_df = _csv_to_dataframe(req(input.class_table_text()))
            except Exception as exc:
                session.notification_show(f"Failed to parse tables: {exc}", type="error")
                return
            if "Source" not in column_df.columns and state.column_table is not None:
                source_map = dict(zip(state.column_table["Order"], state.column_table["Source"]))
                column_df["Source"] = column_df["Order"].map(source_map)
            errors = _validate_tables(column_df, group_df)
            if errors:
                session.notification_show("\n".join(errors), type="warning")
                return
            state.column_table = column_df
            state.group_table = group_df
            state.class_table = class_df
            state.group_colors = _apply_group_colors(len(group_df), input.palette(), state.group_colors)
            _sync_tables_to_inputs(session, state)
            _sync_color_inputs(session, state)

        @reactive.effect
        @reactive.event(input.palette)
        def _palette_changed():
            if state.group_table is None:
                return
            state.group_colors = preset_colors(input.palette(), len(state.group_table))
            _sync_color_inputs(session, state)

        @reactive.effect
        @reactive.event(input.create_viewer)
        def _create_viewer():
            req(state.column_table is not None and state.group_table is not None)
            if not state.metadata:
                session.notification_show("Load data before creating a viewer", type="error")
                return
            dataset_name = req(input.dataset_name())
            notes = req(input.dataset_notes())
            if not dataset_name or not notes:
                session.notification_show("Dataset name and notes are required", type="warning")
                return
            column_df = state.column_table[state.column_table["Use"]]
            column_ids = [int(v) for v in column_df["Order"].tolist()]
            metadata = state.metadata.copy()
            species_input = getattr(input, "species", None)
            metadata.update({
                "Dataset": dataset_name,
                "Notes": notes,
                "NumGroups": len(state.group_table),
                "AreReps": bool(input.are_reps()),
                "Species": species_input() if callable(species_input) else None,
                "Username": state.client.username if state.client else "",
            })
            colors = _read_color_inputs(input, state)
            if colors:
                state.group_colors = colors
            viewer = ViewerMetadata(
                qid=metadata["QID"],
                username=metadata.get("Username", ""),
                dataset=dataset_name,
                notes=notes,
                num_samples=len(column_df),
                num_groups=metadata["NumGroups"],
                plex_level=metadata["PlexLevel"],
                num_plexes=metadata["NumPlexes"],
                server=metadata["Server"],
                species=metadata.get("Species"),
                initial_id=metadata["InitialID"],
                are_reps=metadata.get("AreReps", False),
                is_site_quant=metadata["IsSiteQuant"],
                group_names=state.group_table["Group Name"].tolist(),
                group_colors=state.group_colors,
                group_ids=[int(v) for v in column_df["Group"].tolist()],
                column_classes=column_df["Class"].tolist(),
                column_names=column_df["Name"].tolist(),
                column_ids=column_ids,
            )
            key = save_metadata(database_path, viewer)
            state.viewer_key = key
            origin = viewer.server or "Upload"
            cache_key = mosaic_cache_key(viewer.qid, viewer.is_site_quant, "proteins", origin)
            if state.proteins is not None:
                source_cols = column_df["Source"].tolist()
                data = state.proteins.quant[source_cols].copy()
                rename_map = dict(zip(source_cols, column_df["Name"].tolist()))
                data = data.rename(columns=rename_map)
                cache.write(data, cache_key, overwrite=True)
            session.notification_show("Viewer created", type="message")

        @output
        @render.ui
        def viewer_link():
            if not state.viewer_key:
                return ui.div()
            href = f"/{config.viewer_subdir}?sessid={state.viewer_key}"
            return ui.div(ui.a("Open Viewer", href=href, target="_blank"))

        def _viewer_rows():
            input.refresh_viewers()
            data = admin_table(database_path)
            if not data:
                return pd.DataFrame()
            df = pd.DataFrame(data)
            if config.super_user_password and input.admin_password() != config.super_user_password:
                username = state.client.username if state.client else None
                if username:
                    df = df[df["Username"] == username]
                else:
                    df = pd.DataFrame()
            return df

        @output
        @render.data_frame
        def viewer_table():
            df = _viewer_rows()
            if df.empty:
                return pd.DataFrame({"Message": ["No viewers found"]})
            df = df.assign(Link=[f"/{config.viewer_subdir}?sessid={key}" for key in df["Key"]])
            return df

    return server


def create_builder_app() -> App:
    config = runtime.get_config()
    return App(_app_ui(config), _builder_server(config))
