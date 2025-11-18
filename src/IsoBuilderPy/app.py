from __future__ import annotations

import io
from pathlib import Path
from typing import Iterable

import pandas as pd
import plotly.express as px
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_plotly

APP_DIR = Path(__file__).parent
SAMPLE_DATA = APP_DIR / "data" / "sample_protein_quant.csv"

CLASS_OPTIONS = [
    "Sample",
    "Control",
    "Reference",
]

COLOR_PRESETS = {
    "Basic": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"],
    "Steve": [
        "#F7FFA1",
        "#756027",
        "#45FF42",
        "#FF9100",
        "#1244D4",
        "#FF0000",
        "#1F1F1F",
        "#FFA1FC",
        "#9903A1",
        "#BBFAF5",
        "#0E661C",
    ],
    "Pastel": ["#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6"],
    "Rainbow": ["#FF0000", "#FFA500", "#FFFF00", "#008000", "#0000FF", "#4B0082"],
    "Warm": ["#FEE8C8", "#FDBB84", "#E34A33", "#B30000"],
    "Earth": ["#8C510A", "#D8B365", "#F6E8C3", "#C7EAE5", "#5AB4AC", "#01665E"],
    "Vibrant": ["#00429D", "#73A2C6", "#E73F74", "#FCCD59", "#3B9AB2"],
    "Light": ["#F0F9E8", "#BAE4BC", "#7BCCC4", "#43A2CA", "#0868AC"],
    "Neon": ["#E6FB04", "#00FF00", "#FF00CC", "#9D00FF", "#00FFFF", "#FF0000"],
    "Dark": ["#8A2BE2", "#000000", "#008B8B", "#006400", "#483D8B", "#B8C709"],
    "Jade": ["#94E8B4", "#72BDA3", "#5E8C61", "#4E6151", "#3B322C", "#CCAD99"],
    "Beach ball": ["#25CED1", "#FCEADE", "#FF8A5B", "#EA526F", "#FFFFFF", "#FFC300"],
}


def load_sample_data() -> pd.DataFrame:
    return pd.read_csv(SAMPLE_DATA)


def load_uploaded_file(file_info: dict) -> pd.DataFrame:
    path = Path(file_info["datapath"])
    ext = path.suffix.lower()
    if ext == ".tsv":
        sep = "\t"
    elif ext == ".csv":
        sep = ","
    else:
        raise ValueError("Please upload a CSV or TSV file.")
    return pd.read_csv(path, sep=sep)


def numeric_sample_columns(df: pd.DataFrame) -> list[str]:
    columns: list[str] = []
    for name in df.columns:
        if pd.api.types.is_numeric_dtype(df[name]):
            columns.append(name)
    return columns


def guess_group(sample_name: str) -> str:
    tokens = sample_name.replace("-", "_").split("_")
    if not tokens:
        return "Group"
    if tokens[0].lower() == "sample" and len(tokens) > 1:
        return tokens[1]
    return tokens[0]


def build_metadata(df: pd.DataFrame) -> pd.DataFrame:
    sample_cols = numeric_sample_columns(df)
    records = []
    for idx, column in enumerate(sample_cols, start=1):
        records.append({
            "Order": idx,
            "Name": column,
            "Class": "Sample",
            "Group": guess_group(column),
        })
    if not records:
        raise ValueError("No numeric intensity columns were found in the dataset.")
    return pd.DataFrame.from_records(records)


def build_group_colors(metadata: pd.DataFrame, preset: str = "Basic") -> pd.DataFrame:
    groups = metadata["Group"].dropna().unique().tolist()
    palette = COLOR_PRESETS.get(preset, COLOR_PRESETS["Basic"])
    colors: list[str] = []
    while len(colors) < len(groups):
        colors.extend(palette)
    trimmed = colors[: len(groups)]
    return pd.DataFrame({"Group": groups, "Color": trimmed})


def table_to_csv(df: pd.DataFrame) -> str:
    buffer = io.StringIO()
    df.to_csv(buffer, index=False)
    return buffer.getvalue().strip()


def parse_metadata(text: str) -> pd.DataFrame:
    if not text.strip():
        raise ValueError("Provide sample metadata before saving.")
    df = pd.read_csv(io.StringIO(text))
    required = {"Order", "Name", "Class", "Group"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns in metadata: {', '.join(sorted(missing))}")
    ordered = df.sort_values("Order").reset_index(drop=True)
    ordered["Order"] = ordered["Order"].astype(int)
    return ordered


def parse_colors(text: str) -> pd.DataFrame:
    if not text.strip():
        raise ValueError("Provide at least one group color before saving.")
    df = pd.read_csv(io.StringIO(text))
    required = {"Group", "Color"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns in color table: {', '.join(sorted(missing))}")
    df["Color"] = df["Color"].fillna("#666666")
    return df


INITIAL_DATA = load_sample_data()
INITIAL_METADATA = build_metadata(INITIAL_DATA)
INITIAL_COLORS = build_group_colors(INITIAL_METADATA)


app_ui = ui.page_navbar(
    ui.nav_panel(
        "Data Source",
        ui.layout_column_wrap(
            2,
            ui.input_radio_buttons(
                "data_mode",
                "Data source",
                {"sample": "Use example data", "upload": "Upload CSV/TSV"},
                selected="sample",
            ),
            ui.input_text("qid", "ProteinQuant ID", placeholder="i.e. 16097"),
            ui.input_checkbox("is_site_quant", "Use site quantification", value=False),
            ui.input_file(
                "quant_file",
                "Protein Quant",
                multiple=False,
                accept=[".csv", ".tsv"],
            ),
            ui.input_action_button("load_data", "Check Data", class_="btn-primary"),
        ),
        ui.hr(),
        ui.card(
            ui.card_header("Detected columns"),
            ui.output_data_frame("data_preview"),
        ),
    ),
    ui.nav_panel(
        "Display Names/Order",
        ui.row(
            ui.column(
                6,
                ui.h4("Sample metadata"),
                ui.input_text_area(
                    "metadata_editor",
                    "Copy/paste metadata (Order, Name, Class, Group)",
                    value=table_to_csv(INITIAL_METADATA),
                    rows=10,
                ),
                ui.input_action_button("save_names", "Save Names", class_="btn-success"),
            ),
            ui.column(
                6,
                ui.h4("Preview"),
                ui.output_data_frame("metadata_table"),
            ),
        ),
    ),
    ui.nav_panel(
        "Colors",
        ui.layout_column_wrap(
            2,
            ui.input_select(
                "color_palette",
                "Palette",
                list(COLOR_PRESETS.keys()),
                selected="Basic",
            ),
            ui.input_action_button("apply_palette", "Apply Preset Colors"),
        ),
        ui.row(
            ui.column(
                6,
                ui.input_text_area(
                    "color_editor",
                    "Group colors (Group, Color)",
                    value=table_to_csv(INITIAL_COLORS),
                    rows=8,
                ),
                ui.input_action_button("save_colors", "Save Colors", class_="btn-success"),
            ),
            ui.column(
                6,
                ui.h4("Color table"),
                ui.output_data_frame("color_table"),
                ui.hr(),
                ui.h4("Plot preview"),
                output_widget("color_plot"),
            ),
        ),
    ),
    ui.nav_panel(
        "Notes",
        ui.layout_column_wrap(
            2,
            ui.input_text("dataset_name", "Dataset Name", placeholder="ex. MCT vs WT"),
            ui.input_checkbox("are_reps", "Perfect replicates", value=False),
        ),
        ui.input_text_area(
            "notes",
            "Notes",
            placeholder="Provide detailed context for this dataset...",
            rows=6,
        ),
    ),
    ui.nav_panel(
        "Finalize",
        ui.h4("Review"),
        ui.row(
            ui.column(6, ui.output_ui("review_information_1")),
            ui.column(6, ui.output_ui("review_information_2")),
        ),
        ui.hr(),
        ui.card(
            ui.card_header("Plot preview"),
            output_widget("review_plot"),
        ),
    ),
    title="IsoBuilder (Python Shiny)",
)


def server(input, output, session):
    dataset = reactive.Value(INITIAL_DATA.copy())
    metadata_state = reactive.Value(INITIAL_METADATA.copy())
    color_state = reactive.Value(INITIAL_COLORS.copy())

    def update_preview_tables() -> None:
        ui.update_text_area("metadata_editor", value=table_to_csv(metadata_state()))
        ui.update_text_area("color_editor", value=table_to_csv(color_state()))

    @reactive.effect
    @reactive.event(input.load_data)
    def _load_data() -> None:
        try:
            if input.data_mode() == "upload":
                file_info = input.quant_file()
                if not file_info:
                    raise ValueError("Upload a CSV or TSV file before loading data.")
                data = load_uploaded_file(file_info[0])
            else:
                data = load_sample_data()
        except ValueError as exc:
            ui.notification_show(str(exc), type="error")
            return

        try:
            metadata = build_metadata(data)
        except ValueError as exc:
            ui.notification_show(str(exc), type="error")
            return

        dataset.set(data)
        metadata_state.set(metadata)
        color_state.set(build_group_colors(metadata, input.color_palette()))
        update_preview_tables()
        ui.notification_show(
            f"Loaded dataset with {len(data)} proteins and {metadata.shape[0]} samples.",
            type="message",
        )

    @reactive.effect
    @reactive.event(input.save_names)
    def _save_names() -> None:
        try:
            parsed = parse_metadata(input.metadata_editor() or "")
        except ValueError as exc:
            ui.notification_show(str(exc), type="error")
            return
        metadata_state.set(parsed)
        ui.notification_show("Saved sample metadata.", type="message")

    @reactive.effect
    @reactive.event(input.save_colors)
    def _save_colors() -> None:
        try:
            parsed = parse_colors(input.color_editor() or "")
        except ValueError as exc:
            ui.notification_show(str(exc), type="error")
            return
        color_state.set(parsed)
        ui.notification_show("Saved color assignments.", type="message")

    @reactive.effect
    @reactive.event(input.apply_palette)
    def _apply_palette() -> None:
        color_state.set(build_group_colors(metadata_state(), input.color_palette()))
        ui.update_text_area("color_editor", value=table_to_csv(color_state()))

    @output
    @render.data_frame
    def data_preview():
        df = dataset()
        return render.DataTable(df.head(10), filters=True, height="350px")

    @output
    @render.data_frame
    def metadata_table():
        df = metadata_state()
        return render.DataTable(df, filters=True, height="400px")

    @output
    @render.data_frame
    def color_table():
        df = color_state()
        return render.DataTable(df, filters=True, height="250px")

    def plot_dataframe() -> pd.DataFrame:
        df = dataset()
        metadata = metadata_state()
        sample_names = [name for name in metadata["Name"] if name in df.columns]
        if not sample_names:
            return pd.DataFrame(columns=["Group", "Intensity"])
        melted = df[sample_names].melt(var_name="Sample", value_name="Intensity")
        merged = melted.merge(metadata, left_on="Sample", right_on="Name", how="left")
        grouped = (
            merged.groupby("Group", dropna=True)["Intensity"].mean().reset_index()
        )
        grouped = grouped[grouped["Group"].notna()]
        return grouped

    def color_map() -> dict[str, str]:
        df = color_state()
        mapping = {}
        for _, row in df.iterrows():
            mapping[str(row["Group"])] = str(row["Color"])
        return mapping

    @output
    @render_plotly
    def color_plot():
        grouped = plot_dataframe()
        if grouped.empty:
            return px.bar(title="No samples available")
        colors = color_map()
        fig = px.bar(
            grouped,
            x="Group",
            y="Intensity",
            color="Group",
            color_discrete_map=colors,
        )
        fig.update_layout(showlegend=False, height=320, margin=dict(t=30, l=20, r=20, b=30))
        return fig

    @output
    @render_plotly
    def review_plot():
        return color_plot()

    @output
    @render.ui
    def review_information_1():
        df = dataset()
        metadata = metadata_state()
        return ui.tags.div(
            ui.tags.p(ui.tags.b("Proteins:"), f" {len(df)}"),
            ui.tags.p(ui.tags.b("Samples:"), f" {metadata.shape[0]}"),
            ui.tags.p(ui.tags.b("Classes:"), ", ".join(sorted(metadata["Class"].unique()))),
        )

    @output
    @render.ui
    def review_information_2():
        colors = color_state()
        dataset_name = input.dataset_name() or "Untitled dataset"
        notes = input.notes() or "No additional notes provided."
        return ui.tags.div(
            ui.tags.p(ui.tags.b("Dataset:"), f" {dataset_name}"),
            ui.tags.p(ui.tags.b("Groups:"), ", ".join(colors["Group"].astype(str))),
            ui.tags.p(ui.tags.b("Notes:"), f" {notes}"),
        )


app = App(app_ui, server)
