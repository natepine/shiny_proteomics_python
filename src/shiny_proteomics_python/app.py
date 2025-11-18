from __future__ import annotations

from starlette.applications import Starlette
from starlette.routing import Mount

from shiny_proteomics_python.builder.app import create_builder_app
from shiny_proteomics_python.viewer.app import create_viewer_app


builder_app = create_builder_app()
viewer_app = create_viewer_app()

app = Starlette(
    routes=[
        Mount("/tmtmosaic/IsoBuilder", builder_app),
        Mount("/tmtmosaic/IsoParser", viewer_app),
        Mount("/", builder_app),
    ]
)
