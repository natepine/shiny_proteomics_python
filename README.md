# ShinyProteomics (Python Edition)

This repository contains the Python Shiny implementation of the TMT Mosaic editor (IsoBuilder) and viewer (IsoParser). The app can pull data directly from MassPike instances or ingest CSV/TSV exports, guide users through metadata curation, and publish interactive viewers backed by SQLite.

## Prerequisites

* Docker / Docker Compose v2
* Access to at least one MassPike instance with API credentials
* A `conf.yml` file derived from [`conf.tmp`](conf.tmp)

## Configuration

1. Create a working directory for persistent data (database, cache, configuration).
2. Copy `conf.tmp` to `<data>/conf.yml` and edit the MassPike servers, super-user password, and species list as needed.
3. Mount the directory into the container (see `docker-compose.yml`). The application expects the file at `/srv/shiny-proteomics/data/conf.yml`.

## Building the container

```
cd path/to/shiny_proteomics_python
docker compose build
```

## Running the application

### Prerequisites

The container listens on port `3838` by default (configurable via `.env`). Navigate to `http://localhost:3838/tmtmosaic/IsoBuilder` to build viewers or `http://localhost:3838/tmtmosaic/IsoParser?sessid=<viewer-key>` to open a saved viewer.

### Logs

```
docker logs tmtmosaic-app-1
```

### Environment variables

The `.env` file can override runtime defaults:

```
PORT=5000
```

## Development workflow

The dev compose file mounts the source tree into the container. Create a symlink (`docker-compose.override.yml -> docker-compose.dev.yml`) and run:

```
docker compose up -d
docker compose attach shiny-dev
```

Within the container you can run `uvicorn shiny_proteomics_python.app:app --reload --host 0.0.0.0 --port 3838`.

## Directory structure

* `src/shiny_proteomics_python/builder` – Python Shiny IsoBuilder implementation
* `src/shiny_proteomics_python/viewer` – Python Shiny IsoParser implementation
* `src/shiny_proteomics_python/data_sources` – MassPike and CSV ingestion helpers
* `src/shiny_proteomics_python/database.py` – SQLite persistence utilities
* `conf.tmp` – configuration template

## Notes

* Viewer links embed `sessid` query parameters. If a browser request omits the parameter you can paste the key into the viewer UI.
* Cached datasets are stored under `/srv/shiny-proteomics/data/cache` using pickle serialization.
* Authentication uses MassPike REST API keys; ensure MassPike accounts have access to the datasets you plan to publish.
