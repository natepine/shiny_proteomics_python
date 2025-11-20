# ShinyProteomics (Python Edition)

This repository contains the Python Shiny implementation of the TMT Mosaic editor (IsoBuilder) and viewer (IsoParser). The app can pull data directly from MassPike instances or ingest CSV/TSV exports, guide users through metadata curation, and publish interactive viewers backed by SQLite.

---

## Quickstart (no Docker experience required)

1. **Install Docker Desktop** (Windows/macOS) or Docker Engine (Linux). Once installed you only need the `docker compose` command.
2. **Download this repository**
   ```
   git clone https://github.com/<your-org>/shiny_proteomics_python.git
   cd shiny_proteomics_python
   ```
   > If you downloaded a ZIP, unzip it and `cd` into the extracted folder instead.
3. **Prepare a data folder and starter config**
   ```
   ./scripts/bootstrap_data_dir.sh
   ```
   This creates `local-data/` next to the compose file and copies [`conf.tmp`](conf.tmp) to `local-data/conf.yml`.
4. **Edit `local-data/conf.yml`** in any text editor:
   * list the MassPike servers you can access
   * paste the REST API keys the builder should use
   * set the admin password that protects the IsoBuilder UI
   * update the species list to match your lab’s datasets
5. **Build (first run only) and start the stack**
   ```
   docker compose up -d --build
   ```
   Docker downloads the base image, installs Python dependencies, and starts the Shiny server. Expect the first build to take a few minutes.
6. **Confirm it is running**
   ```
   docker compose ps
   docker logs -f tmtmosaic-app-1   # stop tailing with Ctrl+C
   ```
7. **Exercise the UIs**
   * Builder: <http://localhost:3838/tmtmosaic/IsoBuilder>
     1. Enter the admin password you configured above.
     2. Choose a MassPike server or upload CSV/TSV exports.
     3. Walk through metadata cleanup, color selection, and publishing.
   * Viewer: <http://localhost:3838/tmtmosaic/IsoParser>
     1. Paste the `sessid` from a builder publish step.
     2. Confirm datasets load, plots render, and filters react.
8. **Stop everything when you are done**
   ```
   docker compose down
   ```
   The `local-data/` folder persists your cache, SQLite DB, and `conf.yml` between runs.

### Common troubleshooting steps

| Symptom | Fix |
| --- | --- |
| `docker compose up` complains about missing `conf.yml` | Re-run `./scripts/bootstrap_data_dir.sh` and make sure `local-data/conf.yml` exists and is readable. |
| Browser cannot connect to `localhost:3838` | Make sure Docker Desktop is running and that no VPN/firewall blocks the port. Re-run `docker compose ps` to ensure the container is `running`. |
| Builder rejects the admin password | Open `local-data/conf.yml`, verify the `builder_admin_password` field, and restart the stack. |

---

## Reference

### Prerequisites

* Docker / Docker Compose v2
* Access to at least one MassPike instance with API credentials
* A `conf.yml` file derived from [`conf.tmp`](conf.tmp) – created for you by `./scripts/bootstrap_data_dir.sh`

### Configuration layout

The application reads `/srv/shiny-proteomics/data/conf.yml` inside the container. With the default compose file this path is backed by the host folder `./local-data`. You can point to another directory by editing `docker-compose.yml` or by running `./scripts/bootstrap_data_dir.sh /absolute/path/to/data`.

### Environment variables

The `.env` file can override runtime defaults. Example:

### Prerequisites

### Logs

```
docker logs tmtmosaic-app-1
```

### Environment variables

The `.env` file can override runtime defaults:

```
PORT=5000
```

### Development workflow

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
