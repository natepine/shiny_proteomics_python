# ShinyProteomics

## Build a new image for ShinyProteomics
You'll need docker/docker compose V2 installed on your machine.

Navigate to the shinyproteomics folder and build the image:

```
cd path\to\ShinyProteomics
docker compose build
```

## Create a config file
ShinyProteomics pulls data directly from MassPike instances. ShinyProteomics depends on newer versions of MassPike (2022+) that include API authentication.

Follow the interactive setup using the appropriate variables based on your instance(s) of MassPike.

```
docker compose run --rm -it app genConf
```

Note: Provide URLs of MassPike servers as domains; protocol and path are appended by TMT Mosaic. For example, MassPike web login happens at `https://<domain>/gfy/www/auth/login.php`.

## Starting the webapp
The main container can be started via `docker compose`.

```
docker compose up -d
```

This command starts the ShinyServer. You can now navigate to:

```
localhost:3838/
```
You should see the default shiny page here which will tell you that Shiny is working properly.

If not, make sure that the container started. It should be called something like `tmtmosaic-app-1`.
```
docker container list -a
```

Or check the docker logs.
```
docker logs <container>
```

## Available build options
The `.env` file allows configuration of `docker-compose.yml` files.

```
PORT=3838
```
PORT defines which host port docker serves the app on.

## Using an App
To get started with the apps, substitute the values you supplied in `.env`. Navigate to:
```
http://localhost:${PORT}/tmtmosaic/<IsoBuilder|IsoParser>
```

# Administrative scripts
ShinyProteomics includes several command line scripts to configure and maintain the app. To execute a script run:
```
docker compose run --rm -it app <script>
```

## Avaliable scripts

- `genConf`: Convenience script for generating the `conf.yml` file.
- `annotations`: Accesses various web resources to update `annotations.db` with the latest data from GO, KEGG, etc. Several sub-commands are available, but the most common use is:
```
docker compose run --rm -it app annotations update [organisms:Human/Mouse etc (this needs to match the name in conf.yml)]
```
- `migrateDB`: Checks for old database schemas and updates to the latest schema. (Runs during container startup)
- `rename`: Updates the database, cache directory names, and `conf.yml` to rename Server or Species entries. To add or remove servers, edit `conf.yml` directly or use the `genConf` script.

To manually create or update `conf.yml`, refer to the [template] (conf.tmp)

## Dev environment
Create a symbolic link in the project root directory so that docker-compose.override.yml -> docker-compose.dev.yml.
Then drop into a docker shell with:
```
docker compose up -d
docker compose attach shiny-dev
```
And start a shiny app
```
R
shiny::runApp("IsoBuilder")
```
Edit source files on your machine, not within docker.

