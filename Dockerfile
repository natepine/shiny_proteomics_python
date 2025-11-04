FROM rocker/shiny:4.3.2 AS tmtmosaic

## install packages
# debian
RUN set -x \
    && apt-get update -qq \
    && apt-get -y --no-install-recommends install \
        libsqlite3-dev \
        libglpk-dev \
        libssl-dev \
    && apt-get remove --purge --auto-remove -y \
    && rm -rf /var/lib/apt/lists/*

# R
RUN install2.r --error --skipinstalled \
#Shiny packages
   #shiny \ already installed
    shinydashboard \
    shinyjs \
    shinymanager \
    shinyWidgets \
#Plotting packages
    ggdendro \
    ggplot2 \
    igraph \
    plotly \
    svglite \
    visNetwork \
#Color
    RColorBrewer \
    colourpicker \
#Data formats
    htmlwidgets \
    jsonlite \
    openxlsx \
    rhandsontable \
#Data manipulation
    DT \
    readr \
    RSQLite \
    stringr \
    && rm -rf /tmp/downloaded_packages

## configure Shiny base image
RUN rm -rf /srv/shiny-server/* /bin/shiny-server
COPY shiny-shim /bin/shiny-server
COPY --chown=shiny:shiny shiny.conf /etc/shiny-server/shiny-server.conf
RUN chmod +x /usr/bin/shiny-server

# setup tmtmosaic environment
ENV APP_DIR="/srv/shiny-server/tmtmosaic"
ENV PATH="$PATH:$APP_DIR/scripts"

WORKDIR $APP_DIR

# prepare filesystem and borrow s6 service's run script
RUN mv /etc/services.d/shiny-server scripts \
    && mkdir --mode=700 data \
    && chown shiny:shiny data

# copy webapp src
COPY src/ .

RUN chmod -R +x scripts

CMD ["/init", "entrypoint.sh"]
