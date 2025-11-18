FROM python:3.11-slim AS tmtmosaic

ENV APP_DIR="/srv/shiny-proteomics"
ENV DATA_DIR="$APP_DIR/data"
ENV PORT=3838

WORKDIR $APP_DIR
ENV PYTHONPATH="$APP_DIR/src:$PYTHONPATH"

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends build-essential libssl-dev libsqlite3-dev && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY src ./src
COPY README.md ./README.md
COPY conf.tmp ./conf.tmp

RUN mkdir -p $DATA_DIR && \
    chown -R root:root $APP_DIR

EXPOSE 3838

CMD ["uvicorn", "shiny_proteomics_python.app:app", "--host", "0.0.0.0", "--port", "3838"]
