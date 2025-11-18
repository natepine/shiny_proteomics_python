#!/usr/bin/env bash
set -euo pipefail
DATA_DIR=${1:-$(pwd)/local-data}
if [[ "$DATA_DIR" != /* ]]; then
  DATA_DIR="$(pwd)/$DATA_DIR"
fi
mkdir -p "$DATA_DIR"
CONF_TEMPLATE="$(pwd)/conf.tmp"
TARGET_CONF="$DATA_DIR/conf.yml"
if [[ ! -f "$TARGET_CONF" ]]; then
  cp "$CONF_TEMPLATE" "$TARGET_CONF"
  CREATED_MSG="Copied conf.tmp to $TARGET_CONF. Please edit this file before starting the app."
else
  CREATED_MSG="$TARGET_CONF already exists. Skipping copy."
fi
cat <<MSG
Local data directory prepared at: $DATA_DIR
$CREATED_MSG

Edit conf.yml with your MassPike servers, API keys, species, and IsoBuilder admin password before running:\n  open "$TARGET_CONF" in your text editor of choice.

After editing conf.yml, start the stack with:\n  docker compose up -d
MSG
