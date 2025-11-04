#!/bin/bash

set -e

chown -R shiny:shiny data/

scripts/validateBuild

EXEC='scripts/run'
if [ ! -f $EXEC ]
then
    EXEC='shiny-server'
fi

exec $EXEC

