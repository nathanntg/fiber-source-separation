#!/usr/bin/env bash

# set -o pipefail # get errors even when piped
# set -o xtrace # debug
set -o nounset # prevent using undeclared variables
set -o errexit # exit on command fail; allow failure: || true

clang tMCimgLOT.c -O3 -Wall -DANGLE_FROM_NA -DSOURCE_RADIUS -o tMCimgLOT
