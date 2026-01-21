#!/usr/bin/env bash

set -euoo -pipefail

rm -rf buildTestY
mkdir -p buildTestY
cmake -S . -B buildTestY
cmake --build buildTestY

