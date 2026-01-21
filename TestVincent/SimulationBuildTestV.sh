#!/usr/bin/env bash

set -euoo -pipefail

rm -rf buildTestV
mkdir -p buildTestV
cmake -S . -B buildTestV
cmake --build buildTestV

