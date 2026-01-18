#!/usr/bin/env bash

set -euoo -pipefail

rm -rf build
mkdir -p build
cmake -S . -B build
cmake --build build

