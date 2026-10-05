#!/usr/bin/env bash
set -e

echo "Building PrusaSlicer..."
cmake --build build -j4

echo "Build finished successfully."

