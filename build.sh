#!/usr/bin/env bash
set -e

echo "Building PrusaSlicer..."
cmake --build build --target prusa-slicer -j4

echo "Build finished successfully."
