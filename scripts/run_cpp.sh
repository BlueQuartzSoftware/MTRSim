#!/usr/bin/env bash
# run_cpp.sh — Run the C++ MTR simulation for a given JSON config.
#
# Usage:
#   ./scripts/run_cpp.sh <config.json> [output_dir]
#
# Example:
#   ./scripts/run_cpp.sh configs/test_tiny_default.json output/test_tiny_default
#
# Requirements:
#   - C++ mtrsim binary built at the expected location
#   - Working directory should be the project root (MTRSim/)
#
# Timing:
#   Writes cpp_timing.json to the output directory with wall-clock elapsed time.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="/Users/mjackson/Workspace5/Build/mtrsim-Rel"
MTRSIM_BIN="$BUILD_DIR/app/mtrsim"

if [[ ! -x "$MTRSIM_BIN" ]]; then
    echo "ERROR: mtrsim binary not found at $MTRSIM_BIN"
    echo "Build it first with:"
    echo "  cd $BUILD_DIR && cmake --build . --target all"
    exit 1
fi

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <config.json> [output_dir]"
    exit 1
fi

CONFIG_PATH="$1"
# Make config path absolute if it isn't already
if [[ "$CONFIG_PATH" != /* ]]; then
    CONFIG_PATH="$PROJECT_ROOT/$CONFIG_PATH"
fi

# Default output dir: output/<config_basename_without_extension>
CONFIG_BASENAME="$(basename "$CONFIG_PATH" .json)"
OUTPUT_DIR="${2:-$PROJECT_ROOT/output/${CONFIG_BASENAME}}"

# Make output dir absolute if it isn't already
if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="$PROJECT_ROOT/$OUTPUT_DIR"
fi

mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "  C++ Simulation"
echo "  Config:  $CONFIG_PATH"
echo "  Output:  $OUTPUT_DIR"
echo "========================================"

# Capture wall-clock time
WALL_START=$(python3 -c "import time; print(time.time())")

# Run from project root so relative paths in config (e.g. data/simulation_ODF.h5) resolve
cd "$PROJECT_ROOT"
"$MTRSIM_BIN" --config "$CONFIG_PATH" --output "$OUTPUT_DIR"

WALL_END=$(python3 -c "import time; print(time.time())")
WALL_ELAPSED=$(python3 -c "print(round($WALL_END - $WALL_START, 4))")

# Rename outputs with cpp_ prefix
if [[ -f "$OUTPUT_DIR/sim_results.csv" ]]; then
    mv "$OUTPUT_DIR/sim_results.csv" "$OUTPUT_DIR/cpp_sim_results.csv"
fi
if [[ -f "$OUTPUT_DIR/sim_IPF_map.png" ]]; then
    mv "$OUTPUT_DIR/sim_IPF_map.png" "$OUTPUT_DIR/cpp_sim_IPF_map.png"
fi

# Write timing JSON
python3 -c "
import json
data = {
    'implementation': 'cpp',
    'config': '$CONFIG_BASENAME',
    'wall_clock_seconds': $WALL_ELAPSED,
    'simulation_seconds': $WALL_ELAPSED,
    'notes': 'C++ has negligible startup overhead; wall_clock ~ simulation time'
}
with open('$OUTPUT_DIR/cpp_timing.json', 'w') as f:
    json.dump(data, f, indent=2)
"

echo ""
echo "C++ simulation complete. Output in: $OUTPUT_DIR"
echo "  Wall-clock time: ${WALL_ELAPSED}s"
ls -la "$OUTPUT_DIR"/cpp_*
