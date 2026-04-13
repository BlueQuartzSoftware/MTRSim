#!/usr/bin/env bash
# run_matlab.sh — Run the MATLAB MTR simulation for a given JSON config.
#
# Usage:
#   ./scripts/run_matlab.sh <config.json> [output_dir]
#
# Example:
#   ./scripts/run_matlab.sh configs/test_tiny_default.json output/test_tiny_default
#
# Requirements:
#   - MATLAB R2025b installed at /Applications/MATLAB_R2025b.app
#   - Working directory should be the project root (MTRSim/)
#
# Timing:
#   Writes matlab_timing.json to the output directory with wall-clock and
#   simulation-only elapsed times.

set -euo pipefail

MATLAB_BIN="/Applications/MATLAB_R2025b.app/bin/matlab"
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

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
echo "  MATLAB Simulation"
echo "  Config:  $CONFIG_PATH"
echo "  Output:  $OUTPUT_DIR"
echo "========================================"

# Capture wall-clock time for the entire MATLAB invocation (includes startup)
WALL_START=$(python3 -c "import time; print(time.time())")

"$MATLAB_BIN" -nodisplay -nosplash -batch \
    "addpath('$PROJECT_ROOT/matlab'); addpath('$PROJECT_ROOT/data'); run_simulation_from_json('$CONFIG_PATH', '$OUTPUT_DIR')" \
    2>&1 | tee "$OUTPUT_DIR/matlab_console.log"

WALL_END=$(python3 -c "import time; print(time.time())")
WALL_ELAPSED=$(python3 -c "print(round($WALL_END - $WALL_START, 4))")

# Parse simulation-only time from MATLAB console output
# (grep -oP not available on macOS, so use python)
SIM_TIME=$(python3 -c "
import re
try:
    text = open('$OUTPUT_DIR/matlab_console.log').read()
    m = re.search(r'Total simulation time:\s+([\d.]+)', text)
    print(m.group(1) if m else '0')
except: print('0')
")

# Write timing JSON
python3 -c "
import json, sys
data = {
    'implementation': 'matlab',
    'config': '$CONFIG_BASENAME',
    'wall_clock_seconds': $WALL_ELAPSED,
    'simulation_seconds': float('$SIM_TIME'),
    'notes': 'wall_clock includes MATLAB startup/shutdown overhead; simulation_seconds is algorithm-only'
}
with open('$OUTPUT_DIR/matlab_timing.json', 'w') as f:
    json.dump(data, f, indent=2)
"

echo ""
echo "MATLAB simulation complete. Output in: $OUTPUT_DIR"
echo "  Wall-clock time: ${WALL_ELAPSED}s  (simulation only: ${SIM_TIME}s)"
ls -la "$OUTPUT_DIR"/matlab_*
