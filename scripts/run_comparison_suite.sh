#!/usr/bin/env bash
# run_comparison_suite.sh — Run all test configs through both MATLAB and C++,
# then compare outputs statistically and report timing.
#
# Usage:
#   ./scripts/run_comparison_suite.sh [--matlab-only | --cpp-only | --compare-only]
#
# Options:
#   --matlab-only    Run only the MATLAB simulations
#   --cpp-only       Run only the C++ simulations
#   --compare-only   Skip simulations, just run comparison on existing outputs
#   --configs <pat>  Glob pattern for configs (default: configs/test_*.json)
#
# Example:
#   ./scripts/run_comparison_suite.sh
#   ./scripts/run_comparison_suite.sh --configs "configs/test_tiny*.json"
#   ./scripts/run_comparison_suite.sh --compare-only

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT_DIR="$PROJECT_ROOT/scripts"
OUTPUT_ROOT="$PROJECT_ROOT/output"

# Parse arguments
RUN_MATLAB=true
RUN_CPP=true
RUN_COMPARE=true
CONFIG_PATTERN="configs/test_*.json"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --matlab-only)
            RUN_CPP=false
            RUN_COMPARE=false
            shift
            ;;
        --cpp-only)
            RUN_MATLAB=false
            RUN_COMPARE=false
            shift
            ;;
        --compare-only)
            RUN_MATLAB=false
            RUN_CPP=false
            shift
            ;;
        --configs)
            CONFIG_PATTERN="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

cd "$PROJECT_ROOT"

# Collect config files
CONFIG_FILES=( $CONFIG_PATTERN )

if [[ ${#CONFIG_FILES[@]} -eq 0 ]]; then
    echo "No config files found matching: $CONFIG_PATTERN"
    exit 1
fi

echo "============================================================"
echo "  MTRSim Comparison Suite"
echo "  Configs: ${#CONFIG_FILES[@]} files matching $CONFIG_PATTERN"
echo "  Run MATLAB: $RUN_MATLAB"
echo "  Run C++:    $RUN_CPP"
echo "  Compare:    $RUN_COMPARE"
echo "============================================================"
echo ""

# Track pass/fail
TOTAL=0
MATLAB_FAIL=0
CPP_FAIL=0

for CONFIG in "${CONFIG_FILES[@]}"; do
    CONFIG_NAME="$(basename "$CONFIG" .json)"
    OUT_DIR="$OUTPUT_ROOT/$CONFIG_NAME"
    TOTAL=$((TOTAL + 1))

    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "  [$TOTAL/${#CONFIG_FILES[@]}]  $CONFIG_NAME"
    echo "────────────────────────────────────────────────────────────"

    # Run MATLAB
    if $RUN_MATLAB; then
        echo ""
        echo ">>> Running MATLAB..."
        if "$SCRIPT_DIR/run_matlab.sh" "$CONFIG" "$OUT_DIR" 2>&1; then
            echo "  MATLAB: OK"
        else
            echo "  MATLAB: FAILED"
            MATLAB_FAIL=$((MATLAB_FAIL + 1))
        fi
    fi

    # Run C++
    if $RUN_CPP; then
        echo ""
        echo ">>> Running C++..."
        if "$SCRIPT_DIR/run_cpp.sh" "$CONFIG" "$OUT_DIR" 2>&1; then
            echo "  C++: OK"
        else
            echo "  C++: FAILED"
            CPP_FAIL=$((CPP_FAIL + 1))
        fi
    fi
done

# Run comparison
if $RUN_COMPARE; then
    echo ""
    echo "============================================================"
    echo "  Running Statistical Comparisons"
    echo "============================================================"

    COMPARE_SCRIPT="$SCRIPT_DIR/compare_results.py"

    if [[ ! -f "$COMPARE_SCRIPT" ]]; then
        echo "ERROR: Comparison script not found: $COMPARE_SCRIPT"
        exit 1
    fi

    COMPARE_FAIL=0
    for CONFIG in "${CONFIG_FILES[@]}"; do
        CONFIG_NAME="$(basename "$CONFIG" .json)"
        OUT_DIR="$OUTPUT_ROOT/$CONFIG_NAME"

        MATLAB_CSV="$OUT_DIR/matlab_sim_results.csv"
        CPP_CSV="$OUT_DIR/cpp_sim_results.csv"

        if [[ ! -f "$MATLAB_CSV" ]] || [[ ! -f "$CPP_CSV" ]]; then
            echo "  [$CONFIG_NAME] SKIP — missing output files"
            continue
        fi

        echo ""
        echo "  [$CONFIG_NAME] Comparing..."
        if python3 "$COMPARE_SCRIPT" "$MATLAB_CSV" "$CPP_CSV" "$CONFIG" "$OUT_DIR"; then
            echo "  [$CONFIG_NAME] PASS"
        else
            echo "  [$CONFIG_NAME] FAIL"
            COMPARE_FAIL=$((COMPARE_FAIL + 1))
        fi
    done

    echo ""
    echo "============================================================"
    echo "  Comparison Summary"
    echo "  Total configs:    $TOTAL"
    echo "  MATLAB failures:  $MATLAB_FAIL"
    echo "  C++ failures:     $CPP_FAIL"
    echo "  Compare failures: $COMPARE_FAIL"
    echo "============================================================"
fi

# ── Timing summary ──────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  Timing Summary"
echo "============================================================"

# Collect timing data from all output dirs and print a formatted table
python3 -c "
import json, os, sys

output_root = '$OUTPUT_ROOT'
configs = []

for name in sorted(os.listdir(output_root)):
    d = os.path.join(output_root, name)
    if not os.path.isdir(d):
        continue

    matlab_file = os.path.join(d, 'matlab_timing.json')
    cpp_file = os.path.join(d, 'cpp_timing.json')

    m_wall = m_sim = c_wall = None
    if os.path.isfile(matlab_file):
        with open(matlab_file) as f:
            m = json.load(f)
            m_wall = m.get('wall_clock_seconds', 0)
            m_sim = m.get('simulation_seconds', 0)
    if os.path.isfile(cpp_file):
        with open(cpp_file) as f:
            c = json.load(f)
            c_wall = c.get('wall_clock_seconds', 0)

    if m_wall is not None or c_wall is not None:
        configs.append((name, m_wall, m_sim, c_wall))

if not configs:
    print('  No timing data found.')
    sys.exit(0)

# Header
print(f\"\"\"  {'Config':<30s} {'MATLAB wall':>12s} {'MATLAB sim':>12s} {'C++ wall':>10s} {'Speedup':>10s}\"\"\")
print(f\"\"\"  {'-'*30} {'-'*12} {'-'*12} {'-'*10} {'-'*10}\"\"\")

for name, m_wall, m_sim, c_wall in configs:
    m_wall_s = f'{m_wall:.3f}s' if m_wall is not None else '-'
    m_sim_s  = f'{m_sim:.3f}s'  if m_sim  is not None else '-'
    c_wall_s = f'{c_wall:.3f}s' if c_wall is not None else '-'

    # Speedup: MATLAB sim time / C++ time (algorithm-to-algorithm comparison)
    if m_sim is not None and c_wall is not None and c_wall > 0:
        speedup = m_sim / c_wall
        speedup_s = f'{speedup:.1f}x'
    else:
        speedup_s = '-'

    print(f'  {name:<30s} {m_wall_s:>12s} {m_sim_s:>12s} {c_wall_s:>10s} {speedup_s:>10s}')

print()
print('  MATLAB wall  = total wall-clock including MATLAB startup/shutdown')
print('  MATLAB sim   = algorithm-only time (PGRF + ODF sampling)')
print('  C++ wall     = total wall-clock (negligible startup overhead)')
print('  Speedup      = MATLAB sim / C++ wall')
" 2>&1 || echo "  (timing summary requires python3)"

# Write consolidated timing CSV
python3 -c "
import json, os, csv

output_root = '$OUTPUT_ROOT'
csv_path = os.path.join(output_root, 'timing_summary.csv')

rows = []
for name in sorted(os.listdir(output_root)):
    d = os.path.join(output_root, name)
    if not os.path.isdir(d):
        continue

    matlab_file = os.path.join(d, 'matlab_timing.json')
    cpp_file = os.path.join(d, 'cpp_timing.json')

    m_wall = m_sim = c_wall = ''
    if os.path.isfile(matlab_file):
        with open(matlab_file) as f:
            m = json.load(f)
            m_wall = m.get('wall_clock_seconds', '')
            m_sim = m.get('simulation_seconds', '')
    if os.path.isfile(cpp_file):
        with open(cpp_file) as f:
            c = json.load(f)
            c_wall = c.get('wall_clock_seconds', '')

    if m_wall or c_wall:
        speedup = ''
        if m_sim and c_wall and float(c_wall) > 0:
            speedup = round(float(m_sim) / float(c_wall), 2)
        rows.append([name, m_wall, m_sim, c_wall, speedup])

if rows:
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['config', 'matlab_wall_s', 'matlab_sim_s', 'cpp_wall_s', 'speedup_x'])
        w.writerows(rows)
    print(f'  Timing CSV written: {csv_path}')
" 2>&1 || true

echo ""
echo "Done. All outputs are in: $OUTPUT_ROOT/"
