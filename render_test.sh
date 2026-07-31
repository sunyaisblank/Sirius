#!/bin/bash
# Two-image operator smoke test. A render counts only when the selected binary
# exits successfully and creates a non-empty output.

set -uo pipefail

SIRIUS="${SIRIUS_BINARY:-./bin/linux-gcc/src/sirius/app/sirius}"
if [ ! -x "$SIRIUS" ]; then
    echo "Error: Sirius binary is not executable: $SIRIUS" >&2
    echo "Build first: cmake --preset linux-gcc && cmake --build --preset linux-gcc" >&2
    exit 1
fi

OUT="${SIRIUS_RENDER_TEST_DIR:-renders}"
mkdir -p "$OUT"
PASS=0
FAIL=0
W="${SIRIUS_RENDER_TEST_WIDTH:-512}"
H="${SIRIUS_RENDER_TEST_HEIGHT:-512}"
SPP="${SIRIUS_RENDER_TEST_SAMPLES:-32}"

run_render() {
    local name="$1"
    local output="$2"
    shift 2
    local log="${output}.log"
    rm -f "$output"
    echo "Rendering $name..."
    if "$SIRIUS" render "$@" --output "$output" >"$log" 2>&1 && [ -s "$output" ]; then
        echo "  -> $output created"
        PASS=$((PASS + 1))
    else
        echo "  -> Failed to render $name" >&2
        tail -n 20 "$log" >&2
        FAIL=$((FAIL + 1))
    fi
}

echo "=== Sirius Black Hole Render Test ==="

run_render "Schwarzschild black hole" "$OUT/schwarzschild.png" \
    --metric Schwarzschild --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 80 --fov 60

run_render "Kerr black hole (spin=0.9)" "$OUT/kerr.png" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 80 --fov 60

echo "=== Results: $PASS passed, $FAIL failed out of 2 ==="
if [ "$FAIL" -ne 0 ]; then
    exit 1
fi
