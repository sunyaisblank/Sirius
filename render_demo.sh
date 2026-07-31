#!/bin/bash
# =============================================================================
# Sirius 13-Image Demonstration Suite
# Validates spectral emission pipeline, exotic metrics, and cinematic presets.
# A case passes only when the renderer succeeds and its declared output exists.
# =============================================================================

set -uo pipefail

SIRIUS="${SIRIUS_BINARY:-./bin/linux-gcc/src/sirius/app/sirius}"
if [ ! -x "$SIRIUS" ]; then
    echo "Error: Sirius binary is not executable: $SIRIUS" >&2
    echo "Build first: cmake --preset linux-gcc && cmake --build --preset linux-gcc" >&2
    exit 1
fi

OUT="${SIRIUS_DEMO_DIR:-renders/demo}"
mkdir -p "$OUT"

PASS=0
FAIL=0
W="${SIRIUS_DEMO_WIDTH:-1920}"
H="${SIRIUS_DEMO_HEIGHT:-1080}"
SPP="${SIRIUS_DEMO_SAMPLES:-64}"
REFERENCE_SPP="${SIRIUS_DEMO_REFERENCE_SAMPLES:-256}"

render() {
    local name="$1"
    local output="$2"
    shift 2
    local log="${output}.log"
    rm -f "$output"
    echo -n "  [$((PASS + FAIL + 1))/13] $name... "
    if "$SIRIUS" render "$@" --output "$output" >"$log" 2>&1 && [ -s "$output" ]; then
        echo "passed"
        PASS=$((PASS + 1))
    else
        echo "FAILED" >&2
        tail -n 20 "$log" >&2
        FAIL=$((FAIL + 1))
    fi
}

echo ""
echo "=== Sirius Demonstration Suite (13 images) ==="
echo "    Output: $OUT/"
echo "    Resolution: ${W}x${H}, ${SPP} spp"
echo ""

echo "── Spectral Emission ──"

render "Schwarzschild NT planar disk" "$OUT/01_schwarzschild_nt.ppm" \
    --metric Schwarzschild --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 80 --fov 60 \
    --temperature-model NovikovThorne --disk-temperature 50000

render "Kerr a=0.9 NT planar (Doppler asymmetry)" "$OUT/02_kerr09_nt.ppm" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 80 --fov 60 \
    --temperature-model NovikovThorne --disk-temperature 50000

render "Kerr a=0.998 NT planar (near-extremal)" "$OUT/03_kerr0998_nt.ppm" \
    --metric Kerr --spin 0.998 --width "$W" --height "$H" --samples "$SPP" \
    --distance 20 --inclination 75 --fov 45 \
    --temperature-model NovikovThorne --disk-temperature 50000

render "Schwarzschild Shakura-Sunyaev (comparison)" "$OUT/04_schwarzschild_ss.ppm" \
    --metric Schwarzschild --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 80 --fov 60 \
    --temperature-model ShakuraSunyaev --disk-temperature 6500

echo "── Viewing Angles ──"

render "Kerr a=0.9 edge-on (i=89 deg)" "$OUT/05_kerr09_edgeon.ppm" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 89 --fov 60

render "Kerr a=0.9 face-on (i=15 deg)" "$OUT/06_kerr09_faceon.ppm" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 15 --fov 60

echo "── Volumetric Disk ──"

render "Kerr a=0.9 volumetric (turbulence + corona)" "$OUT/07_kerr09_volumetric.ppm" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$SPP" \
    --distance 25 --inclination 75 --fov 50 \
    --volumetric --h-over-r 0.15 --turbulence --corona

render "Schwarzschild volumetric (thick disk)" "$OUT/08_schwarzschild_volumetric.ppm" \
    --metric Schwarzschild --width "$W" --height "$H" --samples "$SPP" \
    --distance 30 --inclination 80 --fov 60 \
    --volumetric --h-over-r 0.2 --turbulence

echo "── Cinematic Preset ──"

render "Kerr a=0.9 cinematic (bloom, film)" "$OUT/09_kerr09_cinematic.ppm" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$SPP" \
    --distance 25 --inclination 75 --fov 50 \
    --cinematic --film --film-preset Interstellar

render "Kerr a=0.998 cinematic close-up" "$OUT/10_kerr0998_cinematic.ppm" \
    --metric Kerr --spin 0.998 --width "$W" --height "$H" --samples "$SPP" \
    --distance 15 --inclination 70 --fov 35 \
    --cinematic --film --film-preset Interstellar

echo "── Exotic Metrics ──"

render "Morris-Thorne wormhole (b0=1.0)" "$OUT/11_wormhole.ppm" \
    --metric Morris-Thorne --throat-radius 1.0 --width "$W" --height "$H" \
    --samples "$SPP" --distance 15 --inclination 80 --fov 60 --no-disk

render "Alcubierre warp drive (vs=0.5)" "$OUT/12_alcubierre.ppm" \
    --metric Alcubierre --warp-velocity 0.5 --bubble-radius 1.0 \
    --width "$W" --height "$H" --samples "$SPP" --distance 20 \
    --inclination 80 --fov 90 --no-disk

echo "── Reference ──"

render "Kerr a=0.9 high-quality reference (256 spp)" "$OUT/13_kerr09_reference.ppm" \
    --metric Kerr --spin 0.9 --width "$W" --height "$H" --samples "$REFERENCE_SPP" \
    --distance 25 --inclination 75 --fov 50 \
    --cinematic --film --film-preset Interstellar

echo ""
echo "=== Results: $PASS passed, $FAIL failed out of 13 ==="
echo "    Output directory: $OUT/"
if [ "$FAIL" -ne 0 ]; then
    exit 1
fi
