#!/bin/bash
# =============================================================================
# Sirius 13-Image Demonstration Suite
# Validates spectral emission pipeline, exotic metrics, and cinematic presets
# =============================================================================

SIRIUS="./bin/Sirius.Build/bin/Sirius"

if [ ! -f "$SIRIUS" ]; then
    echo "Error: Sirius binary not found at $SIRIUS"
    echo "Build first: cmake -B bin/Sirius.Build && cmake --build bin/Sirius.Build"
    exit 1
fi

OUT="renders/demo"
mkdir -p "$OUT"

PASS=0
FAIL=0
W=1920
H=1080
SPP=64

render() {
    local name="$1"
    shift
    echo -n "  [$((PASS+FAIL+1))/13] $name... "
    if $SIRIUS render "$@" 2>&1 | tail -1; then
        PASS=$((PASS+1))
    else
        echo "FAILED"
        FAIL=$((FAIL+1))
    fi
}

echo ""
echo "=== Sirius Demonstration Suite (13 images) ==="
echo "    Output: $OUT/"
echo "    Resolution: ${W}x${H}, ${SPP} spp"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Category 1: Spectral Emission (blackbody colour, NT temperature)
# ─────────────────────────────────────────────────────────────────────────────
echo "── Spectral Emission ──"

render "Schwarzschild NT planar disk" \
    -m Schwarzschild -w $W -h $H -s $SPP -d 30 -i 80 --fov 60 \
    --temperature-model NovikovThorne --disk-temperature 50000 \
    -o "$OUT/01_schwarzschild_nt.ppm"

render "Kerr a=0.9 NT planar (Doppler asymmetry)" \
    -m Kerr -a 0.9 -w $W -h $H -s $SPP -d 30 -i 80 --fov 60 \
    --temperature-model NovikovThorne --disk-temperature 50000 \
    -o "$OUT/02_kerr09_nt.ppm"

render "Kerr a=0.998 NT planar (near-extremal)" \
    -m Kerr -a 0.998 -w $W -h $H -s $SPP -d 20 -i 75 --fov 45 \
    --temperature-model NovikovThorne --disk-temperature 50000 \
    -o "$OUT/03_kerr0998_nt.ppm"

render "Schwarzschild Shakura-Sunyaev (comparison)" \
    -m Schwarzschild -w $W -h $H -s $SPP -d 30 -i 80 --fov 60 \
    --temperature-model ShakuraSunyaev --disk-temperature 6500 \
    -o "$OUT/04_schwarzschild_ss.ppm"

# ─────────────────────────────────────────────────────────────────────────────
# Category 2: Viewing Angles
# ─────────────────────────────────────────────────────────────────────────────
echo "── Viewing Angles ──"

render "Kerr a=0.9 edge-on (i=89 deg)" \
    -m Kerr -a 0.9 -w $W -h $H -s $SPP -d 30 -i 89 --fov 60 \
    -o "$OUT/05_kerr09_edgeon.ppm"

render "Kerr a=0.9 face-on (i=15 deg)" \
    -m Kerr -a 0.9 -w $W -h $H -s $SPP -d 30 -i 15 --fov 60 \
    -o "$OUT/06_kerr09_faceon.ppm"

# ─────────────────────────────────────────────────────────────────────────────
# Category 3: Volumetric Disk
# ─────────────────────────────────────────────────────────────────────────────
echo "── Volumetric Disk ──"

render "Kerr a=0.9 volumetric (turbulence + corona)" \
    -m Kerr -a 0.9 -w $W -h $H -s $SPP -d 25 -i 75 --fov 50 \
    --volumetric --h-over-r 0.15 --turbulence --corona \
    -o "$OUT/07_kerr09_volumetric.ppm"

render "Schwarzschild volumetric (thick disk)" \
    -m Schwarzschild -w $W -h $H -s $SPP -d 30 -i 80 --fov 60 \
    --volumetric --h-over-r 0.2 --turbulence \
    -o "$OUT/08_schwarzschild_volumetric.ppm"

# ─────────────────────────────────────────────────────────────────────────────
# Category 4: Cinematic Preset
# ─────────────────────────────────────────────────────────────────────────────
echo "── Cinematic Preset ──"

render "Kerr a=0.9 cinematic (bloom, film)" \
    -m Kerr -a 0.9 -w $W -h $H -s $SPP -d 25 -i 75 --fov 50 \
    --cinematic --film --film-preset Interstellar \
    -o "$OUT/09_kerr09_cinematic.ppm"

render "Kerr a=0.998 cinematic close-up" \
    -m Kerr -a 0.998 -w $W -h $H -s $SPP -d 15 -i 70 --fov 35 \
    --cinematic --film --film-preset Interstellar \
    -o "$OUT/10_kerr0998_cinematic.ppm"

# ─────────────────────────────────────────────────────────────────────────────
# Category 5: Exotic Metrics
# ─────────────────────────────────────────────────────────────────────────────
echo "── Exotic Metrics ──"

render "Morris-Thorne wormhole (b0=1.0)" \
    -m Morris-Thorne --throat-radius 1.0 -w $W -h $H -s $SPP \
    -d 15 -i 80 --fov 60 \
    -o "$OUT/11_wormhole.ppm"

render "Alcubierre warp drive (vs=0.5)" \
    -m Alcubierre --warp-velocity 0.5 --bubble-radius 1.0 \
    -w $W -h $H -s $SPP -d 20 -i 80 --fov 90 \
    -o "$OUT/12_alcubierre.ppm"

# ─────────────────────────────────────────────────────────────────────────────
# Category 6: High-Quality Reference
# ─────────────────────────────────────────────────────────────────────────────
echo "── Reference ──"

render "Kerr a=0.9 high-quality reference (256 spp)" \
    -m Kerr -a 0.9 -w $W -h $H -s 256 -d 25 -i 75 --fov 50 \
    --cinematic --film --film-preset Interstellar \
    -o "$OUT/13_kerr09_reference.ppm"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== Results: $PASS passed, $FAIL failed out of 13 ==="
echo "    Output directory: $OUT/"
ls -lh "$OUT"/*.ppm 2>/dev/null | awk '{print "    " $5 "  " $NF}'
