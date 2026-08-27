#!/usr/bin/env bash
# Physical-GPU validation runbook (exceedance E3, docs/SPECIFICATION.md).
#
# Lavapipe is covered by the required CI runtime profile; only silicon can
# prove physical-device performance and driver behaviour. The Radeon 780M
# (2 GiB governed budget) is the named floor. This script refuses software
# devices and, by default, any device whose name does not contain "Radeon
# 780M". Override SIRIUS_EXPECT_DEVICE_PATTERN only to attest another named
# physical target.
#
# WSL2 route (validated 2026-07-28): stock Ubuntu Mesa ships no Dozen ICD, so
# the physical GPU is invisible to Vulkan even though /dev/dxg is present.
# Build Mesa with -Dvulkan-drivers=microsoft-experimental and point the loader
# at the built manifest:
#   VK_DRIVER_FILES=<mesa>/build/src/microsoft/vulkan/dzn_devenv_icd.x86_64.json
# Dozen presents the GPU as Vulkan 1.2 over D3D12 (SPIR-V 1.5 kernels load
# as-is) and, on the 780M, exposes shaderFloat64, so all three rungs run on
# silicon. Performance through the D3D12 translation is indicative, not
# native; treat native-driver numbers (Windows AMD Vulkan) as the benchmark
# of record when both are available.
#
# Pass criteria (a JSON attestation and complete log are written on success):
#   1. Full ctest estate green.
#   2. The governed P3/P5 scene completes inside a 2048 MiB budget cap at both
#      1920x1080 and the specification's exact 5616x4096 IMAX frame, with each
#      tile plan and wall time logged.
#   3. Each precision rung renders (fp64 may decline if the device lacks
#      shaderFloat64 — a loud decline is a pass for that rung, silence is not).
#   4. The CPU/Vulkan parity suites (part of the estate) pass on the real
#      driver, not just Lavapipe.

set -euo pipefail

PRESET="${SIRIUS_PRESET:-linux-ci}"
EXPECTED_DEVICE="${SIRIUS_EXPECT_DEVICE_PATTERN:-Radeon 780M}"
ATTESTATION_DOMAINS="${SIRIUS_ATTESTATION_DOMAINS:-physical-radeon-780m,physical-imax-5616x4096}"
OUT_ROOT="${SIRIUS_HARDWARE_OUTPUT_DIR:-renders/hardware-validation}"
RUN_ID="$(date +%Y%m%d-%H%M%S)"

if [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
    echo "source tree is dirty; commit or remove every tracked and untracked change before attestation"
    exit 1
fi
SOURCE_REVISION="$(git rev-parse HEAD)"

assert_source_identity() {
    if [[ "$(git rev-parse HEAD)" != "$SOURCE_REVISION" ]] ||
       [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
        echo "source revision or cleanliness changed during hardware attestation"
        exit 1
    fi
}

mkdir -p "$OUT_ROOT"
OUT="$(mktemp -d "$OUT_ROOT/$RUN_ID-XXXXXX")"
OUT="$(cd "$OUT" && pwd -P)"
RUN_ID="$(basename "$OUT")"
LOG="$OUT/validation-$RUN_ID.log"
DEVICE_JSON="$OUT/device-$RUN_ID.json"
TEST_REPORT="$OUT/hardware-tests.xml"
CTEST_INVENTORY="$OUT/ctest-inventory.json"
ATTESTATION="$OUT/attestation-$RUN_ID.json"

say() { echo "== $*" | tee -a "$LOG"; }

say "Sirius physical-GPU validation ($(date -u +%FT%TZ))"
say "host: $(uname -a)"
say "source revision: $SOURCE_REVISION"
say "expected physical device pattern: $EXPECTED_DEVICE"
say "attestation domains: $ATTESTATION_DOMAINS"

say "[1/6] configure + strict qualification build ($PRESET)"
cmake --preset "$PRESET" \
    -DSIRIUS_ALIGNMENT_MODE=qualification \
    -DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON 2>&1 | tee -a "$LOG"
cmake --build --preset "$PRESET" -j"$(nproc)" 2>&1 | tee -a "$LOG"

SIRIUS="./bin/$PRESET/src/sirius/app/sirius"
[ -x "$SIRIUS" ] || { say "sirius binary missing at $SIRIUS"; exit 1; }
QUALIFICATION_EXECUTABLE="$OUT/qualification-sirius.bin"
cp "$SIRIUS" "$QUALIFICATION_EXECUTABLE"
ALIGNMENT_RECEIPT_SOURCE="bin/$PRESET/generated/sirius/alignment_receipt.json"
[ -s "$ALIGNMENT_RECEIPT_SOURCE" ] || {
    say "compiled-authority alignment receipt missing at $ALIGNMENT_RECEIPT_SOURCE"
    exit 1
}
cp "$ALIGNMENT_RECEIPT_SOURCE" "$OUT/alignment_receipt.json"
MANDATORY_GATE_SOURCE="bin/$PRESET/generated/sirius/mandatory_gate.json"
[ -s "$MANDATORY_GATE_SOURCE" ] || {
    say "qualification build-gate receipt missing at $MANDATORY_GATE_SOURCE"
    exit 1
}
cp "$MANDATORY_GATE_SOURCE" "$OUT/mandatory_gate.json"
for gate_artifact in mandatory_gate_junit.xml mandatory_gate_ctest.log; do
    source_path="bin/$PRESET/generated/sirius/$gate_artifact"
    [ -s "$source_path" ] || {
        say "qualification gate artifact missing at $source_path"
        exit 1
    }
done
cp "bin/$PRESET/generated/sirius/mandatory_gate_junit.xml" \
    "$OUT/qualification-gate-junit"
cp "bin/$PRESET/generated/sirius/mandatory_gate_ctest.log" \
    "$OUT/qualification-gate-log"
QUALIFICATION_RESOURCE_ROOT="$(dirname "$SIRIUS")/resources"
QUALIFICATION_PRODUCTS=(
    "operating_model=model/operating_model.json"
    "starfield=assets/Starfield.png"
    "trace_fp32comp_spv=kernels/trace_fp32comp.spv"
    "trace_fp64_spv=kernels/trace_fp64.spv"
    "trace_spv=kernels/trace.spv"
    "viewer_rdsd003a_fragment=shaders/RDSD003A.frag"
    "viewer_rdsd003a_vertex=shaders/RDSD003A.vert"
    "viewer_rdsd004a_fragment=shaders/RDSD004A.frag"
    "viewer_rdsd004a_vertex=shaders/RDSD004A.vert"
    "viewer_rdsd005a_fragment=shaders/RDSD005A.frag"
    "viewer_rdsd005a_vertex=shaders/RDSD005A.vert"
)
for product in "${QUALIFICATION_PRODUCTS[@]}"; do
    logical_name="${product%%=*}"
    relative_path="${product#*=}"
    source_path="$QUALIFICATION_RESOURCE_ROOT/$relative_path"
    [ -s "$source_path" ] || {
        say "qualification product missing at $source_path"
        exit 1
    }
    cp "$source_path" "$OUT/qualification-product-$logical_name"
done

say "[2/6] physical-device and readiness precondition"
DEVICE_INFO=$("$SIRIUS" --json info system)
SELECTION=$(
    DEVICE_INFO="$DEVICE_INFO" EXPECTED_DEVICE="$EXPECTED_DEVICE" python3 - <<'PY'
import json
import os
import sys

data = json.loads(os.environ["DEVICE_INFO"])
devices = data.get("backends", {}).get("vulkan", {}).get("devices", [])
physical = [d for d in devices if d.get("kind") in {"integrated", "discrete"}]
expected = os.environ["EXPECTED_DEVICE"].casefold()
matched = [d for d in physical if expected in str(d.get("name", "")).casefold()]
if not physical:
    sys.exit("no physical Vulkan device was reported; software validation is not an attestation")
if not matched:
    names = ", ".join(str(d.get("name", "<unnamed>")) for d in physical)
    sys.exit(f"no physical device matches {os.environ['EXPECTED_DEVICE']!r}; found: {names}")
selected = matched[0]
print(f"{selected['index']}\t{selected['name']}\t{int(bool(selected.get('supports_fp64')))}")
PY
)
IFS=$'\t' read -r SELECTED_INDEX SELECTED_DEVICE FP64_SUPPORTED <<<"$SELECTION"
export SIRIUS_VULKAN_DEVICE="$SELECTED_INDEX"

# Re-report under the selector that every subsequent test and render inherits.
DEVICE_INFO=$("$SIRIUS" --json info system)
printf '%s\n' "$DEVICE_INFO" | tee -a "$LOG"
printf '%s\n' "$DEVICE_INFO" > "$DEVICE_JSON"
say "selected device[$SELECTED_INDEX]: $SELECTED_DEVICE"

set +e
READINESS=$("$SIRIUS" --json info readiness)
READINESS_STATUS=$?
set -e
if [[ "$READINESS_STATUS" -ne 0 && "$READINESS_STATUS" -ne 1 ]]; then
    say "readiness command failed unexpectedly with status $READINESS_STATUS"
    exit 1
fi
READINESS="$READINESS" SELECTED_INDEX="$SELECTED_INDEX" python3 - <<'PY'
import json
import os
import sys

data = json.loads(os.environ["READINESS"])
vulkan = data.get("backends", {}).get("vulkan", {})
if not data.get("evidence_generation_ready") or not vulkan.get("ready"):
    sys.exit(f"Vulkan runtime is not ready: {vulkan.get('reason', 'unreported reason')}")
if vulkan.get("selected_device_index") != int(os.environ["SELECTED_INDEX"]):
    sys.exit("readiness selected a different Vulkan device")
PY
say "readiness: evidence_generation=true vulkan=true selected_device_index=$SELECTED_INDEX"

say "[3/6] full test estate on selected driver"
ctest --test-dir "bin/$PRESET" --show-only=json-v1 >"$CTEST_INVENTORY"
(cd "bin/$PRESET" && ctest -j"$(nproc)" --output-on-failure \
    --output-junit "$TEST_REPORT") 2>&1 | tee -a "$LOG"

say "[4/6] precision ladder on the device"
for rung in fp32 fp32-comp fp64; do
    say "rung: $rung"
    if SIRIUS_PRECISION="$rung" "$SIRIUS" render --backend vulkan \
        -m Kerr -a 0.9 -w 512 -h 512 -s 4 -d 30 -i 80 --fov 60 \
        -o "$OUT/kerr_$rung.png" 2>&1 | tee -a "$LOG"; then
        say "rung $rung: rendered"
    elif [[ "$rung" == "fp64" && "$FP64_SUPPORTED" -eq 0 ]]; then
        say "rung fp64: declined (device inventory reports no shaderFloat64 support)"
    else
        say "rung $rung: failed"
        exit 1
    fi
done

say "[5/6] governed P3/P5 renders: 2048 MiB budget, 1080p and 5616x4096, Kerr a=0.9"
RENDER_1080_START_SECONDS="$SECONDS"
SIRIUS_MEMORY_BUDGET_MB=2048 "$SIRIUS" render --backend vulkan \
    -m Kerr -a 0.9 -w 1920 -h 1080 -s 4 -d 30 -i 80 --fov 60 --no-disk \
    --beams --starfield point --camera-beta 0.1,0.02,-0.01 \
    --lens ThinLens --focal-length 50 --aperture 2.8 --focus-distance 30 \
    -o "$OUT/kerr_governed_1080p.png" 2>&1 | tee -a "$LOG"
RENDER_1080_WALL_SECONDS="$((SECONDS - RENDER_1080_START_SECONDS))"

RENDER_START_SECONDS="$SECONDS"
SIRIUS_MEMORY_BUDGET_MB=2048 "$SIRIUS" render --backend vulkan \
    -m Kerr -a 0.9 -w 5616 -h 4096 -s 4 -d 30 -i 80 --fov 60 --no-disk \
    --beams --starfield point --camera-beta 0.1,0.02,-0.01 \
    --lens ThinLens --focal-length 50 --aperture 2.8 --focus-distance 30 \
    -o "$OUT/kerr_governed.png" 2>&1 | tee -a "$LOG"
RENDER_WALL_SECONDS="$((SECONDS - RENDER_START_SECONDS))"

say "[6/6] output integrity and machine-readable attestation"
for image in "$OUT/kerr_governed_1080p.png" \
             "$OUT/kerr_governed.png" \
             "$OUT/kerr_fp32.png" \
             "$OUT/kerr_fp32-comp.png"; do
    [ -s "$image" ] || { say "missing or empty required output: $image"; exit 1; }
    sha256sum "$image" | tee -a "$LOG"
done
if [[ "$FP64_SUPPORTED" -eq 1 ]]; then
    [ -s "$OUT/kerr_fp64.png" ] ||
        { say "missing or empty fp64 output on an fp64-capable device"; exit 1; }
    sha256sum "$OUT/kerr_fp64.png" | tee -a "$LOG"
fi
assert_source_identity

# Freeze a canonical evidence transcript before hashing it. Interactive
# progress redraws can leave the final completion record after terminal escape
# text on the same physical line; strip only that presentation prefix so the
# verifier sees the renderer's machine-readable completion event verbatim.
# Subsequent verifier and summary messages continue in LOG without invalidating
# the admitted artifact.
EVIDENCE_LOG="$OUT/validation-transcript-$RUN_ID.log"
sed -E 's/^.*(\[Session\] Vulkan render complete:)/\1/' "$LOG" > "$EVIDENCE_LOG"

python3 - "$DEVICE_JSON" "$TEST_REPORT" "$ATTESTATION" "$EVIDENCE_LOG" \
    "$SELECTED_INDEX" "$SOURCE_REVISION" "$PRESET" "$FP64_SUPPORTED" "$ATTESTATION_DOMAINS" \
    "$RENDER_1080_WALL_SECONDS" "$RENDER_WALL_SECONDS" "$SIRIUS" <<'PY'
import hashlib
import json
import pathlib
import struct
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

(
    device_path,
    test_report_path,
    output_path,
    log_path,
    selected_index,
    commit,
    preset,
    fp64,
    domain_text,
    wall_seconds_1080,
    wall_seconds,
    candidate_executable_path,
) = sys.argv[1:]
root = pathlib.Path(output_path).parent
names = [
    "kerr_governed_1080p.png",
    "kerr_governed.png",
    "kerr_fp32.png",
    "kerr_fp32-comp.png",
]
if fp64 == "1":
    names.append("kerr_fp64.png")
names.extend([
    pathlib.Path(device_path).name,
    pathlib.Path(test_report_path).name,
    pathlib.Path(log_path).name,
    "alignment_receipt.json",
    "mandatory_gate.json",
    "ctest-inventory.json",
    "qualification-sirius.bin",
    "qualification-gate-junit",
    "qualification-gate-log",
])
names.extend([
    "qualification-product-operating_model",
    "qualification-product-starfield",
    "qualification-product-trace_fp32comp_spv",
    "qualification-product-trace_fp64_spv",
    "qualification-product-trace_spv",
    "qualification-product-viewer_rdsd003a_fragment",
    "qualification-product-viewer_rdsd003a_vertex",
    "qualification-product-viewer_rdsd004a_fragment",
    "qualification-product-viewer_rdsd004a_vertex",
    "qualification-product-viewer_rdsd005a_fragment",
    "qualification-product-viewer_rdsd005a_vertex",
])
inventory = json.loads(pathlib.Path(device_path).read_text())
devices = inventory["backends"]["vulkan"]["devices"]
device = next(item for item in devices if item["index"] == int(selected_index))
test_cases = list(ET.parse(test_report_path).getroot().iter("testcase"))
test_report = {
    "cases": len(test_cases),
    "failures": sum(case.find("failure") is not None for case in test_cases),
    "errors": sum(case.find("error") is not None for case in test_cases),
    "skipped": sum(case.find("skipped") is not None for case in test_cases),
}
candidate_payload = pathlib.Path(candidate_executable_path).read_bytes()
if (test_report["cases"] == 0 or test_report["failures"] != 0
        or test_report["errors"] != 0 or test_report["skipped"] != 0):
    raise SystemExit(f"hardware JUnit report is not a complete green estate: {test_report}")
outputs = {}
for name in names:
    path = root / name
    outputs[name] = {
        "path": name,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
with (root / "kerr_governed.png").open("rb") as stream:
    signature = stream.read(24)
if signature[:8] != b"\x89PNG\r\n\x1a\n":
    raise SystemExit("governed output is not a PNG")
width, height = struct.unpack(">II", signature[16:24])
if (width, height) != (5616, 4096):
    raise SystemExit(f"governed PNG dimensions are {width}x{height}, expected 5616x4096")
with (root / "kerr_governed_1080p.png").open("rb") as stream:
    signature_1080 = stream.read(24)
if signature_1080[:8] != b"\x89PNG\r\n\x1a\n":
    raise SystemExit("governed 1080p output is not a PNG")
width_1080, height_1080 = struct.unpack(">II", signature_1080[16:24])
if (width_1080, height_1080) != (1920, 1080):
    raise SystemExit(
        f"governed 1080p PNG dimensions are {width_1080}x{height_1080}, expected 1920x1080"
    )

scene = {
    "metric": "Kerr",
    "spin": 0.9,
    "samples_per_pixel": 4,
    "field_of_view": 60.0,
    "disk_enabled": False,
    "ray_bundles": True,
    "point_starfield": True,
    "star_catalogue_minimum": 100000,
    "point_brightness_scale": 100.0,
    "camera_beta": [0.1, 0.02, -0.01],
    "lens": "ThinLens",
    "focal_length": 50.0,
    "aperture": 2.8,
    "focus_distance": 30.0,
}

attestation = {
    "schema_version": 1,
    "status": "pass",
    "completed_utc": datetime.now(timezone.utc).isoformat(),
    "source_revision": commit,
    "domains": [item for item in domain_text.split(",") if item],
    "platform": {
        "os": inventory["platform"],
        "wsl2": inventory["wsl2"],
    },
    "device": device,
    "claims": {
        "preset": preset,
        "log_path": log_path,
        "test_estate_passed": True,
        "test_report": test_report,
        "qualification_executable": {
            "artifact": "qualification-sirius.bin",
            "bytes": len(candidate_payload),
            "sha256": hashlib.sha256(candidate_payload).hexdigest(),
        },
        "runtime_ready": True,
        "fp64_supported": fp64 == "1",
        "p3_render_1080p": {
            "width": width_1080,
            "height": height_1080,
            "memory_budget_mib": 2048,
            "wall_seconds": max(1, int(wall_seconds_1080)),
            "scene": scene,
        },
        "imax_render": {
            "width": width,
            "height": height,
            "memory_budget_mib": 2048,
            "wall_seconds": max(1, int(wall_seconds)),
            "scene": scene,
        },
    },
    "artifacts": outputs,
}
pathlib.Path(output_path).write_text(json.dumps(attestation, indent=2) + "\n")
PY

python3 scripts/verify-attestation.py "$ATTESTATION" 2>&1 | tee -a "$LOG"
assert_source_identity

say "PASS; log: $LOG"
say "PASS; attestation: $ATTESTATION"
