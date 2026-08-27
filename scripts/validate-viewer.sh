#!/usr/bin/env bash
# Produce revision-bound evidence that a real GLFW window received keyboard and
# pointer events through the host window system. Headless handler tests are a
# prerequisite, but they are intentionally not accepted as this attestation.

set -euo pipefail

PRESET="${SIRIUS_PRESET:-linux-ci}"
EXPECTED_DEVICE="${SIRIUS_EXPECT_DEVICE_PATTERN:-Radeon 780M}"
OUT_ROOT="${SIRIUS_VIEWER_OUTPUT_DIR:-renders/viewer-validation}"
RUN_ID="$(date +%Y%m%d-%H%M%S)"

if [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
    echo "source tree is dirty; commit or remove every tracked and untracked change before attestation"
    exit 1
fi
SOURCE_REVISION="$(git rev-parse HEAD)"

assert_source_identity() {
    if [[ "$(git rev-parse HEAD)" != "$SOURCE_REVISION" ]] ||
       [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
        echo "source revision or cleanliness changed during viewer attestation"
        exit 1
    fi
}

for command in cmake python3 timeout xdotool; do
    command -v "$command" >/dev/null || {
        echo "required viewer-attestation command is unavailable: $command"
        exit 1
    }
done

mkdir -p "$OUT_ROOT"
OUT="$(mktemp -d "$OUT_ROOT/$RUN_ID-XXXXXX")"
OUT="$(cd "$OUT" && pwd -P)"
RUN_ID="$(basename "$OUT")"
LOG="$OUT/validation-$RUN_ID.log"
TRANSCRIPT="$OUT/input-transcript-$RUN_ID.log"
EVIDENCE_LOG="$OUT/validation-transcript-$RUN_ID.log"
DEVICE_JSON="$OUT/device-$RUN_ID.json"
TEST_REPORT="$OUT/viewer-tests.xml"
CTEST_INVENTORY="$OUT/ctest-inventory.json"
ATTESTATION="$OUT/attestation-$RUN_ID.json"

say() { echo "== $*" | tee -a "$LOG"; }

VIEWER_PID=""
stop_viewer() {
    if [[ -n "$VIEWER_PID" ]] && kill -0 "$VIEWER_PID" 2>/dev/null; then
        kill "$VIEWER_PID" 2>/dev/null || true
        wait "$VIEWER_PID" 2>/dev/null || true
    fi
}
trap stop_viewer EXIT

say "Sirius native viewer-input validation ($(date -u +%FT%TZ))"
say "host: $(uname -a)"
say "source revision: $SOURCE_REVISION"

say "[1/6] configure + strict qualification build ($PRESET)"
cmake --preset "$PRESET" \
    -DSIRIUS_ALIGNMENT_MODE=qualification \
    -DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON 2>&1 | tee -a "$LOG"
cmake --build --preset "$PRESET" -j"$(nproc)" 2>&1 | tee -a "$LOG"

SIRIUS="./bin/$PRESET/src/sirius/app/sirius"
[[ -x "$SIRIUS" ]] || { say "sirius binary missing at $SIRIUS"; exit 1; }
QUALIFICATION_EXECUTABLE="$OUT/qualification-sirius.bin"
cp "$SIRIUS" "$QUALIFICATION_EXECUTABLE"
ALIGNMENT_RECEIPT_SOURCE="bin/$PRESET/generated/sirius/alignment_receipt.json"
[[ -s "$ALIGNMENT_RECEIPT_SOURCE" ]] || {
    say "compiled-authority alignment receipt missing at $ALIGNMENT_RECEIPT_SOURCE"
    exit 1
}
ALIGNMENT_RECEIPT="$OUT/alignment_receipt.json"
cp "$ALIGNMENT_RECEIPT_SOURCE" "$ALIGNMENT_RECEIPT"
MANDATORY_GATE_SOURCE="bin/$PRESET/generated/sirius/mandatory_gate.json"
[[ -s "$MANDATORY_GATE_SOURCE" ]] || {
    say "qualification build-gate receipt missing at $MANDATORY_GATE_SOURCE"
    exit 1
}
MANDATORY_GATE="$OUT/mandatory_gate.json"
cp "$MANDATORY_GATE_SOURCE" "$MANDATORY_GATE"
for gate_artifact in mandatory_gate_junit.xml mandatory_gate_ctest.log; do
    source_path="bin/$PRESET/generated/sirius/$gate_artifact"
    [[ -s "$source_path" ]] || {
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
    [[ -s "$source_path" ]] || {
        say "qualification product missing at $source_path"
        exit 1
    }
    cp "$source_path" "$OUT/qualification-product-$logical_name"
done

say "[2/6] physical-device and readiness precondition"
SYSTEM_JSON=$("$SIRIUS" --json info system)
SELECTION=$(
    SYSTEM_JSON="$SYSTEM_JSON" EXPECTED_DEVICE="$EXPECTED_DEVICE" python3 - <<'PY'
import json
import os

data = json.loads(os.environ["SYSTEM_JSON"])
devices = data.get("backends", {}).get("vulkan", {}).get("devices", [])
physical = [item for item in devices if item.get("kind") in {"integrated", "discrete"}]
expected = os.environ["EXPECTED_DEVICE"].casefold()
matched = [
    item for item in physical
    if expected in str(item.get("name", "")).casefold()
]
if not matched:
    names = ", ".join(str(item.get("name", "<unnamed>")) for item in physical)
    raise SystemExit(
        f"no physical viewer device matches {os.environ['EXPECTED_DEVICE']!r}; found: {names}"
    )
selected = matched[0]
print(f"{selected['index']}\t{selected['name']}")
PY
)
IFS=$'\t' read -r SELECTED_INDEX SELECTED_DEVICE <<<"$SELECTION"
export SIRIUS_VULKAN_DEVICE="$SELECTED_INDEX"
SYSTEM_JSON=$("$SIRIUS" --json info system)
printf '%s\n' "$SYSTEM_JSON" | tee -a "$LOG"
printf '%s\n' "$SYSTEM_JSON" >"$DEVICE_JSON"
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

data = json.loads(os.environ["READINESS"])
vulkan = data.get("backends", {}).get("vulkan", {})
if not data.get("evidence_generation_ready") or not vulkan.get("ready"):
    raise SystemExit(f"Vulkan viewer backend is not ready: {vulkan.get('reason', '')}")
if vulkan.get("selected_device_index") != int(os.environ["SELECTED_INDEX"]):
    raise SystemExit("viewer readiness selected a different Vulkan device")
PY
say "readiness: evidence_generation=true vulkan=true selected_device_index=$SELECTED_INDEX"

say "[3/6] full source-available test estate on selected driver"
ctest --test-dir "bin/$PRESET" --show-only=json-v1 >"$CTEST_INVENTORY"
(cd "bin/$PRESET" && ctest -j"$(nproc)" --output-on-failure \
    --output-junit "$TEST_REPORT") 2>&1 | tee -a "$LOG"

say "[4/6] create a native Vulkan window, publish a frame, and inject host input events"
# GLFW otherwise prefers native Wayland under WSLg. This evidence producer
# deliberately selects XWayland so xdotool can discover the window and deliver
# events through the X server's native event queue.
WAYLAND_DISPLAY= SIRIUS_VIEWER_INPUT_TRANSCRIPT="$TRANSCRIPT" \
    "$SIRIUS" view --backend vulkan --width 1024 --height 768 >>"$LOG" 2>&1 &
VIEWER_PID="$!"

WINDOW_IDS="$(timeout 30s xdotool search --sync --onlyvisible --pid "$VIEWER_PID" \
    --name '^Sirius -')" || {
    say "no visible Sirius viewer window appeared within 30 seconds"
    exit 1
}
read -r WINDOW_ID _ <<<"$WINDOW_IDS"
[[ "$WINDOW_ID" =~ ^[0-9]+$ ]] || { say "xdotool returned no numeric window id"; exit 1; }

FRAME_READY=0
for ((attempt = 0; attempt < 600; ++attempt)); do
    if grep -q '^frame-published backend=Vulkan ' "$TRANSCRIPT" 2>/dev/null; then
        FRAME_READY=1
        break
    fi
    kill -0 "$VIEWER_PID" 2>/dev/null || break
    sleep 0.1
done
[[ "$FRAME_READY" -eq 1 ]] || {
    say "the native Vulkan viewer published no frame within 60 seconds"
    exit 1
}

xdotool windowactivate --sync "$WINDOW_ID"
xdotool keydown --window "$WINDOW_ID" w
xdotool keyup --window "$WINDOW_ID" w
xdotool mousemove --window "$WINDOW_ID" 120 120
xdotool mousedown --window "$WINDOW_ID" 1
xdotool mousemove --window "$WINDOW_ID" 180 160
xdotool mouseup --window "$WINDOW_ID" 1
xdotool click --window "$WINDOW_ID" 4
xdotool key --window "$WINDOW_ID" Escape

if ! wait "$VIEWER_PID"; then
    VIEWER_PID=""
    say "viewer exited unsuccessfully after native input delivery"
    exit 1
fi
VIEWER_PID=""

say "[5/6] callback and progressive-frame transcript integrity"
[[ -s "$TRANSCRIPT" ]] || { say "viewer input transcript is missing or empty"; exit 1; }
grep -q '^window-created ' "$TRANSCRIPT" || { say "window creation was not recorded"; exit 1; }
grep -q '^frame-published backend=Vulkan ' "$TRANSCRIPT" || {
    say "no progressive Vulkan frame was recorded"
    exit 1
}
grep -q '^keyboard-callback ' "$TRANSCRIPT" || { say "no keyboard callback was recorded"; exit 1; }
grep -q '^pointer-callback ' "$TRANSCRIPT" || { say "no pointer callback was recorded"; exit 1; }
tee -a "$LOG" <"$TRANSCRIPT"
assert_source_identity

say "[6/6] machine-readable attestation"
# Freeze the evidence log before hashing it. Verifier output is appended only
# to the mutable operator log and therefore cannot invalidate the admitted file.
cp "$LOG" "$EVIDENCE_LOG"
python3 - "$ATTESTATION" "$TRANSCRIPT" "$EVIDENCE_LOG" "$DEVICE_JSON" \
    "$TEST_REPORT" "$ALIGNMENT_RECEIPT" "$MANDATORY_GATE" "$CTEST_INVENTORY" \
    "$SOURCE_REVISION" "$SYSTEM_JSON" "$SELECTED_INDEX" "$PRESET" "$SIRIUS" <<'PY'
import hashlib
import json
import pathlib
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

(
    output_path,
    transcript_path,
    log_path,
    device_path,
    test_report_path,
    alignment_receipt_path,
    mandatory_gate_path,
    ctest_inventory_path,
    revision,
    system_text,
    selected_index,
    preset,
    candidate_executable_path,
) = sys.argv[1:]
output = pathlib.Path(output_path)
root = output.parent
system = json.loads(system_text)
devices = system["backends"]["vulkan"]["devices"]
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
    raise SystemExit(f"viewer JUnit report is not a complete green estate: {test_report}")
transcript_lines = pathlib.Path(transcript_path).read_text(encoding="utf-8").splitlines()
frame_count = sum(line.startswith("frame-published backend=Vulkan ")
                  for line in transcript_lines)
if frame_count < 1:
    raise SystemExit("viewer transcript contains no published Vulkan frame")
artifacts = {}
for source_text in (
    transcript_path,
    log_path,
    device_path,
    test_report_path,
    alignment_receipt_path,
    mandatory_gate_path,
    ctest_inventory_path,
    str(root / "qualification-sirius.bin"),
    str(root / "qualification-gate-junit"),
    str(root / "qualification-gate-log"),
    str(root / "qualification-product-operating_model"),
    str(root / "qualification-product-starfield"),
    str(root / "qualification-product-trace_fp32comp_spv"),
    str(root / "qualification-product-trace_fp64_spv"),
    str(root / "qualification-product-trace_spv"),
    str(root / "qualification-product-viewer_rdsd003a_fragment"),
    str(root / "qualification-product-viewer_rdsd003a_vertex"),
    str(root / "qualification-product-viewer_rdsd004a_fragment"),
    str(root / "qualification-product-viewer_rdsd004a_vertex"),
    str(root / "qualification-product-viewer_rdsd005a_fragment"),
    str(root / "qualification-product-viewer_rdsd005a_vertex"),
):
    source = pathlib.Path(source_text)
    payload = source.read_bytes()
    artifacts[source.name] = {
        "path": source.name,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }

document = {
    "schema_version": 1,
    "status": "pass",
    "completed_utc": datetime.now(timezone.utc).isoformat(),
    "source_revision": revision,
    "domains": ["viewer-native-window-input"],
    "platform": {"os": system["platform"], "wsl2": system["wsl2"]},
    "device": device,
    "claims": {
        "preset": preset,
        "test_estate_passed": True,
        "test_report": test_report,
        "qualification_executable": {
            "artifact": "qualification-sirius.bin",
            "bytes": len(candidate_payload),
            "sha256": hashlib.sha256(candidate_payload).hexdigest(),
        },
        "runtime_ready": True,
        "viewer": {
            "backend": "Vulkan",
            "window_created": True,
            "frame_observed": True,
            "published_frame_count": frame_count,
            "keyboard_event_observed": True,
            "pointer_event_observed": True,
        },
    },
    "artifacts": artifacts,
}
output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
PY

python3 scripts/verify-attestation.py "$ATTESTATION" 2>&1 | tee -a "$LOG"
assert_source_identity
say "PASS; log: $LOG"
say "PASS; attestation: $ATTESTATION"
