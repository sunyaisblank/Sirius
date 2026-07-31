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
#   VK_ICD_FILENAMES=<mesa>/build/src/microsoft/vulkan/dzn_devenv_icd.x86_64.json
# Dozen presents the GPU as Vulkan 1.2 over D3D12 (SPIR-V 1.5 kernels load
# as-is) and, on the 780M, exposes shaderFloat64, so all three rungs run on
# silicon. Performance through the D3D12 translation is indicative, not
# native; treat native-driver numbers (Windows AMD Vulkan) as the benchmark
# of record when both are available.
#
# Pass criteria (a JSON attestation and complete log are written on success):
#   1. Full ctest estate green.
#   2. The governed render completes inside a 2048 MiB budget cap at the
#      specification's exact 5616x4096 IMAX frame, with the tile plan and wall
#      time logged.
#   3. Each precision rung renders (fp64 may decline if the device lacks
#      shaderFloat64 — a loud decline is a pass for that rung, silence is not).
#   4. The CPU/Vulkan parity suites (part of the estate) pass on the real
#      driver, not just Lavapipe.

set -euo pipefail

PRESET="${SIRIUS_PRESET:-linux-gcc}"
EXPECTED_DEVICE="${SIRIUS_EXPECT_DEVICE_PATTERN:-Radeon 780M}"
ATTESTATION_DOMAINS="${SIRIUS_ATTESTATION_DOMAINS:-physical-radeon-780m,physical-imax-5616x4096}"
OUT_ROOT="${SIRIUS_HARDWARE_OUTPUT_DIR:-renders/hardware-validation}"
RUN_ID="$(date +%Y%m%d-%H%M%S)"

if [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
    echo "source tree is dirty; commit or remove every tracked and untracked change before attestation"
    exit 1
fi

mkdir -p "$OUT_ROOT"
OUT="$(mktemp -d "$OUT_ROOT/$RUN_ID-XXXXXX")"
RUN_ID="$(basename "$OUT")"
LOG="$OUT/validation-$RUN_ID.log"
DEVICE_JSON="$OUT/device-$RUN_ID.json"
ATTESTATION="$OUT/attestation-$RUN_ID.json"

say() { echo "== $*" | tee -a "$LOG"; }

say "Sirius physical-GPU validation ($(date -u +%FT%TZ))"
say "host: $(uname -a)"
say "expected physical device pattern: $EXPECTED_DEVICE"
say "attestation domains: $ATTESTATION_DOMAINS"

say "[1/6] configure + build ($PRESET)"
cmake --preset "$PRESET" 2>&1 | tee -a "$LOG"
cmake --build --preset "$PRESET" -j"$(nproc)" 2>&1 | tee -a "$LOG"

SIRIUS="./bin/$PRESET/src/sirius/app/sirius"
[ -x "$SIRIUS" ] || { say "sirius binary missing at $SIRIUS"; exit 1; }

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

READINESS=$("$SIRIUS" --json info readiness)
READINESS="$READINESS" python3 - <<'PY'
import json
import os
import sys

data = json.loads(os.environ["READINESS"])
vulkan = data.get("backends", {}).get("vulkan", {})
if not data.get("ready") or not vulkan.get("ready"):
    sys.exit(f"Vulkan runtime is not ready: {vulkan.get('reason', 'unreported reason')}")
PY

say "[3/6] full test estate on selected driver"
(cd "bin/$PRESET" && ctest -j"$(nproc)" --output-on-failure) 2>&1 | tee -a "$LOG"

say "[4/6] governed render: 2048 MiB budget, 5616x4096 IMAX frame, Kerr a=0.9"
RENDER_START_SECONDS="$SECONDS"
SIRIUS_MEMORY_BUDGET_MB=2048 "$SIRIUS" render --backend vulkan \
    -m Kerr -a 0.9 -w 5616 -h 4096 -s 4 -d 30 -i 80 --fov 60 \
    -o "$OUT/kerr_governed.png" 2>&1 | tee -a "$LOG"
RENDER_WALL_SECONDS="$((SECONDS - RENDER_START_SECONDS))"

say "[5/6] precision ladder on the device"
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

say "[6/6] output integrity and machine-readable attestation"
for image in "$OUT/kerr_governed.png" \
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

# Freeze the evidence transcript before hashing it. Subsequent verifier and
# summary messages continue in LOG without invalidating the admitted artifact.
EVIDENCE_LOG="$OUT/validation-transcript-$RUN_ID.log"
cp "$LOG" "$EVIDENCE_LOG"

python3 - "$DEVICE_JSON" "$ATTESTATION" "$EVIDENCE_LOG" "$SELECTED_INDEX" \
    "$(git rev-parse HEAD)" "$PRESET" "$FP64_SUPPORTED" "$ATTESTATION_DOMAINS" \
    "$RENDER_WALL_SECONDS" <<'PY'
import hashlib
import json
import pathlib
import struct
import sys
from datetime import datetime, timezone

(
    device_path,
    output_path,
    log_path,
    selected_index,
    commit,
    preset,
    fp64,
    domain_text,
    wall_seconds,
) = sys.argv[1:]
root = pathlib.Path(output_path).parent
names = ["kerr_governed.png", "kerr_fp32.png", "kerr_fp32-comp.png"]
if fp64 == "1":
    names.append("kerr_fp64.png")
names.extend([pathlib.Path(device_path).name, pathlib.Path(log_path).name])
inventory = json.loads(pathlib.Path(device_path).read_text())
devices = inventory["backends"]["vulkan"]["devices"]
device = next(item for item in devices if item["index"] == int(selected_index))
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
        "runtime_ready": True,
        "fp64_supported": fp64 == "1",
        "imax_render": {
            "width": width,
            "height": height,
            "memory_budget_mib": 2048,
            "wall_seconds": max(1, int(wall_seconds)),
        },
    },
    "artifacts": outputs,
}
pathlib.Path(output_path).write_text(json.dumps(attestation, indent=2) + "\n")
PY

python3 scripts/verify-attestation.py "$ATTESTATION" 2>&1 | tee -a "$LOG"

say "PASS; log: $LOG"
say "PASS; attestation: $ATTESTATION"
