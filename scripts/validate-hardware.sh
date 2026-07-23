#!/usr/bin/env bash
# Physical-GPU validation runbook (exceedance E3, docs/SPECIFICATION.md).
#
# Everything in the test estate already runs on Lavapipe software Vulkan; what
# only silicon can prove is performance and driver behaviour on the real
# device — the Radeon 780M (2 GiB) is the named floor target. Run this script
# on the target machine from the repository root; it builds, runs the full
# gate, then exercises the governed-budget and precision-ladder paths on the
# physical device and records everything under renders/hardware-validation/.
#
# Pass criteria (record the log with the machine name):
#   1. Full ctest estate green.
#   2. The governed render completes inside a 2048 MiB budget cap at IMAX-class
#      resolution, with the tile plan and wall time logged.
#   3. Each precision rung renders (fp64 may decline if the device lacks
#      shaderFloat64 — a loud decline is a pass for that rung, silence is not).
#   4. The CPU/Vulkan parity suites (part of the estate) pass on the real
#      driver, not just Lavapipe.

set -euo pipefail

PRESET="${SIRIUS_PRESET:-linux-gcc}"
OUT="renders/hardware-validation"
LOG="$OUT/validation-$(date +%Y%m%d-%H%M%S).log"
mkdir -p "$OUT"

say() { echo "== $*" | tee -a "$LOG"; }

say "Sirius physical-GPU validation ($(date -u +%FT%TZ))"
say "host: $(uname -a)"

say "[1/5] configure + build ($PRESET)"
cmake --preset "$PRESET" 2>&1 | tail -2 | tee -a "$LOG"
cmake --build --preset "$PRESET" -j"$(nproc)" 2>&1 | tail -2 | tee -a "$LOG"

SIRIUS="./bin/$PRESET/src/sirius/app/sirius"
[ -x "$SIRIUS" ] || { say "sirius binary missing at $SIRIUS"; exit 1; }

say "[2/5] full test estate"
(cd "bin/$PRESET" && ctest -j"$(nproc)" --output-on-failure) 2>&1 | tail -4 | tee -a "$LOG"

say "[3/5] device inventory"
DEVICE_INFO=$("$SIRIUS" info system --json 2>&1)
printf '%s\n' "$DEVICE_INFO" | tee -a "$LOG"
FP64_SUPPORTED=0
if grep -q '"supports_fp64": true' <<<"$DEVICE_INFO"; then
    FP64_SUPPORTED=1
fi

say "[4/5] governed render: 2048 MiB budget, IMAX-class frame, Kerr a=0.9"
SIRIUS_MEMORY_BUDGET_MB=2048 "$SIRIUS" render --backend vulkan \
    -m Kerr -a 0.9 -w 4096 -h 2864 -s 4 -d 30 -i 80 --fov 60 \
    -o "$OUT/kerr_governed.png" 2>&1 | tee -a "$LOG"

say "[5/5] precision ladder on the device"
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

say "done; log at $LOG"
