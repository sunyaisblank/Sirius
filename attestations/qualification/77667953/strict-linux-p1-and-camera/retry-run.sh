#!/usr/bin/env bash
set -euo pipefail
floor=out/strict-p1-camera/floor
root=$(pwd)
mkdir "$floor/cache" "$floor/tmp"
export VK_DRIVER_FILES=/usr/share/vulkan/icd.d/lvp_icd.json
export SIRIUS_VULKAN_DEVICE_INDEX=0
export MESA_SHADER_CACHE_DIR="$root/$floor/cache"
export TMPDIR="$root/$floor/tmp"
/usr/bin/timeout --kill-after=10s 240s "$floor/replay" "$floor/retry-probe.spv" "$floor/sample.bin" >"$floor/retry.json" 2>"$floor/retry.log"
python3 - "$floor" <<'PY'
from pathlib import Path
import json
import sys

floor = Path(sys.argv[1])
stock, retry = [json.loads((floor/f"{variant}.json").read_text()) for variant in ('stock','retry')]
assert stock[1] == 0 and retry[1] == 1 and stock[2] > 1 >= retry[2]
assert stock[3] < stock[17] and stock[4] < stock[17]
assert retry[3] < retry[17] and retry[4] < retry[17]
assert stock[16:28] == retry[16:28]
PY
set +e
/usr/bin/timeout --kill-after=10s 900s "$floor/render/sirius_render_tests" --gtest_filter=VulkanRenderSession.KerrNearExtremalBardeenBoundaryAt1080p --gtest_output="xml:$floor/retry-render.xml" >"$floor/retry-render.log" 2>&1
status=$?
set -e
printf '%s\n' "$status" >"$floor/retry-render.exit"
exit "$status"
