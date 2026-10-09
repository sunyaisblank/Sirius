"""Run only after coordinator approval; one bounded dispatch per offline module."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

env = os.environ.copy()
env["VK_DRIVER_FILES"] = "/usr/share/vulkan/icd.d/lvp_icd.json"
env["VK_ICD_FILENAMES"] = env["VK_DRIVER_FILES"]
env["SIRIUS_VULKAN_DEVICE"] = "0"
env["MESA_SHADER_CACHE_DIR"] = str(ROOT / "mesa-shader-cache")
processes = subprocess.check_output(["ps", "-eo", "pid,ppid,etime,args"], text=True)
live = [line for line in processes.splitlines()
        if any(needle in line for needle in ("/sirius render", "/sirius_backend_tests",
                                            "/sirius_render_tests", "./gpu_probe "))]
if live:
    raise SystemExit("another possible GPU diagnostic is live: " + repr(live))
report = {"scope": "actual llvmpipe integer primitive dispatch only; not retained/product qualification",
          "pre_dispatch_gpu_process_matches": live, "modules": {},
          "artifact_report_sha256": digest(ROOT / "artifact-report.json"),
          "gpu_probe_sha256": digest(ROOT / "gpu_probe"),
          "gpu_probe_source_sha256": digest(ROOT / "gpu_probe.cpp"),
          "portable_source_sha256": digest(ROOT / "portable_f32.h"),
          "input_sha256": digest(ROOT / "shader-input.bin"),
          "expected_sha256": digest(ROOT / "shader-expected.bin"),
          "icd_manifest_sha256": digest(env["VK_DRIVER_FILES"]),
          "driver_sha256": digest("/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so"),
          "backend_archive_sha256": digest(REPO / "bin/linux-gcc/src/sirius/backend/libsirius_backend.a")}
for stem in ("raw_word_probe", "raw_word_probe_o3"):
    start = time.monotonic()
    args = [str(ROOT / "gpu_probe"), str(ROOT / (stem + ".spv")),
            str(ROOT / "shader-input.bin"), str(ROOT / "shader-expected.bin"), "--dispatch"]
    completed = subprocess.run(args, env=env, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=30, cwd=ROOT)
    (ROOT / (stem + "-gpu.log")).write_text(completed.stdout)
    report["modules"][stem] = {"terminal_exit": completed.returncode,
                              "elapsed_seconds": time.monotonic() - start,
                              "module_sha256": digest(ROOT / (stem + ".spv")),
                              "log_sha256": digest(ROOT / (stem + "-gpu.log")),
                              "output": completed.stdout}
    (ROOT / "gpu-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(stem + ": " + completed.stdout, flush=True)
    if completed.returncode:
        raise SystemExit(completed.returncode)
print("both bounded diagnostics terminal; no further dispatch scheduled", flush=True)
