"""Inspect terminal native scene evidence and reuse identical EXR decoding."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
revision = sys.argv[1]
assert len(revision) == 40
current = ROOT / "attestations/production-acceptance" / revision[:7]
name = "radeon-moving-thinlens-timing"
identity = json.loads((current / (name + "-identity.json")).read_text())
assert identity["source_revision"] == revision and identity["exit_code"] == 0
assert identity["skipped_testcases"] == 0
assert identity["estate"]["tests"] == "1"
assert all(identity["estate"][key] == "0" for key in ("failures", "errors", "disabled"))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

for suffix, key in ((".log", "log_sha256"), (".xml", "xml_sha256")):
    assert digest(current / (name + suffix)) == identity[key]
records = []
for line in (current / (name + ".log")).read_text().splitlines():
    start = line.find('{"schema":')
    if start >= 0:
        records.append(json.loads(line[start:]))
scene_by_backend = {}
for scene in (record for record in records if record["schema"] == "sirius-render-scene-v1"):
    backend = scene["backend"]
    if backend in scene_by_backend:
        assert scene_by_backend[backend] == scene, "contradictory session/source scene records"
    scene_by_backend[backend] = scene
scenes = list(scene_by_backend.values())
renders = [record for record in records if record["schema"] == "sirius-vulkan-render-v1"]
assert len(scenes) == 2 and len(renders) == 1
render = renders[0]
assert render["route"] == "retained" and render["device_name"] == "AMD Radeon 780M Graphics"
assert (render["width"], render["height"]) == (4, 2)
assert render["allocated_bytes"] <= render["usable_bytes"] <= render["budget_bytes"] == 2147483648
assert all(render[key] == 0 for key in ("target_overshoots", "batch_subdivisions", "safety_fallbacks"))
timing = render["retained_timing"]
assert sum(timing["batch_row_counts"]) == timing["batches"]
assert sum(i * count for i, count in enumerate(timing["batch_row_counts"])) == timing["interval_rows"] + timing["camera_rows"]
assert timing["acceleration_calls"] == render["accepted_intervals"]
baseline_path = ROOT / "attestations/production-acceptance/09a263e/radeon-moving-thinlens-timing-output.json"
baseline = json.loads(baseline_path.read_text())
assert scenes == baseline["scene_evidence"], "scene configuration changed"
outputs = identity["outputs"]
assert len(outputs) == 2
for path, description in outputs.items():
    assert (ROOT / path).stat().st_size == description["bytes"] and digest(ROOT / path) == description["sha256"]
record = {
    "source_revision": revision,
    "scope": "one completed unchanged native4x2/CPU scene; host wall observations and independently checked output; no full-frame/interactive/release qualification",
    "scene_evidence": scenes, "render_evidence": renders, "output_identity": outputs,
    "estate": identity["estate"], "skipped_testcases": identity["skipped_testcases"],
}
prior_decode_path = ROOT / "attestations/production-acceptance/bba118b/radeon-moving-thinlens-output.json"
prior_decode = json.loads(prior_decode_path.read_text())
gpu_path = next(path for path in outputs if Path(path).name == "sirius_retained_moving_kerr.exr")
cpu_path = next(path for path in outputs if Path(path).name == "sirius_cpu_moving_kerr.exr")
if outputs[gpu_path]["sha256"] == prior_decode["gpu"]["sha256"] and outputs[cpu_path]["sha256"] == prior_decode["cpu"]["sha256"]:
    record["reused_independent_decode"] = {
        "path": str(prior_decode_path.relative_to(ROOT)), "sha256": digest(prior_decode_path),
        "decoder_sha256": prior_decode["decoder_sha256"],
        "all_finite": prior_decode["gpu"]["all_finite"] and prior_decode["cpu"]["all_finite"],
        "width": prior_decode["gpu"]["width"], "height": prior_decode["gpu"]["height"],
        "positive_rgb_channels": prior_decode["gpu"]["positive_rgb_channels"],
        "cpu_rgb_signal_l1": prior_decode["cpu"]["rgb_signal_l1"],
        "relative_linear_rgb_error": prior_decode["relative_linear_rgb_error"],
    }
else:
    decoder_path = prior_decode_path.parent / "decode_scanline_exr.py"
    spec = importlib.util.spec_from_file_location("independent_exr", decoder_path)
    decoder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(decoder)
    gpu, cpu = decoder.decode(ROOT / gpu_path), decoder.decode(ROOT / cpu_path)
    assert (gpu["width"], gpu["height"], cpu["width"], cpu["height"]) == (4, 2, 4, 2)
    error = sum(abs(a - b) for channel in "RGB" for a, b in zip(gpu["values"][channel], cpu["values"][channel]))
    relative = error / cpu["rgb_signal_l1"]
    assert relative < .02
    for decoded in (gpu, cpu):
        del decoded["values"]
    record["independent_decode"] = {"decoder_sha256": digest(decoder_path), "gpu": gpu, "cpu": cpu, "relative_linear_rgb_error": relative}
before = baseline["render_evidence"][0]
record["comparison"] = {
    "scope": "two individual renderer steady-clock observations; not a repeated performance distribution or full-frame extrapolation",
    "baseline_path": str(baseline_path.relative_to(ROOT)), "baseline_sha256": digest(baseline_path),
    "baseline_source_revision": baseline["source_revision"],
    "wall_seconds_before": before["wall_seconds"], "wall_seconds_after": render["wall_seconds"],
    "wall_fraction_after_before": render["wall_seconds"] / before["wall_seconds"],
    "accepted_intervals_before": before["accepted_intervals"], "accepted_intervals_after": render["accepted_intervals"],
    "dispatches_before": before["dispatches"], "dispatches_after": render["dispatches"],
    "gather_ms_before": before["retained_timing"]["coalescing_wait_ms"], "gather_ms_after": timing["coalescing_wait_ms"],
    "execute_ms_before": before["retained_timing"]["execute_ms"], "execute_ms_after": timing["execute_ms"],
    "submit_wait_ms_before": before["retained_timing"]["submit_wait_ms"], "submit_wait_ms_after": timing["submit_wait_ms"],
}
assert not (current / (name + "-output.json")).exists()
(current / (name + "-output.json")).write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps({"comparison": record["comparison"], "decode": record.get("reused_independent_decode", record.get("independent_decode"))}, indent=2))
