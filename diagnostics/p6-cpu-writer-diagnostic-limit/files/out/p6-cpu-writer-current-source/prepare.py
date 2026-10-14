#!/usr/bin/env python3
"""Prepare the five original CPU writer/film controls; execute no tests."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
CASE = "original-five-p6-cpu-writer-controls"
NAMES = [
    "PNGWriterTest.DecodeRoundTripMatchesSRGBEncoding",
    "EXRRoundTripTests.HDRGradientSurvivesWriteAndRead",
    "EXRRoundTripTests.PpmRgbaBoundaryAppliesExactlyOneSrgbEncode",
    "RenderSessionProbe.CpuKerrRenderProducesValidPpmThroughTheOwnedWriter",
    "RenderSessionProbe.FilmAffectsDisplayOutputButNeverLinearExr",
]


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT).decode().strip()


def identity(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + "\n")


assert not (WORK / "preparation.json").exists() and not (WORK / CASE).exists() and not (WORK / "temporary").exists()
assert git("status", "--porcelain", "--untracked-files=normal") == ""
runner = git("rev-parse", "HEAD")
allowed = [".github/workflows/ci.yml", "scripts/diagnose-macos-retained.py",
           "scripts/governance/build_policy.py", "scripts/macos_metal_probe.py"]
assert sorted(git("diff", "--name-only", PRODUCT, "HEAD").splitlines()) == allowed
previous = ROOT / "attestations/scientific-core/4bcd987/xyz-srgb-authority-current-source/files/out/xyz-srgb-authority-current-source"
for name in ("run_owned.py", "owned_guard.py"):
    shutil.copyfile(previous / name, WORK / name)
sources = [ROOT / "tests/render/png_writer_test.cpp", ROOT / "tests/render/exr_writer_test.cpp",
           ROOT / "tests/render/render_session_probe_test.cpp", ROOT / "tests/support/scoped_temporary_directory.h"]
(WORK / "reference").mkdir()
for path in sources:
    assert path.read_bytes() == subprocess.check_output(["git", "show", PRODUCT + ":" + key(path)], cwd=ROOT)
    shutil.copyfile(path, WORK / "reference" / path.name)
review_path = ROOT / "attestations/fresh-session-review/4bcd987/2026-10-14/review.json"
review = json.loads(review_path.read_text())
p6 = next(row for row in review["model_requirement_ledger"] if row["id"] == "P6")
assert set(NAMES) <= set(p6["not_proved_by_reviewed_current_bundles"])
assert not any(set(NAMES) & set(row["passed_names"]) and "/4bcd987/" in row["path"] for row in review["raw_junit"])
accept_path = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation/files/out/gcc-revision-synchronisation/accepted-result.json"
accept = json.loads(accept_path.read_text())
assert accept["passed"] and accept["source_revision"] == PRODUCT
render = ROOT / "bin/linux-gcc/src/sirius/app/sirius_render_tests"
application = render.parent / "sirius"
recoveries = []
for name, path in (("sirius_render_tests", render), ("sirius", application)):
    prior = review["current_live_byte_bindings"]["tested_artifacts"]["records"][name]
    assert identity(path) == {key: prior[key] for key in ("bytes", "sha256")}
    record = next(row for row in accept["current_version_recovery_records"] if row.get("artifact") == name)
    assert identity(path) == {key: record[key] for key in ("bytes", "sha256")}
    recovery = ROOT / record["recovery"]["path"]
    assert identity(recovery) == identity(path)
    recoveries.append(recovery)
resources = render.parent / "resources"
assert resources.is_dir()
resource_files = sorted(path for path in resources.rglob("*") if path.is_file())
assert resource_files and not any(path.is_symlink() for path in resources.rglob("*"))
cache_path = ROOT / "bin/linux-gcc/CMakeCache.txt"
cache = {}
for line in cache_path.read_text().splitlines():
    match = re.match(r"([^#/:][^:]*):[^=]*=(.*)$", line)
    if match:
        cache[match[1]] = match[2]
expected = {"BUILD_TESTS": "ON", "SIRIUS_MANDATORY_TESTS": "ON", "SIRIUS_REQUIRE_VULKAN_RUNTIME": "ON",
            "SIRIUS_ALIGNMENT_MODE": "qualification", "SIRIUS_CONTRACT_MODE": "2", "SIRIUS_WERROR": "ON",
            "SIRIUS_BUILD_VIEWER": "ON", "SIRIUS_SANITIZERS": "none", "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_HOME_DIRECTORY": str(ROOT), "CMAKE_CXX_COMPILER": "g++-14"}
assert all(cache.get(name) == value for name, value in expected.items())
commands = json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
selected_commands = [next(row for row in commands if Path(row["file"]) == source) for source in sources[:3]]
for command in selected_commands:
    assert all(flag in command["command"] for flag in ("-std=c++2c", "-Werror", "-O3"))
    assert "-ffast-math" not in command["command"]
registration = ROOT / "bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake"
selected = []
for name in NAMES:
    matching = [line for line in registration.read_text().splitlines() if line.startswith("add_test(") and name in line]
    assert len(matching) == 1 and str(render) in matching[0] and "--gtest_filter=" + name in matching[0]
    selected += matching
session_header = ROOT / "src/sirius/render/session/render_session.h"
assert "RenderBackend backend = RenderBackend::Cpu;" in session_header.read_text()
inputs = [render, application, *recoveries, *resource_files, cache_path, ROOT / "bin/linux-gcc/build.ninja",
          ROOT / "bin/linux-gcc/compile_commands.json", registration, review_path, accept_path,
          Path(cache["CMAKE_CTEST_COMMAND"]), Path(shutil.which("python3")).resolve(), *sources,
          ROOT / "src/sirius/render/exr_writer.cpp", ROOT / "src/sirius/render/png_writer.cpp",
          ROOT / "src/sirius/render/session/render_session.cpp", session_header,
          ROOT / "tests/render/CMakeLists.txt", ROOT / "tests/operating_model.json"]
for glob in ("CTestTestfile.cmake", "*_tests.cmake", "*_include.cmake"):
    inputs += list((ROOT / "bin/linux-gcc").rglob(glob))
inputs = list(dict.fromkeys(path.resolve() for path in inputs))
local = [WORK / name for name in ("prepare.py", "run_writers.py", "run_owned.py", "owned_guard.py")]
local += [WORK / "reference" / path.name for path in sources]
write("actions.json", {CASE: {
    "revision": runner, "kind": "check", "expected_exit": 0,
    "limits": {"seconds": 300, "rss_kib": 4 * 1024 * 1024, "reserve_kib": 2 * 1024 * 1024},
    "scope": "Five original CPU writer/film cases with unchanged scene sizes, samples, assertions and cleanup. Preserve logs/XML and natural private-TMPDIR residuals. No GPU, forensic pixel identity, retained release-image, full P6, Mandatory, performance or release claim.",
    "argv": [cache["CMAKE_CTEST_COMMAND"], "--test-dir", "bin/linux-gcc", "--output-on-failure", "--no-tests=error",
             "--output-junit", "{folder}/ctest.xml", "-R", "^(" + "|".join(re.escape(name) for name in NAMES) + ")$"],
    "frozen_inputs": [key(path) for path in inputs + local] + [key(WORK / "preparation.json")],
}})
write("preparation.json", {
    "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "product_source_revision": PRODUCT, "runner_source_revision": runner, "diagnostic_only_difference": allowed,
    "executed": False, "scope": "One original five-case CPU writer/film action; no full P6, forensic identity, retained release-image or Mandatory qualification.",
    "finite_acceptance_scope": {
        "png_roundtrip": "Original 8x4 RGB gradient, stbi decode, dimensions, and maximum per-channel error <=1 byte value against the separate existing IEC encoding formula.",
        "hdr_roundtrip": "Original 16x8 gradient (RGB maxima 3.75/14/0.5), Kerr and Sirius linear RGB metadata, TinyEXR decode/dimensions, and max normalized RGB error <1e-6. Writer and decoder use TinyEXR; no independent-format oracle is inferred.",
        "ppm_boundary": "Original 2x1 RGBA input, exact P6 header and six RGB payload bytes [0,188,255,0,137,255]; original clipping and alpha omission assertions.",
        "kerr_ppm_scene": "Original RenderTo/ProbeConfig: 64x64, four SPP, Kerr mass1/spin0.9, tile8, parallel worker pool; all other SessionConfig defaults unchanged. Require Complete, exact P6 64x64/255 header, full RGB payload, EOF and nonconstant bytes.",
        "film_branch_scene": "Four original 16x12/one-SPP serial Minkowski renders with disk/bloom off and Interstellar film toggled: GetDisplayBuffer snapshots selected by PNG output must differ by >1e-4 in total, and snapshots selected by EXR output must compare exactly equal. This test does not decode its EXR files; the separate HDR roundtrip decodes a file.",
    },
    "reference": {
        "kind": "Unchanged original writer expectations and live-session completion/branch comparisons; no new numerical oracle or stronger gate.",
        "original_test_sources": {key(path): {"git_blob": git("rev-parse", PRODUCT + ":" + key(path)), **identity(path)} for path in sources},
        "source_snapshots": {key(WORK / "reference" / path.name): identity(WORK / "reference" / path.name) for path in sources},
        "cleanup_boundary": "The original tests delete their successful image files and scoped directories. Do not intercept unlink or alter cleanup; successful image payload preservation is not claimed. Keep stdout/XML and any natural residuals.",
    },
    "private_tmpdir": {"path": str(WORK / "temporary"), "created_only_at_execution": True,
        "mode": "0700", "owner_uid": os.getuid(),
        "environment_binding": "run_writers.py sets TMPDIR before the unchanged preserved owner copies its child environment", "cleanup_interception": False},
    "source_and_registration_join": {"test_sources": [key(path) for path in sources], "render_producer": identity(render),
        "application_producer": identity(application), "producer_recovery_archives": [key(path) for path in recoveries],
        "runtime_resource_volume": key(resources), "runtime_resource_files": {key(path): identity(path) for path in resource_files},
        "accepted_source_build_evidence": key(accept_path), "fresh_session_review": key(review_path),
        "cache_required": expected, "original_compile_commands": selected_commands,
        "original_registration_file": key(registration), "original_five_registration_lines": selected,
        "expected_names": NAMES, "cpu_backend_default_source": key(session_header)},
    "ledger_scope": "Five current-source names remain unproved in the reviewed ledger. Historical film/PPM passes retain their original different render-producer bounds; they are not relabelled as current execution.",
    "owner_reuse": {name: {"source": key(previous / name), **identity(WORK / name)} for name in ("run_owned.py", "owned_guard.py")},
    "pins": {key(path): identity(path) for path in inputs + local + [WORK / "actions.json"]},
})
for path in local:
    if path.suffix == ".py":
        ast.parse(path.read_text(), filename=str(path))
print(json.dumps({"prepared_without_execution": True, "cases": len(NAMES), "render_producer": identity(render),
                  "runtime_resource_files": len(resource_files), "preparation": identity(WORK / "preparation.json"),
                  "runner": runner, "command": ["python3", "-B", key(WORK / "run_writers.py")]}, indent=2))
