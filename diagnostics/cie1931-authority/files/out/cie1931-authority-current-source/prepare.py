#!/usr/bin/env python3
"""Prepare the fixed CPU CIE witness action; never execute its tests."""
from pathlib import Path
import ast
import csv
import datetime
from decimal import Decimal
import hashlib
import io
import json
import re
import shutil
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
CASE = "original-three-cie1931-controls"
NAMES = [
    "Cie1931ObserverAuthority.FitIsBoundedAgainstOfficialCieDatasetAcrossRepresentedBand",
    "Cie1931ObserverAuthority.HostSpectralFacadesEvaluateTheirActualWavelengths",
    "Cie1931ObserverAuthority.UnrepresentedWavelengthsFailClosed",
]


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT).decode().strip()


def identity(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + "\n")


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


assert not (WORK / "preparation.json").exists() and not (WORK / CASE).exists()
assert git("status", "--porcelain", "--untracked-files=normal") == ""
runner = git("rev-parse", "HEAD")
allowed = [".github/workflows/ci.yml", "scripts/diagnose-macos-retained.py",
           "scripts/governance/build_policy.py", "scripts/macos_metal_probe.py"]
assert sorted(git("diff", "--name-only", PRODUCT, "HEAD").splitlines()) == allowed
previous = ROOT / "attestations/scientific-core/4bcd987/polarised-scientific-witnesses/files/out/cpu-polarised-scientific-witnesses"
for filename in ("run_owned.py", "owned_guard.py"):
    shutil.copyfile(previous / filename, WORK / filename)
source = ROOT / "tests/core/cie1931_observer_authority_test.cpp"
assert source.read_bytes() == subprocess.check_output(["git", "show", PRODUCT + ":" + key(source)], cwd=ROOT)
reference_base = "https://files.cie.co.at/Publications-datasets/"
for filename in ("CIE_xyz_1931_2deg.csv", "CIE_xyz_1931_2deg.csv_metadata.json"):
    with urllib.request.urlopen(reference_base + filename, timeout=30) as response:
        (WORK / filename).write_bytes(response.read())
raw = (WORK / "CIE_xyz_1931_2deg.csv").read_bytes()
assert hashlib.md5(raw).hexdigest() == "17cca777db64b17170f06f67ce9d3ab7"
assert hashlib.sha256(raw).hexdigest() == "fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1"
metadata = json.loads((WORK / "CIE_xyz_1931_2deg.csv_metadata.json").read_text())
assert metadata["identifier"]["identifier"] == "10.25039/CIE.DS.xvudnb9b"
assert {x["hashMethod"]: x["checksum"] for x in metadata["checksums"]} == {
    "md5": "17cca777db64b17170f06f67ce9d3ab7",
    "sha256": "fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1"}
table = {Decimal(r[0]): tuple(Decimal(x) for x in r[1:]) for r in csv.reader(io.StringIO(raw.decode()))}
block = source.read_text().split("kOfficialCieSamples = {{", 1)[1].split("}};", 1)[0]
samples = [tuple(Decimal(x.strip()) for x in line.split(",")) for line in re.findall(r"\{([^{}]+)\}", block)]
assert len(samples) == 22 and all(table[x[0]] == x[1:] for x in samples)
review_path = ROOT / "attestations/fresh-session-review/4bcd987/2026-10-14/review.json"
review = json.loads(review_path.read_text())
prior = review["current_live_byte_bindings"]["tested_artifacts"]["records"]["sirius_core_tests"]
core = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests"
assert identity(core) == {k: prior[k] for k in ("bytes", "sha256")}
accept_path = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation/files/out/gcc-revision-synchronisation/accepted-result.json"
accept = json.loads(accept_path.read_text())
assert accept["passed"] and accept["source_revision"] == PRODUCT
record = next(x for x in accept["current_version_recovery_records"] if x.get("artifact") == "sirius_core_tests")
archive = ROOT / record["recovery"]["path"]
assert identity(archive) == identity(core)
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
assert all(cache.get(k) == value for k, value in expected.items())
commands = json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
command = next(c for c in commands if Path(c["file"]) == source)
assert all(flag in command["command"] for flag in ("-std=c++2c", "-Werror", "-O3"))
assert "-ffast-math" not in command["command"]
registration = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake"
selected = []
for name in NAMES:
    matching = [line for line in registration.read_text().splitlines() if line.startswith("add_test(") and name in line]
    assert len(matching) == 1 and str(core) in matching[0] and "--gtest_filter=" + name in matching[0]
    selected += matching
inputs = [core, archive, cache_path, ROOT / "bin/linux-gcc/build.ninja", ROOT / "bin/linux-gcc/compile_commands.json",
          registration, review_path, accept_path, Path(cache["CMAKE_CTEST_COMMAND"]), Path(shutil.which("python3")).resolve(),
          source, ROOT / "src/sirius/core/cie1931_observer.h", ROOT / "src/sirius/core/spectral/blackbody.h",
          ROOT / "src/sirius/core/spectral/spectral_radiance.h", ROOT / "tests/operating_model.json"]
for glob in ("CTestTestfile.cmake", "*_tests.cmake", "*_include.cmake"):
    inputs += list((ROOT / "bin/linux-gcc").rglob(glob))
inputs = list(dict.fromkeys(p.resolve() for p in inputs))
local = [WORK / name for name in ("prepare.py", "run_cie.py", "run_owned.py", "owned_guard.py",
                                 "CIE_xyz_1931_2deg.csv", "CIE_xyz_1931_2deg.csv_metadata.json")]
pattern = "^(" + "|".join(re.escape(name) for name in NAMES) + ")$"
write("actions.json", {CASE: {
    "revision": runner, "kind": "check", "expected_exit": 0,
    "limits": {"seconds": 300, "rss_kib": 4 * 1024 * 1024, "reserve_kib": 2 * 1024 * 1024},
    "scope": "Exactly three original CPU CIE observer authority witnesses at unchanged current GCC core producer bytes; original official reference data and all assertions retained. No renderer, Vulkan dispatch, full P6, Mandatory, performance or release qualification claim.",
    "argv": [cache["CMAKE_CTEST_COMMAND"], "--test-dir", "bin/linux-gcc", "--output-on-failure", "--no-tests=error",
             "--output-junit", "{folder}/ctest.xml", "-R", pattern],
    "frozen_inputs": [key(p) for p in inputs + local] + [key(WORK / "preparation.json")],
}})
pins = {key(p): identity(p) for p in inputs + local + [WORK / "actions.json"]}
write("preparation.json", {
    "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "product_source_revision": PRODUCT, "runner_source_revision": runner, "diagnostic_only_difference": allowed,
    "executed": False, "scope": "Ready three-case original CPU host CIE action only; no product defect or complete P6/Mandatory claim.",
    "finite_acceptance_scope": {
        "official_reference": "Exactly the original 22 embedded official samples, each with x/y/z absolute tolerance 0.014; no bound at every wavelength is established.",
        "facades": "The existing spectral-bin loop verifies actual bin wavelengths and both host facade consumers.",
        "domain": "The original five outside-domain/nonfinite values decline to zero; 380 nm z-bar and 780 nm x-bar remain positive.",
    },
    "preparation_attempts": [{"url": "https://files.cie.co.at/CIE_xyz_1931_2deg.csv", "error": "HTTP 404; no child launched; official page now links Publications-datasets/"}],
    "reference": {"primary_page": "https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer",
                  "csv_url": reference_base + "CIE_xyz_1931_2deg.csv", "metadata_url": reference_base + "CIE_xyz_1931_2deg.csv_metadata.json",
                  "doi": "10.25039/CIE.DS.xvudnb9b", "original_publication": "CIE 018:2019, Table 6",
                  "csv_md5": "17cca777db64b17170f06f67ce9d3ab7",
                  "csv_sha256": "fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1",
                  "embedded_original_samples_exactly_joined": 22, "original_absolute_per_channel_assertion": 0.014,
                  "source_test_introduction": git("rev-parse", "831e89b")},
    "source_and_registration_join": {"test_source": key(source), "test_source_git_blob": git("rev-parse", PRODUCT + ":" + key(source)),
        "core_producer": identity(core), "core_recovery_archive": key(archive), "accepted_source_build_evidence": key(accept_path),
        "fresh_session_review": key(review_path), "cache_required": expected, "original_compile_command": command,
        "original_registration_file": key(registration), "original_three_registration_lines": selected, "expected_names": NAMES},
    "owner_reuse": {name: {"source": key(previous / name), **identity(WORK / name)} for name in ("run_owned.py", "owned_guard.py")},
    "pins": pins,
})
for path in local:
    if path.suffix == ".py":
        ast.parse(path.read_text(), filename=str(path))
print(json.dumps({"prepared": str(WORK), "executed": False, "core_producer": identity(core),
                  "official_samples_joined": 22, "pinned_inputs": len(pins),
                  "command": ["python3", "-B", key(WORK / "run_cie.py")]}, indent=2))
