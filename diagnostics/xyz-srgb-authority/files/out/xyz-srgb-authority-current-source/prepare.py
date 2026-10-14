#!/usr/bin/env python3
"""Prepare six original CPU XYZ/sRGB authority cases; execute no tests."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
CASE = "original-six-host-colour-controls"
NAMES = [
    "XyzSrgbAuthority.MatchesIndependentPrimaryWhiteDerivation",
    "XyzSrgbAuthority.PreservesExtendedGamutAndPropagatesNonfiniteInputs",
    "XyzSrgbAuthority.BlackbodyColourFacadeDelegatesToTheExactTransform",
    "SrgbTransferAuthority.CurveMatchesIndependentIecOracleAcrossBothBranches",
    "SrgbTransferAuthority.EightBitQuantisationClipsAndDeclinesNonfiniteInputs",
    "SrgbTransferAuthority.SpectralFacadeDelegatesWithoutChangingClippingSemantics",
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


assert not (WORK / "preparation.json").exists() and not (WORK / CASE).exists()
assert git("status", "--porcelain", "--untracked-files=normal") == ""
runner = git("rev-parse", "HEAD")
allowed = [".github/workflows/ci.yml", "scripts/diagnose-macos-retained.py",
           "scripts/governance/build_policy.py", "scripts/macos_metal_probe.py"]
assert sorted(git("diff", "--name-only", PRODUCT, "HEAD").splitlines()) == allowed
previous = ROOT / "attestations/scientific-core/4bcd987/cie1931-authority-current-source/files/out/cie1931-authority-current-source"
for filename in ("run_owned.py", "owned_guard.py"):
    shutil.copyfile(previous / filename, WORK / filename)
wrapper = (previous / "run_cie.py").read_text()
wrapper = wrapper.replace("three-case CPU CIE", "six-case CPU XYZ/sRGB")
wrapper = wrapper.replace("original-three-cie1931-controls", CASE)
names_literal = "NAMES = {\n" + "".join("    " + repr(name) + ",\n" for name in NAMES) + "}"
wrapper = re.sub(r"NAMES = \{.*?\n\}", names_literal, wrapper, count=1, flags=re.S)
wrapper = wrapper.replace("Three original CPU CIE authority", "Six original CPU XYZ/sRGB authority")
wrapper = wrapper.replace("len(cases) == 3", "len(cases) == 6")
wrapper = wrapper.replace("three original cases", "six original cases")
assert "Cie1931" not in wrapper and "CPU CIE" not in wrapper
(WORK / "run_colour.py").write_text(wrapper)
sources = [ROOT / "tests/core/xyz_srgb_authority_test.cpp", ROOT / "tests/core/srgb_transfer_authority_test.cpp"]
(WORK / "reference").mkdir()
for path in sources:
    assert path.read_bytes() == subprocess.check_output(["git", "show", PRODUCT + ":" + key(path)], cwd=ROOT)
    shutil.copyfile(path, WORK / "reference" / path.name)
review_path = ROOT / "attestations/fresh-session-review/4bcd987/2026-10-14/review.json"
review = json.loads(review_path.read_text())
p6 = next(x for x in review["model_requirement_ledger"] if x["id"] == "P6")
assert set(NAMES) <= set(p6["not_proved_by_reviewed_current_bundles"])
assert not any(set(NAMES) & set(x["passed_names"]) for x in review["raw_junit"])
core = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests"
prior = review["current_live_byte_bindings"]["tested_artifacts"]["records"]["sirius_core_tests"]
assert identity(core) == {k: prior[k] for k in ("bytes", "sha256")}
accept_path = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation/files/out/gcc-revision-synchronisation/accepted-result.json"
accept = json.loads(accept_path.read_text())
assert accept["passed"] and accept["source_revision"] == PRODUCT
record = next(x for x in accept["current_version_recovery_records"] if x.get("artifact") == "sirius_core_tests")
recovery = ROOT / record["recovery"]["path"]
assert identity(recovery) == identity(core)
cie_disposition = ROOT / "attestations/scientific-core/4bcd987/cie1931-authority-current-source/disposition.json"
cie = json.loads(cie_disposition.read_text())
assert cie["evidence_joins"]["core_producer"] == identity(core)
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
selected_commands = [next(c for c in commands if Path(c["file"]) == source) for source in sources]
for command in selected_commands:
    assert all(flag in command["command"] for flag in ("-std=c++2c", "-Werror", "-O3"))
    assert "-ffast-math" not in command["command"]
registration = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake"
selected = []
for name in NAMES:
    matching = [line for line in registration.read_text().splitlines() if line.startswith("add_test(") and name in line]
    assert len(matching) == 1 and str(core) in matching[0] and "--gtest_filter=" + name in matching[0]
    selected += matching
inputs = [core, recovery, cache_path, ROOT / "bin/linux-gcc/build.ninja", ROOT / "bin/linux-gcc/compile_commands.json",
          registration, review_path, accept_path, cie_disposition, Path(cache["CMAKE_CTEST_COMMAND"]),
          Path(shutil.which("python3")).resolve(), *sources, ROOT / "src/sirius/core/xyz_srgb.h",
          ROOT / "src/sirius/core/srgb_transfer.h", ROOT / "src/sirius/core/constants.h",
          ROOT / "src/sirius/core/spectral/blackbody.h", ROOT / "tests/operating_model.json"]
for glob in ("CTestTestfile.cmake", "*_tests.cmake", "*_include.cmake"):
    inputs += list((ROOT / "bin/linux-gcc").rglob(glob))
inputs = list(dict.fromkeys(path.resolve() for path in inputs))
local = [WORK / name for name in ("prepare.py", "run_colour.py", "run_owned.py", "owned_guard.py")]
local += [WORK / "reference" / path.name for path in sources]
write("actions.json", {CASE: {
    "revision": runner, "kind": "check", "expected_exit": 0,
    "limits": {"seconds": 300, "rss_kib": 4 * 1024 * 1024, "reserve_kib": 2 * 1024 * 1024},
    "scope": "Exactly six original CPU XYZ/sRGB authority witnesses on unchanged GCC core bytes; no CIE repeat, changed assertion, renderer, GPU dispatch, full P6, Mandatory, performance or release claim.",
    "argv": [cache["CMAKE_CTEST_COMMAND"], "--test-dir", "bin/linux-gcc", "--output-on-failure", "--no-tests=error",
             "--output-junit", "{folder}/ctest.xml", "-R", "^(" + "|".join(re.escape(name) for name in NAMES) + ")$"],
    "frozen_inputs": [key(path) for path in inputs + local] + [key(WORK / "preparation.json")],
}})
write("preparation.json", {
    "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "product_source_revision": PRODUCT, "runner_source_revision": runner, "diagnostic_only_difference": allowed,
    "executed": False, "scope": "One six-case original CPU colour-authority action; no full P6 or Mandatory qualification claim.",
    "finite_acceptance_scope": {
        "xyz_matrix": "Exactly six original XYZ vectors against independent primary/white reconstruction, with the existing 64*epsilon_double*max(1,abs(expected channels)) tolerance; no all-input bound.",
        "xyz_domain": "Original green-axis extended-gamut inequalities, NaN-channel propagation, and nonfinite results for one infinite-X input.",
        "xyz_facade": "One original float XYZ vector (0.25f,0.4f,0.1f), judged by original EXPECT_FLOAT_EQ for each blackbody facade channel.",
        "srgb_curve": "Ten original linear inputs, including the breakpoint and nextafter neighbours, with existing 16*epsilon_double*max(1,abs(expected)) tolerance.",
        "srgb_quantisation": "Seven original exact eight-bit expectations and NaN/positive-infinity/negative-infinity refusals.",
        "srgb_facade": "Seven original scalar float inputs, one three-channel clipping vector, and the existing NaN scalar check; original EXPECT_FLOAT_EQ semantics.",
    },
    "reference": {
        "kind": "Reuse the original independently reconstructed IEC sRGB primaries/CIE D65 transform and independent IEC piecewise transfer oracle; no new reference or numerical predicate.",
        "xyz_source_function": "IndependentlyDerivedXyzD65ToLinearSrgb",
        "primary_xy": {"red": ["0.640", "0.330"], "green": ["0.300", "0.600"], "blue": ["0.150", "0.060"], "D65": ["0.3127", "0.3290"]},
        "xyz_derivation": "Construct unit-luminance primary columns, solve their scales against D65, scale the RGB-to-XYZ columns, then independently invert. It does not copy the production rational matrix.",
        "xyz_six_vectors": [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"], ["0.9504559270516716", "1", "1.0890577507598784"], ["0.25", "0.4", "0.1"], ["-0.1", "0.2", "1.5"]],
        "srgb_source_function": "IndependentIecSrgbEncode",
        "srgb_definition": "linear<=0.0031308 ? 12.92*linear : 1.055*pow(linear,1/2.4)-0.055",
        "original_test_sources": {key(path): {"git_blob": git("rev-parse", PRODUCT + ":" + key(path)), **identity(path)} for path in sources},
        "source_snapshots": {key(WORK / "reference" / path.name): identity(WORK / "reference" / path.name) for path in sources},
    },
    "source_and_registration_join": {"test_sources": [key(path) for path in sources], "core_producer": identity(core),
        "core_recovery_archive": key(recovery), "accepted_source_build_evidence": key(accept_path),
        "completed_cie_core_join_reused": key(cie_disposition), "fresh_session_review": key(review_path),
        "cache_required": expected, "original_compile_commands": selected_commands,
        "original_registration_file": key(registration), "original_six_registration_lines": selected, "expected_names": NAMES},
    "ledger_scope": "All six names are explicitly unproved in the existing current-source ledger and absent from its reviewed JUnit passes; no global absence of all historical executions is inferred. Completed CIE evidence is retained and not rerun.",
    "owner_reuse": {name: {"source": key(previous / name), **identity(WORK / name)} for name in ("run_owned.py", "owned_guard.py")},
    "pins": {key(path): identity(path) for path in inputs + local + [WORK / "actions.json"]},
})
for path in local:
    if path.suffix == ".py":
        ast.parse(path.read_text(), filename=str(path))
print(json.dumps({"prepared_without_execution": True, "cases": len(NAMES), "core_producer": identity(core),
                  "preparation": identity(WORK / "preparation.json"), "runner": runner,
                  "command": ["python3", "-B", key(WORK / "run_colour.py")]}, indent=2))
