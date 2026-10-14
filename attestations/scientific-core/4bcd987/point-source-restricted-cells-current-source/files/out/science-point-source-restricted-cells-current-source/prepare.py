"""Prepare three original CPU restricted-cell controls; no test/build/device work."""
from pathlib import Path
import ast
import copy
import difflib
import hashlib
import json
import re
import shutil
import subprocess

ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
COMMIT = "f6b5fbee9e2a3cc3e3be03d00ddd18b3103f4f44"
TAG = "evidence/issue-53-point-source-transfer-current-source"
BASE = ROOT / "attestations/scientific-core/4bcd987/point-source-transfer-current-source"
OLD = BASE / "files/out/science-point-source-transfer-current-source"
PRODUCER_BASE = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation"
PRODUCER = PRODUCER_BASE / "files/out/gcc-revision-synchronisation"
CORE_SHA = "c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58"
CASE = "original-three-cpu-restricted-cell-controls"
NAMES = [
    "PointSourceResponseTest.RestrictedCellsPreserveOriginalGaussianAndNormalization",
    "PointSourceResponseTest.RestrictedCellsOwnSharedEdgesAndKeepParentCutoff",
    "PointSourceResponseTest.RestrictedCellsDeclineInvalidOrUnrepresentedGeometry",
]
LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def key(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def identity(path):
    path = Path(path)
    with path.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"bytes": path.stat().st_size, "sha256": sha}


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def main():
    require(WORK == ROOT / "out/science-point-source-restricted-cells-current-source", "fixed new root")
    require(not (WORK / "preparation.json").exists(), "already sealed")
    require(git("rev-parse", "HEAD").decode().strip() == RUNNER and not git("status", "--porcelain").strip(), "clean pinned runner")
    require(git("rev-parse", TAG + "^{commit}").decode().strip() == COMMIT, "immutable accepted CPU tag")
    manifest = json.loads((BASE / "manifest.json").read_text())["files"]
    old_case = OLD / "original-six-cpu-point-source-transfer-controls"
    old_files = [BASE / "manifest.json", BASE / "result.json", OLD / "preparation.json", OLD / "actions.json",
                 OLD / "run_transfer.py", OLD / "run_owned.py", OLD / "owned_guard.py", OLD / "result.json",
                 old_case / "owner.json", old_case / "ctest.xml", old_case / "ctest.log"]
    historical = {}
    for path in old_files:
        raw = path.read_bytes()
        require(raw == git("show", COMMIT + ":" + key(path)), "immutable prior CPU bytes:" + key(path))
        relative = str(path.relative_to(BASE))
        if relative in manifest:
            require(identity(path) == {k: manifest[relative][k] for k in ("bytes", "sha256")}, "prior CPU manifest")
        historical[key(path)] = {**identity(path), "git_blob": git("rev-parse", COMMIT + ":" + key(path)).decode().strip()}
    old_result = json.loads((OLD / "result.json").read_text())
    old_owner = json.loads((old_case / "owner.json").read_text())
    require(old_result["passed"] and old_owner["passed"] and old_owner["child_exit"] == 0
            and old_owner["source_unchanged"] and old_owner["whole_inputs_unchanged"]
            and old_owner["remaining_owned_processes"] == [] and all(x["absent"] for x in old_owner["observed_births"]), "accepted mature CPU pattern")
    original = json.loads((OLD / "preparation.json").read_text())
    core = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests"
    core_recovery = PRODUCER_BASE / "prior-versions" / CORE_SHA
    require(identity(core) == identity(core_recovery) == {"bytes": 5199120, "sha256": CORE_SHA}, "exact core/recovery")
    authority = json.loads((PRODUCER / "current-product-authority.json").read_text())
    record = authority["actual_artifact_records"]["tested_artifacts"]["sirius_core_tests"]
    require(authority["passed"] and authority["source_revision"] == PRODUCT
            and all(record[k] == identity(core)[k] for k in ("bytes", "sha256")), "current accepted core authority")
    accepted = json.loads((PRODUCER / "accepted-result.json").read_text())
    owner = json.loads((PRODUCER / "build/owner.json").read_text())
    require(accepted["passed"] and owner["passed"] and owner["child_exit"] == 0
            and owner["source_unchanged"] and owner["whole_inputs_unchanged"] and owner["remaining_owned_processes"] == [], "existing build evidence")
    before = json.loads((PRODUCER / "build/inputs-before.json").read_text())
    after = json.loads((PRODUCER / "build/inputs-after.json").read_text())
    source = (ROOT / "tests/core/point_source_response_test.cpp").read_text()
    require(["PointSourceResponseTest." + n for n in re.findall(r"TEST\(PointSourceResponseTest, (RestrictedCells\w+)\)", source)] == NAMES, "only original three source names")
    require("constexpr int radii = 1024, angles = 8;" in source and "EXPECT_NEAR(integral, 1.0, 1e-6);" in source, "original finite quadrature")
    model = json.loads((ROOT / "tests/operating_model.json").read_text())
    p3 = next(x for x in model["acceptance_criteria"] if x["id"] == "P3")
    require(set(NAMES).issubset({x.get("name") for x in p3["evidence"]}), "original model names")
    reg = (ROOT / "bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake").read_text().splitlines()
    registrations = {}
    for name in NAMES:
        add = next(x for x in reg if x.startswith("add_test(") and "[=[" + name + "]=]" in x)
        properties = next(x for x in reg if x.startswith("set_tests_properties(") and "[=[" + name + "]=]" in x)
        require(str(core) in add and "--gtest_filter=" + name in add and "DISABLED" not in properties and "TIMEOUT" not in properties, "original enabled registration")
        registrations[name] = {"add_test": add, "properties": properties}
    command = next(x for x in json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
                   if x["file"].endswith("/tests/core/point_source_response_test.cpp"))
    require(all(x in command["command"] for x in ("SIRIUS_CONTRACT_MODE=2", "SIRIUS_RELEASE_RESOURCE_LOCKED=1", "-O3 -DNDEBUG", "-std=c++2c", "-Werror", "tests/support/fail_on_skip.h")), "strict producer recipe")
    link_line = next(x for x in (ROOT / "bin/linux-gcc/build.ninja").read_text().splitlines()
                     if x.startswith("build tests/core/sirius_core_tests ") and "CXX_EXECUTABLE_LINKER__sirius_core_tests_Release" in x)
    require("point_source_response_test.cpp.o" in link_line and "src/sirius/core/libsirius_core.a" in link_line, "actual core link topology")
    pending = [ROOT / n for n in ("tests/core/point_source_response_test.cpp", "tests/support/fail_on_skip.h",
        "tests/core/CMakeLists.txt", "tests/CMakeLists.txt", "cmake/sirius_testing.cmake", "tests/operating_model.json")]
    paths = set()
    while pending:
        path = pending.pop().resolve()
        if path in paths:
            continue
        paths.add(path)
        if path.suffix in (".cpp", ".h", ".c"):
            for name in re.findall(r'^\s*#\s*include\s*"([^\"]+)"', path.read_text(), re.M):
                target = next((p.resolve() for p in (path.parent / name, ROOT / "src" / name, ROOT / "tests" / name) if p.is_file()), None)
                if target is not None and target.is_relative_to(ROOT) and not target.is_relative_to(ROOT / "bin"):
                    pending.append(target)
    snapshots = {}
    for path in sorted(paths):
        name = key(path)
        raw = git("show", PRODUCT + ":" + name)
        require(raw == path.read_bytes() and before[name] == after[name]
                and all(before[name][k] == identity(path)[k] for k in ("bytes", "sha256")), "original source/producer continuity:" + name)
        dest = WORK / "reference" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
        snapshots[name] = {**identity(path), "revision": PRODUCT, "git_blob": git("rev-parse", PRODUCT + ":" + name).decode().strip(),
                           "snapshot": key(dest), "build_before_after_equal": True}
    helpers = {}
    for filename, sha in (("run_owned.py", "b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d"),
                          ("owned_guard.py", "747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056")):
        require(identity(OLD / filename)["sha256"] == sha, "mature CPU helper identity")
        shutil.copy2(OLD / filename, WORK / filename)
        helpers[filename] = {"source": key(OLD / filename), **identity(OLD / filename), "copied_byte_exact": True}
    old_reader = (OLD / "run_transfer.py").read_text()
    reader = old_reader.replace('"""One frozen six-case original CPU point-source transfer action."""',
                               '"""One frozen three-case original CPU restricted-cell response action."""')
    reader = reader.replace('CASE = "original-six-cpu-point-source-transfer-controls"', 'CASE = "' + CASE + '"')
    start, end = reader.index("NAMES = {"), reader.index("\nLIMITS = ")
    reader = reader[:start] + "NAMES = {\n" + "".join("    " + repr(n) + ",\n" for n in NAMES) + "}" + reader[end:]
    require(reader.count("len(NAMES) == 6") == 1 and reader.count("len(cases) == 6") == 1, "existing exact-count guards")
    reader = reader.replace("len(NAMES) == 6", "len(NAMES) == 3").replace("len(cases) == 6", "len(cases) == 3")
    ast.parse(reader)
    (WORK / "run_restricted_cells.py").write_text(reader)
    (WORK / "reader-selection.patch").write_text("".join(difflib.unified_diff(old_reader.splitlines(True), reader.splitlines(True),
        fromfile="accepted/run_transfer.py", tofile="restricted-cells/run_restricted_cells.py")))
    context = {"schema_version": 1, "prepared_only": True, "executed": False,
        "product_source_revision": PRODUCT, "runner_source_revision": RUNNER,
        "diagnostic_only_difference": git("diff", "--name-only", PRODUCT, "HEAD").decode().splitlines(),
        "expected_names": NAMES, "limits": LIMITS, "scope": "Three original PointSourceResponseTest restricted-cell CPU controls only, finite P3 with P2 normalization overlap; no fullP2/P3/P5/science/Mandatory/domain/release or performance claim.",
        "original_assertions_changed": False,
        "finite_acceptance_scope": {
            NAMES[0]: "Original skew affine map, four cells, five density samples and direct truncated-Gaussian formula8e-14*max(expected,1); independent original-coordinate1024x8 polar midpoint quadrature with Jacobian applied once, integral1 within1e-6.",
            NAMES[1]: "Original shared-edge half-open ownership, inclusive parent4-sigma cutoff and excluded corner/exterior; skew inverse-rounding NeedsRefinement, exact original-root ownership, rotated-corner and narrow-cell outward query-radius coverage fixtures unchanged.",
            NAMES[2]: "Original out-of-parent/empty/singular/lost-coordinate domains decline; nonfinite density/root inputs return original InvalidInput with no fabricated density."},
        "reference": {"original_source_snapshots": snapshots,
            "independent_reference": "Test-local RestrictedSourceOffset uses direct forward affine multiplication. Gaussian density is defined directly with absolute map determinant and -expm1(-8); independent polar quadrature is in original coordinates and applies the source-area Jacobian once. No production inverse/normalizer supplies expected values.",
            "composition_boundary": "Restriction partitions the original response support without renormalizing child cells. Source algebra and finite fixtures are not universal response/integration or all-map proofs; original tests/tolerances remain sole numerical acceptance."},
        "source_and_registration_join": {"test_sources": ["tests/core/point_source_response_test.cpp"],
            "production_headers": ["src/sirius/core/point_source_response.h"], "source_snapshot_count": len(snapshots),
            "core_producer": identity(core), "core_recovery_archive": key(core_recovery),
            "accepted_source_build_evidence": key(PRODUCER / "accepted-result.json"), "existing_build_owner": key(PRODUCER / "build/owner.json"),
            "source_before_after": [key(PRODUCER / "build/inputs-before.json"), key(PRODUCER / "build/inputs-after.json")],
            "original_compile_commands": [command], "original_link_ninja_lines": [link_line], "registrations": registrations,
            "runtime_resources_bound": {}, "external_image_or_device_resources_required": False, "selected_core_cases": NAMES,
            "selected_render_cases": [], "original_test_TIMEOUT_added": False,
            "compiler_provenance_limit": "Reuse existing accepted current core/source authority; current recipe/runtime bytes frozen, no new build or retroactive compiler-byte attestation."},
        "owner_reuse": {"archive": key(BASE), "revision": COMMIT, "tag": TAG,
            "run_owned.py": helpers["run_owned.py"], "owned_guard.py": helpers["owned_guard.py"],
            "reader_baseline": key(OLD / "run_transfer.py"), "reader_change": "Only component/selection/cardinality3 instead of6; original source/pin/owner/XML/TMPDIR checks and300s/4GiB/2GiB bounds unchanged.",
            "immutable_prior_files": historical},
        "private_tmpdir": {"path": key(WORK / "temporary"), "mode": "0700", "images_expected": False,
            "diagnostics_only": True, "natural_residuals_preserved": True},
        "xml_scope": "One CTest JUnit with exactly3 original registered cases and zero failures/errors/skips/disabled; original per-case GTest output in verbose CTest/stdout, GTEST_OUTPUT cleared to prevent last-child overwrite. No synthetic XML.",
        "scheduling": "Root-only one-shot CPU action after independent review; no dependencies on old/fresh compiler or host workspaces and no preparation-time launches."}
    write(WORK / "context.json", context)
    assets = set(paths) | set(old_files) | {core, core_recovery, PRODUCER_BASE / "manifest.json"}
    assets.update(PRODUCER / n for n in ("accepted-result.json", "current-product-authority.json", "source-bindings.json", "build/owner.json", "build/inputs-before.json", "build/inputs-after.json"))
    assets.update(ROOT / "bin/linux-gcc" / n for n in ("CMakeCache.txt", "build.ninja", "compile_commands.json", "CTestTestfile.cmake", "tests/core/CTestTestfile.cmake", "tests/core/sirius_core_tests[1]_include.cmake", "tests/core/sirius_core_tests[1]_tests.cmake"))
    assets.update(Path(n) for n in original["pins"] if n.startswith(("/usr/lib/", "/lib/")) or n in ("/usr/bin/ctest", "/usr/bin/python3.12"))
    assets.update(p for p in WORK.rglob("*") if p.is_file())
    frozen = sorted({key(p) for p in assets} | {key(WORK / name) for name in ("preparation.json", "actions.json")})
    action = {"revision": RUNNER, "kind": "check", "expected_exit": 0, "limits": LIMITS,
        "scope": context["scope"], "frozen_inputs": frozen,
        "argv": ["/usr/bin/ctest", "--test-dir", "bin/linux-gcc", "-V", "--output-on-failure", "--no-tests=error", "--output-junit", "{folder}/ctest.xml", "--output-log", "{folder}/ctest.log", "-R", "^(" + "|".join(re.escape(n) for n in sorted(NAMES)) + ")$"]}
    write(WORK / "actions.json", {CASE: action})
    assets.add(WORK / "actions.json")
    context["pins"] = {key(p): identity(p) for p in sorted(assets)}
    write(WORK / "preparation.json", context)
    require(set(frozen) == set(context["pins"]) | {key(WORK / "preparation.json")}, "complete owner closure")
    require(not any(n.startswith("out/") and not n.startswith(key(WORK) + "/") for n in context["pins"]), "no other live scratch inputs")
    print(json.dumps({"prepared_only": True, "runtime_actions": 0, "cases": 3, "pins": len(context["pins"]),
        "source_snapshots": len(snapshots), "core": identity(core), "preparation": identity(WORK / "preparation.json"),
        "reader": identity(WORK / "run_restricted_cells.py"), "actions": identity(WORK / "actions.json")}, indent=2))


if __name__ == "__main__":
    main()
