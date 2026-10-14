"""Seal six original CPU transfer controls; never build or execute a test."""
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
OLD_BUNDLE = ROOT / "attestations/scientific-core/4bcd987/p4-cpu-radiation-current-source"
OLD = OLD_BUNDLE / "files/out/p4-cpu-radiation-current-source"
PRODUCER_BUNDLE = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation"
PRODUCER = PRODUCER_BUNDLE / "files/out/gcc-revision-synchronisation"
CORE_SHA = "c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58"
CASE = "original-six-cpu-point-source-transfer-controls"
NAMES = [
    "PointSourceTransfer.ReferenceFrequencyPreservesCatalogueColourAndFlux",
    "PointSourceTransfer.ShiftedSpectrumMatchesIndependentFixedBandPlanckLaw",
    "PointSourceTransfer.RejectsReciprocalExtraGainAndShiftedNormalizerSubstitutions",
    "PointSourceTransfer.MovingObserverUsesBandGainAndAngularAreaExactlyOnce",
    "PointSourceTransfer.ShiftedTemperatureHasNoCatalogueLookupClamp",
    "PointSourceTransfer.InvalidOrOverflowingTransferHasNoPublishedValue",
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
    require(WORK == ROOT / "out/science-point-source-transfer-current-source", "fixed fresh packet")
    require(not (WORK / "preparation.json").exists(), "already sealed")
    require(git("rev-parse", "HEAD").decode().strip() == RUNNER and not git("status", "--porcelain").strip(), "clean pinned runner")
    core = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests"
    recovery = PRODUCER_BUNDLE / "prior-versions" / CORE_SHA
    require(identity(core) == identity(recovery) == {"bytes": 5199120, "sha256": CORE_SHA}, "exact current core/recovery")
    authority = json.loads((PRODUCER / "current-product-authority.json").read_text())
    require(authority["passed"] and authority["source_revision"] == PRODUCT, "accepted producer source")
    record = authority["actual_artifact_records"]["tested_artifacts"]["sirius_core_tests"]
    require(record["path"] == "tests/core/sirius_core_tests" and all(record[k] == identity(core)[k] for k in ("bytes", "sha256")), "actual accepted core record")
    accepted = json.loads((PRODUCER / "accepted-result.json").read_text())
    owner = json.loads((PRODUCER / "build/owner.json").read_text())
    require(accepted["passed"] and owner["passed"] and owner["child_exit"] == 0 and owner["source_unchanged"]
            and owner["whole_inputs_unchanged"] and owner["remaining_owned_processes"] == [], "completed existing build authority")
    before = json.loads((PRODUCER / "build/inputs-before.json").read_text())
    after = json.loads((PRODUCER / "build/inputs-after.json").read_text())

    source = (ROOT / "tests/core/point_source_transfer_test.cpp").read_text()
    require(["PointSourceTransfer." + name for name in re.findall(r"TEST\(PointSourceTransfer, (\w+)\)", source)] == NAMES, "exact six original source cases")
    model = json.loads((ROOT / "tests/operating_model.json").read_text())
    p3 = next(row for row in model["acceptance_criteria"] if row["id"] == "P3")
    require(set(NAMES).issubset({row.get("name") for row in p3["evidence"]}), "original P3 names")
    registered = (ROOT / "bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake").read_text().splitlines()
    registrations = {}
    for name in NAMES:
        add = next(line for line in registered if line.startswith("add_test(") and "[=[" + name + "]=]" in line)
        properties = next(line for line in registered if line.startswith("set_tests_properties(") and "[=[" + name + "]=]" in line)
        require(str(core) in add and "--gtest_filter=" + name in add and "DISABLED" not in properties and "TIMEOUT" not in properties, "original enabled registration")
        registrations[name] = {"add_test": add, "properties": properties}
    command = next(row for row in json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
                   if row["file"].endswith("/tests/core/point_source_transfer_test.cpp"))
    require(all(token in command["command"] for token in ("SIRIUS_CONTRACT_MODE=2", "SIRIUS_RELEASE_RESOURCE_LOCKED=1", "-O3 -DNDEBUG", "-std=c++2c", "-Werror", "tests/support/fail_on_skip.h")), "original strict compilation flags")
    link_line = next(line for line in (ROOT / "bin/linux-gcc/build.ninja").read_text().splitlines()
                     if line.startswith("build tests/core/sirius_core_tests ") and "CXX_EXECUTABLE_LINKER__sirius_core_tests_Release" in line)
    require("point_source_transfer_test.cpp.o" in link_line and "src/sirius/core/libsirius_core.a" in link_line, "actual core topology")

    seeds = ["tests/core/point_source_transfer_test.cpp", "tests/support/point_source_band_reference.h",
             "tests/support/fail_on_skip.h", "tests/core/CMakeLists.txt", "tests/CMakeLists.txt",
             "cmake/sirius_testing.cmake", "tests/operating_model.json"]
    pending = [ROOT / name for name in seeds]
    paths = set()
    while pending:
        path = pending.pop().resolve()
        if path in paths:
            continue
        paths.add(path)
        if path.suffix in (".cpp", ".h", ".c"):
            for name in re.findall(r'^\s*#\s*include\s*"([^\"]+)"', path.read_text(), re.M):
                candidate = next((p.resolve() for p in (path.parent / name, ROOT / "src" / name, ROOT / "tests" / name) if p.is_file()), None)
                if candidate is not None and candidate.is_relative_to(ROOT) and not candidate.is_relative_to(ROOT / "bin"):
                    pending.append(candidate)
    snapshots = {}
    for path in sorted(paths):
        name = key(path)
        raw = git("show", PRODUCT + ":" + name)
        require(path.read_bytes() == raw and before[name] == after[name], "source/build continuity:" + name)
        require(all(before[name][k] == identity(path)[k] for k in ("bytes", "sha256")), "original build input bytes:" + name)
        dest = WORK / "reference" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
        snapshots[name] = {**identity(path), "git_blob": git("rev-parse", PRODUCT + ":" + name).decode().strip(),
                           "revision": PRODUCT, "snapshot": key(dest), "build_before_after_equal": True}

    old_context = json.loads((OLD / "preparation.json").read_text())
    helpers = {}
    for filename, sha in (("run_owned.py", "b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d"),
                          ("owned_guard.py", "747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056")):
        require(identity(OLD / filename)["sha256"] == sha, "mature CPU owner/guard identity")
        shutil.copy2(OLD / filename, WORK / filename)
        helpers[filename] = {"source": key(OLD / filename), **identity(OLD / filename), "copied_byte_exact": True}
    original_reader = (OLD / "run_radiation.py").read_text()
    reader = original_reader.replace('"""One frozen nine-case original CPU radiation action; no build/device/frame work."""',
                                     '"""One frozen six-case original CPU point-source transfer action."""')
    reader = reader.replace('CASE = "original-nine-cpu-radiation-controls"', 'CASE = "' + CASE + '"')
    start, end = reader.index("NAMES = {"), reader.index("\nLIMITS = ")
    reader = reader[:start] + "NAMES = {\n" + "".join("    " + repr(n) + ",\n" for n in NAMES) + "}" + reader[end:]
    reader = reader.replace('"qualification_claimed": False, "p4_complete": False,',
                            '"qualification_claimed": False, "p2_complete": False, "p3_complete": False, "p5_complete": False,')
    require(reader.count("len(NAMES) == 9") == 1 and reader.count("len(cases) == 9") == 1, "original cardinality checks")
    reader = reader.replace("len(NAMES) == 9", "len(NAMES) == 6").replace("len(cases) == 9", "len(cases) == 6")
    ast.parse(reader)
    (WORK / "run_transfer.py").write_text(reader)
    (WORK / "reader-selection.patch").write_text("".join(difflib.unified_diff(original_reader.splitlines(True), reader.splitlines(True),
         fromfile="preserved/run_radiation.py", tofile="point-source/run_transfer.py")))
    nist = OLD / "reference/nist-2022-allascii.txt"
    require(identity(nist)["sha256"] == old_context["reference"]["nist_2022"]["sha256"], "preserved primary constants")
    scope = {
        NAMES[0]: "Original eight catalogue temperatures100..1e6K; reference frequency colour/linear flux-density agreement, adjacent-unit-frequency continuity and original finite-difference smoothness checks.",
        NAMES[1]: "Original3temperatures3000/6500/20000K times4frequency factors. Independent long-double g^5 B_lambda(g lambda,T) reference; original32-bin band response,8192 Simpson/direct identity and4096-to8192 refinement bounds; original2e-6 colour budget.",
        NAMES[2]: "Original6500K,g=.5/2 fixtures reject reciprocal frequency, extra g^4, shifted normalization and truncating the source band, by original100-times-tolerance separation.",
        NAMES[3]: "Original beta=.6/gamma1.25, three source directions and both map parities; independent Lorentz g/area relation and2048-interval two-node quadrature verify band gain times1/absolute angular Jacobian exactly once.",
        NAMES[4]: "Original three shifted-temperature fixtures, bright-source recovery below binary32, represented hot and underflow cases; unchanged finite sample assertions.",
        NAMES[5]: "Original nonfinite/negative/zero input and arithmetic overflow/underflow refusals; represented scaling, dark source and no partial published values.",
    }
    context = {"schema_version": 1, "prepared_only": True, "executed": False,
        "product_source_revision": PRODUCT, "runner_source_revision": RUNNER,
        "diagnostic_only_difference": git("diff", "--name-only", PRODUCT, "HEAD").decode().splitlines(),
        "expected_names": NAMES, "scope": "Six unchanged original PointSourceTransfer CPU controls only, finite P3 transfer with P2/P5 overlap; no fullP2/P3/P5/science/Mandatory/platform/release or performance qualification.",
        "finite_acceptance_scope": scope, "original_assertions_changed": False, "limits": LIMITS,
        "reference": {"original_source_snapshots": snapshots, "independent_reference": "tests/support/point_source_band_reference.h reimplements defining Planck equation with exact SI h,c,k and long-double arithmetic; no call to production Planck/colour/transfer helpers. Original matching-function coefficients and rational XYZ/RGB matrix express the same defined fits; this is not an independent all-wavelength CIE-table accuracy claim.",
            "derivation_scope": "g^5 B_lambda(g lambda,T)=B_lambda(lambda,gT) at fixed observed band; normalize by original source-temperature catalogue factor, then multiply catalogue flux and detector density once. Independent Lorentz angular area is distinct from fixed-band gain; no additional bolometric g^4.",
            "nist_2022": copy.deepcopy(old_context["reference"]["nist_2022"]), "nist_preserved_path": key(nist),
            "nist_scope": "Previously preserved primary coefficients only; no new retrieval or numerical gate. Original executable references and assertions are unchanged.",
            "source_composition_boundary": "Finite source composition review and original sample/quadrature witnesses, not exhaustive or machine-checked proof."},
        "source_and_registration_join": {"test_sources": ["tests/core/point_source_transfer_test.cpp"],
            "independent_reference_source": "tests/support/point_source_band_reference.h", "source_snapshot_count": len(snapshots),
            "core_producer": identity(core), "core_recovery_archive": key(recovery),
            "accepted_source_build_evidence": key(PRODUCER / "accepted-result.json"), "existing_build_owner": key(PRODUCER / "build/owner.json"),
            "source_before_after": [key(PRODUCER / "build/inputs-before.json"), key(PRODUCER / "build/inputs-after.json")],
            "original_compile_commands": [command], "original_link_ninja_lines": [link_line], "registrations": registrations,
            "runtime_resources_bound": {}, "external_image_or_device_resources_required": False,
            "selected_core_cases": NAMES, "selected_render_cases": [], "original_test_TIMEOUT_added": False,
            "compiler_provenance_limit": "Reuse accepted completed current producer/source authority; current recipe/tool bytes pinned now, no new build or retroactive historical compiler-byte attestation."},
        "owner_reuse": {"archive": key(OLD_BUNDLE), "run_owned.py": helpers["run_owned.py"], "owned_guard.py": helpers["owned_guard.py"],
            "reader_baseline": key(OLD / "run_radiation.py"), "reader_change": "Only component names, selected original six names/cardinality and nonqualification metadata; original checks and300s/4GiB/2GiB bounds unchanged."},
        "private_tmpdir": {"path": key(WORK / "temporary"), "mode": "0700", "images_expected": False,
            "diagnostics_only": True, "natural_residuals_preserved": True},
        "xml_scope": "One CTest JUnit contains exactly all6 original registered cases with zero failures/errors/skips/disabled. Original per-case GTest output is captured in verbose CTest/stdout; GTEST_OUTPUT stays cleared to avoid overwriting one file across6 distinct children. No synthetic GTest XML.",
        "scheduling": "Root only, after current scene2 retirement; no automatic launches during preparation. Parity3/4 and candidate compilation remain separate unlaunched actions."}
    write(WORK / "context.json", context)
    assets = set(paths) | {core, recovery, nist, OLD_BUNDLE / "manifest.json", OLD / "preparation.json", OLD / "run_radiation.py", OLD / "run_owned.py", OLD / "owned_guard.py"}
    assets.update(PRODUCER / name for name in ("accepted-result.json", "current-product-authority.json", "source-bindings.json", "build/owner.json", "build/inputs-before.json", "build/inputs-after.json"))
    assets.add(PRODUCER_BUNDLE / "manifest.json")
    assets.update(ROOT / "bin/linux-gcc" / name for name in ("CMakeCache.txt", "build.ninja", "compile_commands.json", "CTestTestfile.cmake", "tests/core/CTestTestfile.cmake", "tests/core/sirius_core_tests[1]_include.cmake", "tests/core/sirius_core_tests[1]_tests.cmake"))
    # Reuse known CPU controller/runtime library inventory, without ldd/tool/device probes.
    assets.update(Path(name) for name in old_context["pins"] if name.startswith(("/usr/lib/", "/lib/")) or name in ("/usr/bin/ctest", "/usr/bin/python3.12"))
    assets.update(p for p in WORK.rglob("*") if p.is_file())
    frozen = sorted({key(p) for p in assets} | {key(WORK / name) for name in ("preparation.json", "actions.json")})
    action = {"revision": RUNNER, "kind": "check", "expected_exit": 0, "limits": LIMITS,
        "scope": context["scope"], "frozen_inputs": frozen,
        "argv": ["/usr/bin/ctest", "--test-dir", "bin/linux-gcc", "-V", "--output-on-failure", "--no-tests=error",
                 "--output-junit", "{folder}/ctest.xml", "--output-log", "{folder}/ctest.log", "-R", "^(" + "|".join(re.escape(name) for name in sorted(NAMES)) + ")$"]}
    write(WORK / "actions.json", {CASE: action})
    assets.add(WORK / "actions.json")
    context["pins"] = {key(p): identity(p) for p in sorted(assets)}
    write(WORK / "preparation.json", context)
    require(set(frozen) == set(context["pins"]) | {key(WORK / "preparation.json")}, "complete owner frozen closure")
    require(not any(name.startswith("out/") and not name.startswith(key(WORK) + "/") for name in context["pins"]), "no other scratch input")
    print(json.dumps({"prepared_only": True, "runtime_actions": 0, "cases": 6, "pins": len(context["pins"]),
          "source_snapshot_count": len(snapshots), "core": identity(core), "preparation": identity(WORK / "preparation.json"),
          "reader": identity(WORK / "run_transfer.py"), "actions": identity(WORK / "actions.json")}, indent=2))


if __name__ == "__main__":
    main()
