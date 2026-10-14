"""Prepare five unlaunched original cases from immutable evidence; no tool/test execution."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil
import subprocess

ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
COMMIT = "0d8beec5fdb1406d35add385cfa6f2d9f62848d2"
BASE = ROOT / "attestations/scientific-core/4bcd987/p4-ordinary-fp32-isolated-partial"
OLD_NAME = "out/science-p4-ordinary-fp32-parity-isolated"
OLD = BASE / "files" / OLD_NAME
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
NAMES = (
    "KernelParity.OpticallyThinGreyLayerAbsorptionMatchesHostAuthority",
    "KernelParity.ArbitraryLatitudeKerrZamoMatchesHostAuthority",
    "KernelParity.EquatorialKerrDiskTransferMatchesHostAuthority",
    "KernelParity.BlackbodyMatchesIntegratedCoreSpectrum",
    "KernelParity.DiskEmissionAppliesExactlyOneGFourthFactor",
)
LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}


def identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(data)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(),
            "resolved_path": str(path.resolve(strict=True))}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def write(path, document):
    path.write_text(json.dumps(document, sort_keys=True, indent=2) + "\n")


def main():
    assert WORK == ROOT / "out/science-p4-ordinary-fp32-parity-cases5-9"
    assert not (WORK / "preparation.json").exists()
    assert git("rev-parse", "HEAD").decode().strip() == RUNNER
    assert not git("status", "--porcelain").strip()
    assert git("rev-parse", "evidence/issue-53-p4-ordinary-fp32-isolated-partial^{commit}").decode().strip() == COMMIT
    manifest_path = BASE / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    entries = {}
    for row in git("ls-tree", "-rz", COMMIT + ":" + key(BASE)).split(b"\0"):
        if row:
            metadata, name = row.split(b"\t", 1)
            mode, kind, blob = metadata.split()
            assert kind == b"blob" and mode in (b"100644", b"100755")
            entries[name.decode()] = blob.decode()
    archived = {str(path.relative_to(BASE)): path for path in BASE.rglob("*") if path.is_file()}
    assert set(archived) == set(entries)
    for name, path in archived.items():
        assert not path.is_symlink()
        raw = path.read_bytes()
        assert hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == entries[name]
    assert len(manifest["files"]) == 291
    for name, expected in manifest["files"].items():
        actual = identity(BASE / name)
        assert all(actual[field] == expected[field] for field in ("bytes", "sha256"))
    original = json.loads((OLD / "preparation.json").read_text())
    original_actions = json.loads((OLD / "actions.json").read_text())
    cases = original["cases"][4:9]
    assert tuple(case["name"] for case in cases) == NAMES
    assert [case["label"] for case in cases] == ["parity-" + str(i) for i in range(5, 10)]
    assert [case["opcode"] for case in cases] == [17, 15, 16, 3, 5]
    assert [case["records"] for case in cases] == [5, 4, 2, 5, 2]
    for ordinal, case in enumerate(cases, 5):
        assert case["minimum_allocated_bytes"] == case["records"] * 320 + 16
        assert case["maximum_allocated_bytes"] == 1048576
        prefix = "files/" + OLD_NAME + "/parity-" + str(ordinal)
        assert not any(name.startswith(prefix + "/") or name == prefix + "-result.json" for name in archived)
    # The failed third owner remains historical incomplete evidence. Its private
    # batch ordering guard is not a dependency of this standalone original test.
    prior = json.loads((OLD / "parity-3/owner.json").read_text())
    assert prior["status"] == "terminal" and prior["stop_reason"] == "diagnostic_timeout"
    assert prior["source_unchanged"] and prior["whole_inputs_unchanged"]
    assert not prior["cleanup_errors"] and not prior["remaining_owned_processes"]
    assert all(row["absent"] for row in prior["observed_births"])
    recovery = WORK / "recovery"
    recovery.mkdir()
    assets = set(archived.values())
    relocations = []
    for name, expected in original["pins"].items():
        path = BASE / "files" / name if name.startswith(OLD_NAME + "/") else ROOT / name
        actual = identity(path)
        assert all(actual[field] == expected[field] for field in ("bytes", "sha256")), name
        if not name.startswith(OLD_NAME + "/"):
            assert actual == expected, name
        assets.add(path)
        relocations.append({"original": name, "current": key(path), "identity": actual})
    assert len(relocations) == 418
    write(recovery / "historical-input-relocations.json", relocations)
    for row in original["source_recovery"]:
        live, recovered = ROOT / row["path"], BASE / "files" / row["recovery"]
        raw = git("show", PRODUCT + ":" + row["path"])
        assert live.read_bytes() == recovered.read_bytes() == raw
        assert hashlib.sha256(raw).hexdigest() == row["sha256"]
    source = (ROOT / "tests/backend/kernel_parity_test.cpp").read_text()
    registration = ROOT / "bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake"
    lines = registration.read_text().splitlines()
    registration_properties = {}
    for case in cases:
        body = source.split("TEST(KernelParity, " + case["name"].split(".")[1] + ") {", 1)[1].split("\nTEST(", 1)[0]
        assert "Fixture f = OpenProbe();" in body and 'OpenProbe("' not in body
        assert case["original_registration"] in lines
        properties = next(line for line in lines if line.startswith("set_tests_properties([=[" + case["name"] + "]=]"))
        assert "RESOURCE_LOCK sirius_vulkan_device" in properties and "SKIP_REGULAR_EXPRESSION" in properties
        registration_properties[case["name"]] = properties
    compiler = next(row for row in json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
                    if row["file"].endswith("/tests/backend/kernel_parity_test.cpp"))
    assert compiler == original["compile_command"]
    assert original["shader_compile_rule"] in (ROOT / "bin/linux-gcc/build.ninja").read_text()
    backend = ROOT / "bin/linux-gcc/tests/backend/sirius_backend_tests"
    shader = backend.parent / "resources/kernels/parity_probe.spv"
    assert identity(backend)["sha256"] == "e5f3b0f120ec6a12490e3e9b6e491d19eb12bfad62c1cdce11722c091a152149"
    assert identity(backend)["bytes"] == 13440952
    assert identity(shader)["sha256"] == "19991852c2bf8025aaa30b4cb3029e7bb980425d74de515efc4879df95eb0891"
    assert identity(shader)["bytes"] == 5562236
    producer_owner = json.loads((ROOT / original["producer"]["build_owner"]).read_text())
    producer_result = json.loads((ROOT / original["producer"]["accepted"]).read_text())
    assert producer_owner["status"] == "terminal" and producer_owner["passed"] and producer_owner["child_exit"] == 0
    assert producer_owner["source_unchanged"] and producer_owner["whole_inputs_unchanged"] and producer_result["passed"]
    for name in ("run_owned.py", "owned_guard.py"):
        shutil.copy2(OLD / name, WORK / name)
    reader = (OLD / "run_parity.py").read_text()
    (recovery / "original-run_parity.py").write_text(reader)
    selected = reader.replace('require(len(sys.argv) == 2 and sys.argv[1] in [str(i) for i in range(1, 10)], "select one ordinal 1..9")',
                              'require(len(sys.argv) == 2 and sys.argv[1] in [str(i) for i in range(5, 10)], "select one original ordinal 5..9")', 1)
    selected = selected.replace('    index = int(sys.argv[1]) - 1\n', '    index = int(sys.argv[1]) - 5\n', 1)
    selected = selected.replace('        temporary.mkdir(mode=0o700)\n',
                                '        require("SIRIUS_PRECISION" not in os.environ, "original ordinary-FP32 precision must remain unset")\n        temporary.mkdir(mode=0o700)\n', 1)
    assert selected != reader
    (WORK / "run_parity.py").write_text(selected)
    (recovery / "reader-selection.patch").write_text("".join(difflib.unified_diff(reader.splitlines(True), selected.splitlines(True),
        fromfile="archived/run_parity.py", tofile="case4/run_parity.py")))
    context = {key: original[key] for key in ("product_source_revision", "runner_source_revision", "diagnostic_only_difference",
               "compile_command", "shader_compile_rule", "spirv_modes", "required_provider_mappings", "provider_policy",
               "producer", "allocation_contract")}
    context.update(schema_version=1, prepared_only=True, cases=cases, original_case_ordinals=list(range(5, 10)),
        finite_scope="Five previously unlaunched original ordinary-FP32 P4 cases5..9; no complete P4/science/Mandatory/platform/release claim.",
        limits=LIMITS, backend=identity(backend), stock_runtime_shader=identity(shader),
        historical_evidence={"commit": COMMIT, "tag": "evidence/issue-53-p4-ordinary-fp32-isolated-partial",
            "bundle": key(BASE), "manifest": identity(manifest_path), "whole_manifest_files": 291,
            "whole_raw_files_Git_joined": len(archived), "historical_inputs_byte_joined": 418,
            "source_recovery_count": len(original["source_recovery"]), "fresh_build_claimed": False},
        scheduling={"cases5_through9_previously_unlaunched": True,
            "case3_remains_incomplete": key(OLD / "parity-3/owner.json"),
            "case3_stop_reason": prior["stop_reason"], "preceding_acceptance_fabricated": False,
            "case4_pass_before5": "Root scheduling decision outside this packet; no live case4 scratch/result dependency",
            "change": "Only copied private ordinal indexing changes to selected5..9; original previous-selected acceptance and raw-result checks remain. Stop on first incomplete/failure.",
            "no_completion_or_performance_forecast": True, "no_rejected_loop_policy_selected": True},
        finite_references={case["name"]: {"reference": case["reference"], "source_line": case["source_line"],
            "records": case["records"], "actual_dispatch_ms_maximum": 1000, "actual_submissions": 1,
            "requested_minimum_bytes": case["minimum_allocated_bytes"], "allocation_maximum_bytes": 1048576}
            for case in cases},
        original_registration_properties=registration_properties, source_recovery=[{**row, "recovery": key(BASE / "files" / row["recovery"])}
            for row in original["source_recovery"]], original_assertions_changed=False, qualification_claimed=False)
    write(WORK / "context.json", context)
    assets.update(path for path in WORK.rglob("*") if path.is_file())
    actions = {}
    frozen = sorted({key(path) for path in assets} | {key(WORK / "preparation.json")})
    for case in cases:
        action = dict(original_actions[case["label"]])
        assert action["limits"] == LIMITS and action["software_provider"] is True
        assert action["argv"][3:5] == ["-R", "^" + re.escape(case["name"]) + "$"]
        action["scope"] = context["finite_scope"] + " Case:" + case["name"]
        action["frozen_inputs"] = frozen
        actions[case["label"]] = action
    write(WORK / "actions.json", actions)
    assets.add(WORK / "actions.json")
    context["pins"] = {key(path): identity(path) for path in sorted(assets)}
    write(WORK / "preparation.json", context)
    print(json.dumps({"prepared_only": True, "cases": NAMES, "pins": len(context["pins"]),
                      "preparation": identity(WORK / "preparation.json")}))


if __name__ == "__main__":
    main()
