#!/usr/bin/env python3
"""Seal one unlaunched original P4 scene; no product/tool execution."""
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
COMMIT = "7f14771a30cfafc9d12b0f8fee316bd1107ab2be"
BUNDLE = ROOT / "attestations/scientific-core/4bcd987/p4-vulkan-scene-incomplete"
OLD_PREFIX = "out/science-p4-vulkan-session-scenes"
ARCHIVE = BUNDLE / "files" / OLD_PREFIX
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
LIMITS = {"seconds": 1200, "rss_kib": 4194304, "reserve_kib": 2097152}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT).decode().strip()


def key(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def identity(path):
    path = Path(path)
    with path.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"bytes": path.stat().st_size, "sha256": sha,
            "resolved_path": str(path.resolve(strict=True))}


def write(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")


def relocate(name):
    if name == OLD_PREFIX or name.startswith(OLD_PREFIX + "/"):
        return key(ARCHIVE) + name[len(OLD_PREFIX):]
    return name


def relocated(value):
    if isinstance(value, str):
        return relocate(value)
    if isinstance(value, list):
        return [relocated(v) for v in value]
    if isinstance(value, dict):
        return {relocate(k): relocated(v) for k, v in value.items()}
    return value


def render_layout(reader):
    replacements = [
        ('    temporary = WORK / (label + "-temporary")',
         '    temporary = ROOT / "renders/p4-doppler-scene2-outer1200"'),
        ('    if sys.argv[1] == "--verify-only":',
         '    require(temporary == ROOT / prepared["owned_render_output"]["root_relative_to_checkout"], "declared render root changed")\n'
         '    require(temporary.parent.resolve(strict=True) == ROOT / "renders"\n'
         '            and not temporary.exists() and not temporary.is_symlink(), "private render root must be absent under canonical renders")\n'
         '    if sys.argv[1] == "--verify-only":'),
        ('                          "outer_diagnostic_limits": LIMITS}, indent=2), flush=True)',
         '                          "outer_diagnostic_limits": LIMITS, "owned_render_root": str(temporary)}, indent=2), flush=True)'),
        ('        outputs = {str(p.relative_to(WORK)): identity(p) for p in sorted(entries)}',
         '        require(temporary.resolve(strict=True) == temporary\n'
         '                and all(p.resolve(strict=True).parent == temporary for p in entries), "render output escaped owned root")\n'
         '        outputs = {str(p.relative_to(ROOT)): identity(p) for p in sorted(entries)}'),
    ]
    for before, after in replacements:
        require(reader.count(before) == 1, "render layout correction locus")
        reader = reader.replace(before, after)
    return reader


def main():
    require(WORK == ROOT / "out/science-p4-doppler-scene2-outer1200", "fixed new root")
    require(not (WORK / "preparation.json").exists(), "already sealed")
    require(git("rev-parse", "HEAD") == RUNNER and not git("status", "--porcelain"), "clean pinned runner")
    require(git("rev-parse", "evidence/issue-53-p4-vulkan-scene-incomplete^{commit}") == COMMIT, "immutable evidence tag")
    manifest = json.loads((BUNDLE / "manifest.json").read_text())
    require(len(manifest["files"]) == 188, "original manifest coverage")
    raw_files = {key(BUNDLE / name): pin for name, pin in manifest["files"].items()}
    for name in ("manifest.json", "result.json", "preserve-finite-result.py"):
        raw_files[key(BUNDLE / name)] = identity(BUNDLE / name)
    require(not any("/scene-2/" in name for name in raw_files), "scene2 must remain unlaunched")

    # Whole archived bytes must join both their manifest and immutable Git blob.
    names = sorted(raw_files)
    process = subprocess.Popen(["git", "cat-file", "--batch"], cwd=ROOT,
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    joins = {}
    for name in names:
        process.stdin.write((COMMIT + ":" + name + "\n").encode())
        process.stdin.flush()
        header = process.stdout.readline().decode().strip().split()
        require(len(header) == 3 and header[1] == "blob", "missing evidence Git blob:" + name)
        raw = process.stdout.read(int(header[2]))
        require(process.stdout.read(1) == b"\n", "Git batch framing")
        actual = identity(ROOT / name)
        expected = raw_files[name]
        require(all(actual[k] == expected[k] for k in ("bytes", "sha256")), "archived manifest changed:" + name)
        require(raw == (ROOT / name).read_bytes(), "archived bytes differ from immutable Git:" + name)
        joins[name] = {**actual, "git_blob": header[0], "revision": COMMIT}
    process.stdin.close()
    require(process.wait() == 0, "Git object read failed")

    original = json.loads((ARCHIVE / "preparation.json").read_text())
    original_actions = json.loads((ARCHIVE / "actions.json").read_text())
    require(len(original["pins"]) == 474 and len(original["source_recovery"]) == 120, "original seal coverage")
    require(len(original["retained_payload"]["all_embedded_arrays"]) == 42, "full embedded arrays")
    require(len(original["runtime_resource_volume"]) == 32, "full resource volume")
    pins = {}
    relocations = {}
    for name, expected in original["pins"].items():
        target = relocate(name)
        require(not target.startswith("out/"), "no consumed/live scratch input")
        actual = identity(ROOT / target)
        require(all(actual[k] == expected[k] for k in ("bytes", "sha256")), "original input changed:" + name)
        if target == name:
            require(actual == expected, "original resolved input identity changed:" + name)
        else:
            require(target in joins, "relocated input lacks immutable Git join:" + name)
            relocations[name] = target
        pins[target] = actual
    for name, join in joins.items():
        pins[name] = {k: join[k] for k in ("bytes", "sha256", "resolved_path")}

    scene = copy.deepcopy(original["cases"][1])
    name = "VulkanRenderSession.ThinAndVolumetricDopplerSuppressionAffectLiveEmission"
    require(scene["name"] == name and scene["label"] == "scene-2" and len(scene["original_frames"]) == 4, "original scene selection")
    source = (ROOT / "tests/render/vulkan_render_test.cpp").read_text()
    body = source.split("TEST(VulkanRenderSession, " + name.split(".")[1] + ") {", 1)[1].split("\nTEST(", 1)[0]
    require(body.count('const auto thin_') == 2 and body.count('const auto volume_') == 2
            and "RenderSessionVulkan(40, 24" in body and "config.volumetric_samples = 4" in body
            and body.count("1.0e-5") == 2, "unchanged four-frame source assertions")
    registration = (ROOT / "bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake").read_text()
    require(scene["original_registration"] in registration, "original registration changed")
    properties = next(line for line in registration.splitlines() if line.startswith("set_tests_properties(") and name in line)
    require("TIMEOUT" not in properties and "sirius_vulkan_device" in properties, "original resource lock/no whole-scene timeout")
    require("inline constexpr int kRenderTraceMaximumAttempts = 20000;" in (ROOT / "src/sirius/render/trace_domain.h").read_text(), "original ray attempt bound")
    require("inline constexpr double kDispatchStopMs = 1000.0;" in (ROOT / "src/sirius/render/dispatch_governor.h").read_text(), "original feedback bound")

    for filename, digest in (("run_owned.py", "492c75609a0dc8ff5eaa0299f7206820bb8a42f85c3ecc5b9b0af98f2496ba9f"),
                             ("owned_guard.py", "747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056")):
        require(identity(ARCHIVE / filename)["sha256"] == digest, "mature owner/guard identity")
        shutil.copyfile(ARCHIVE / filename, WORK / filename)

    reader = (ARCHIVE / "run_scene.py").read_text()
    replacements = [
        ('"seconds": 300', '"seconds": 1200'),
        ('require(Path.cwd().resolve() == ROOT and len(sys.argv) == 2 and sys.argv[1] in ("1", "2"), "run from checkout with ordinal1 or2")',
         'require(Path.cwd().resolve() == ROOT and len(sys.argv) == 2 and sys.argv[1] in ("2", "--verify-only"), "run from checkout with original ordinal2 or --verify-only")'),
        ('index = int(sys.argv[1]) - 1', 'index = 0\n    require(len(prepared["cases"]) == 1 and prepared["cases"][0]["label"] == "scene-2", "only original scene2 selected")'),
    ]
    for before, after in replacements:
        require(reader.count(before) == 1, "narrow reader replacement absent")
        reader = reader.replace(before, after)
    start = reader.index('    for prior in prepared["cases"][:index]:')
    end = reader.index('    result = {"passed": False', start)
    reader = reader[:start] + reader[end:]
    # Validation has no output/TMPDIR creation and never imports the owner.
    start = reader.index('    try:\n        require(git("rev-parse", "HEAD")')
    end = reader.index('        temporary.mkdir(mode=0o700)', start)
    validations = reader[start + len('    try:\n'):end]
    validation_top = '\n'.join(line[4:] if line.startswith('    ') else line for line in validations.split('\n'))
    reader = reader[:start] + '    try:\n' + reader[end:]
    insert_at = reader.index('    require(not any(p.exists()')
    reader = reader[:insert_at] + validation_top + '''    if sys.argv[1] == "--verify-only":
        print(json.dumps({"verified_only": True, "runtime_actions": 0,
                          "case": case["name"], "pins": len(prepared["pins"]),
                          "preparation": identity(WORK / "preparation.json"),
                          "outer_diagnostic_limits": LIMITS}, indent=2), flush=True)
        return 0
''' + reader[insert_at:]
    before = 'finite_scope=prepared["finite_scope"], no_performance_or_cache_warmness_claim=True)'
    after = 'finite_scope=prepared["finite_scope"], no_performance_or_cache_warmness_claim=True,\n                      whole_child_wall_seconds=owner["wall_seconds"],\n                      peak_sampled_rss_kib=owner["peak_sampled_rss_kib"],\n                      minimum_sampled_available_kib=owner["minimum_sampled_available_kib"])'
    require(reader.count(before) == 1, "result extension locus")
    reader = reader.replace(before, after)
    before = '        result["error"] = repr(error)\n        save()'
    after = '''        result["error"] = repr(error)
        if (folder / "owner.json").is_file():
            observed_owner = json.loads((folder / "owner.json").read_text())
            result.update(owner=identity(folder / "owner.json"),
                          whole_child_wall_seconds=observed_owner.get("wall_seconds"),
                          peak_sampled_rss_kib=observed_owner.get("peak_sampled_rss_kib"),
                          minimum_sampled_available_kib=observed_owner.get("minimum_sampled_available_kib"),
                          owner_stop_reason=observed_owner.get("stop_reason"),
                          disposition="INCOMPLETE_DIAGNOSTIC_LIMIT" if observed_owner.get("stop_reason") is not None else "NOT_ACCEPTED",
                          complete_original_numerical_verdict=False)
        save()'''
    require(reader.count(before) == 1, "honest incomplete result locus")
    reader = reader.replace(before, after)
    previous_layout = reader
    reader = render_layout(reader)
    (WORK / "recovery").mkdir(exist_ok=True)
    (WORK / "recovery/render-output-path.patch").write_text("".join(difflib.unified_diff(
        previous_layout.splitlines(True), reader.splitlines(True),
        fromfile="released-diagnostic-tmp/run_scene.py", tofile="renders-owned/run_scene.py")))
    ast.parse(reader)
    (WORK / "run_scene.py").write_text(reader)

    context = relocated({k: v for k, v in original.items() if k != "pins"})
    context["cases"] = [scene]
    context["finite_scope"] = "Only original P4 scene2: four complete40x24/1SPP thin/volume times Doppler frames, original two finite same-size nonempty linear-radiance difference assertions>1e-5. No fullP4/science/Mandatory/release/domain/performance/coldKerr300 claim."
    context["original_internal_controls"]["outer"] = "1200s finite outer diagnostic allowance for this first four-frame scene2 attempt; unchanged4GiB sampled process-tree RSS/2GiB sampled MemAvailable reserve. No whole-scene source timeout or completion/performance guarantee; original coldKerr300 criterion unchanged."
    context["original_internal_controls"]["ray_attempts"] = {"maximum": 20000, "source": "src/sirius/render/trace_domain.h", "unchanged": True}
    context["new_packet_scheduling"] = {"selected_original_ordinals": [2], "scene2_previously_unlaunched": True,
         "scene1_repeated": False, "preceding_scene_acceptance_required": False,
         "reason": "Original packet preceding-case check was private scheduling, not scene2 numerical/product precondition. Root schedules this one first attempt only after the parity matrix retires."}
    context["outer_allowance_basis"] = "Four complete original retained frames with setup/work/teardown included. 1200s is a declared finite operational allowance, not measured completion prediction or an original performance threshold. No scene/settings/sample/attempt/assertion change."
    context["frame_fixture"] = {"width": 40, "height": 24, "samples_per_pixel": 1, "tile_size": 64,
         "metric": "Kerr", "mass": 1.0, "spin": 0.9, "observer_distance": 30.0,
         "observer_inclination_degrees": 80, "camera_fov_degrees": 50,
         "volumetric_samples_if_enabled": 4, "other_SessionConfig_defaults_unchanged": True,
         "numerical_authority": "Original returned GetDisplayBuffer().SnapshotFloatData vectors; separate thin and volume pairs nonempty/equal length/finite and sum absolute represented-float differences>1e-5. EXR header checks only transport metadata; no independent EXR pixel-decode or original tape identity claim."}
    context["owned_render_output"] = {"root_relative_to_checkout": "renders/p4-doppler-scene2-outer1200",
         "private_TMPDIR": True, "required_mode_octal": "0700", "must_be_absent_before_launch": True,
         "output_identity_paths_relative_to": "checkout_root", "required_output_names": scene["outputs"],
         "containment": "Canonical renders parent, no root symlink, regular non-symlink outputs resolving directly inside the owned root.",
         "diagnostics_remain_under": key(WORK), "original_frames_and_assertions_changed": False,
         "reversible_reader_diff": key(WORK / "recovery/render-output-path.patch")}
    context["archive_reuse"] = {"commit": COMMIT, "tag": "evidence/issue-53-p4-vulkan-scene-incomplete",
         "bundle": key(BUNDLE), "original_input_count": len(original["pins"]),
         "relocated_deleted_scratch_inputs": len(relocations), "manifest_files": len(manifest["files"]),
         "whole_git_joined_archive_files": len(joins), "bindings": key(WORK / "recovery/archive-bindings.json"),
         "owner_and_guard_byte_identical": True, "source_recovery_count": len(original["source_recovery"]),
         "generated_recovery_count": 40, "full_embedded_array_count": 42, "full_resource_count": 32}
    (WORK / "recovery").mkdir(exist_ok=True)
    write(WORK / "recovery/archive-bindings.json", {"archive_revision": COMMIT, "files": joins, "original_to_archived_paths": relocations})
    write(WORK / "context.json", context)
    # No restoration of deleted scratch or dependency on another active root.
    for row in context["source_recovery"]:
        require(row["recovery"].startswith(key(ARCHIVE) + "/recovery/source/"), "immutable source recovery")
    frozen = sorted(set(pins) | {key(p) for p in WORK.rglob("*") if p.is_file()}
                    | {key(WORK / filename) for filename in ("preparation.json", "actions.json")})
    action = copy.deepcopy(original_actions["scene-2"])
    action.update(limits=LIMITS, scope=context["finite_scope"], frozen_inputs=frozen)
    require(action["argv"] == original_actions["scene-2"]["argv"] and action["software_provider"] is True,
            "original registered runtime action")
    write(WORK / "actions.json", {"scene-2": action})
    pins.update({key(p): identity(p) for p in WORK.rglob("*") if p.is_file()})
    context["pins"] = dict(sorted(pins.items()))
    write(WORK / "preparation.json", context)
    require(set(action["frozen_inputs"]) == set(pins) | {key(WORK / "preparation.json")}, "owner complete closure")
    require(not any(name.startswith("out/") and not name.startswith(key(WORK) + "/") for name in pins), "no other live scratch dependencies")
    print(json.dumps({"prepared_only": True, "runtime_actions": 0, "pins": len(pins),
                      "preparation": identity(WORK / "preparation.json"),
                      "reader": identity(WORK / "run_scene.py"), "actions": identity(WORK / "actions.json"),
                      "whole_git_joined_archive_files": len(joins)}, indent=2))


if __name__ == "__main__":
    main()
