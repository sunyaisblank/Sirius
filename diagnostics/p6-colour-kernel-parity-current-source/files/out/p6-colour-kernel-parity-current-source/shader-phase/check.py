"""Post-owner original P6 checks with the accepted inventory/OverlayFS join."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import os
import sys

sys.dont_write_bytecode = True
PHASE = Path(__file__).resolve().parent
WORK = PHASE.parent
ROOT = WORK.parents[1]
CONTEXT = json.loads((PHASE / "context.json").read_text())
ACTIONS = json.loads((PHASE / "actions.json").read_text())
REVIEW = WORK / "post-owner-overlay-review"


def identity(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return {"resolved_path": str(path.resolve(strict=True)),
            "bytes": path.stat().st_size, "sha256": h.hexdigest()}


def accepted_inventory():
    assert identity(REVIEW / "read_inventory.py") == CONTEXT["inventory_reader"]
    spec = importlib.util.spec_from_file_location("frozen_p6_overlay_inventory", REVIEW / "read_inventory.py")
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)  # Its main/result writer is not invoked.
    join = prior.preserved_inputs()
    path = REVIEW / "analysis-inventory.json"
    assert identity(path) == CONTEXT["accepted_inventory_analysis"]
    result = json.loads(path.read_text())
    assert result["passed"] is True and result["case"] == "inventory"
    assert result["phase"] == "separate_post_owner_overlay_analysis"
    assert result["qualification_claimed"] is False
    assert result["original_checker_reported_failure_preserved"] is True
    assert result["original_checker_sha256"] == prior.CHECKER_SHA256
    assert result["original_preservation_manifest_sha256"] == prior.MANIFEST_SHA256
    assert result["readonly_backing_join_sha256"] == prior.JOIN_SHA256
    assert result["owner"] == identity(WORK / "inventory/owner.json")
    assert result["selected_device"] == CONTEXT["accepted_device"]
    assert not (WORK / "accepted-inventory.json").exists()
    return prior, join


def reader():
    prior, join = accepted_inventory()
    # Reuse the exact already-reviewed path/inode/hash/device predicate. This
    # extracts one function definition; it does not call the old main/writer.
    tree = ast.parse((REVIEW / "read_inventory.py").read_text())
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    predicate = [n for n in main.body if isinstance(n, ast.FunctionDef) and n.name == "device_identity_matches"]
    assert len(predicate) == 1
    scope = {"os": os, "THREE": prior.THREE, "join": join,
             "identity": prior.identity, "substitutions": []}
    exec(compile(ast.Module(body=predicate, type_ignores=[]), str(REVIEW / "read_inventory.py"), "exec"), scope)
    path = REVIEW / "original/files/check.py"
    assert identity(path) == CONTEXT["original_checker"]
    source = path.read_text()
    device = '            assert tuple(int(x, 16) for x in parts[3].split(":")) == (os.major(stat.st_dev), os.minor(stat.st_dev))'
    inventory = '        validate("inventory")'
    assert source.count(device) == 1 and source.count(inventory) == 1
    evaluated = source.replace(device, '            assert _device_identity_matches(path, parts, stat)')
    evaluated = evaluated.replace(inventory, '        _accepted_inventory()')
    namespace = {"__name__": "frozen_original_p6_shader_reader", "__file__": str(WORK / "check.py"),
                 "_device_identity_matches": scope["device_identity_matches"],
                 "_accepted_inventory": accepted_inventory}
    exec(compile(evaluated, str(path), "exec"), namespace)
    # Relocate only the two action folders and their frozen preparation inputs.
    namespace.update(WORK=PHASE, ROOT=ROOT, CONTEXT=CONTEXT, ACTIONS=ACTIONS)
    return namespace, scope, hashlib.sha256(evaluated.encode()).hexdigest()


if __name__ == "__main__":
    assert len(sys.argv) == 3 and sys.argv[1] in {"--preflight", "--accept"}
    case = sys.argv[2]
    assert case in {"cie", "xyz"}
    original, scope, evaluated_sha256 = reader()
    if sys.argv[1] == "--preflight":
        original["preflight"](case)
        print("Prepared shader inputs and accepted inventory remain exact; no provider/test launched.")
    else:
        result = original["validate"](case)
        result.update(phase="separate_original_p6_shader_phase", qualification_claimed=False,
                      original_checker=CONTEXT["original_checker"],
                      evaluated_reader_sha256=evaluated_sha256,
                      accepted_inventory_analysis=CONTEXT["accepted_inventory_analysis"],
                      exact_overlay_device_join_substitutions=scope["substitutions"])
        with (PHASE / ("accepted-" + case + ".json")).open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
        print(json.dumps(result, indent=2, sort_keys=True))
