#!/usr/bin/env python3
"""Two root-operated one-shot actions; inspect only after owner is terminal."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import runpy
import subprocess
import sys

from hint_transform import derive

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
PHASES = ("validate", "run")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    path = Path(path)
    return {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "resolved_path": str(path.resolve(strict=True))}


def query(*arguments):
    return subprocess.check_output(["git", *arguments], cwd=ROOT, text=True, timeout=30).strip()


def main():
    require(Path.cwd().resolve() == ROOT and len(sys.argv) == 2 and sys.argv[1] in PHASES,
            "select validate or run from the checkout")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"), "owner assertions must be enabled")
    phase = sys.argv[1]
    prepared = json.loads((WORK / "preparation.json").read_text())
    folder, result_path = WORK / phase, WORK / (phase + "-result.json")
    require(not folder.exists() and not result_path.exists(), "one-shot action already attempted")
    result = {"passed": False, "phase": phase, "scope": prepared["scope"],
              "numerical_acceptance": False, "performance_or_adoption_claimed": False,
              "qualification_claimed": False}
    def save():
        result_path.write_text(json.dumps(result, indent=2) + "\n")
    save()
    reader = None
    try:
        require(query("rev-parse", "HEAD") == prepared["runner_source_revision"]
                and not query("status", "--porcelain"), "pinned clean runner changed")
        require(set(query("diff", "--name-only", prepared["product_source_revision"], "HEAD").splitlines())
                == set(prepared["diagnostic_only_difference"]), "product source difference changed")
        for name, expected in prepared["pins"].items():
            require(identity(ROOT / name) == expected, "frozen input changed: " + name)
        candidate, proof = derive((ROOT / prepared["stock_shader"]).read_bytes())
        require(proof == json.loads((WORK / "unroll-bit-diff.json").read_text())
                and candidate == (ROOT / prepared["runtime_shader"]).read_bytes(),
                "exact reversible hint-only transformation changed")
        spec = importlib.util.spec_from_file_location("original_setup_reader", ROOT / prepared["baseline_reader"])
        reader = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reader)
        if phase == "run":
            previous = json.loads((WORK / "validate-result.json").read_text())
            require(previous["passed"] is True, "owned validator has not passed")
            reader.accepted_owner(WORK / "validate/owner.json")
            require(identity(WORK / "validate/owner.json") == previous["owner"], "validator owner changed")
            require("SIRIUS_PRECISION" not in os.environ, "preserve original unset precision")
        temporary = WORK / (phase + "-temporary")
        require(not temporary.exists(), "TMPDIR is not fresh")
        temporary.mkdir(mode=0o700)
        require(temporary.stat().st_uid == os.getuid() and temporary.stat().st_mode & 0o777 == 0o700,
                "TMPDIR ownership/mode differs")
        os.environ["TMPDIR"] = str(temporary)
        previous_argv = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), phase]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = previous_argv
        owner = reader.accepted_owner(folder / "owner.json")
        if phase == "run":
            result["setup_observation"] = reader.review_run(folder, owner, prepared)
        result.update(passed=True, owner=identity(folder / "owner.json"),
                      runtime_shader=identity(ROOT / prepared["runtime_shader"]))
        save()
        print(json.dumps({k: result[k] for k in ("passed", "phase", "scope", "numerical_acceptance", "qualification_claimed")}), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        if phase == "run" and reader is not None and (folder / "stderr.log").is_file():
            try:
                result["partial_setup_observation"] = reader.read_markers(folder)
            except BaseException as nested:
                result["marker_read_error"] = repr(nested)
        save()
        raise


if __name__ == "__main__":
    raise SystemExit(main())
