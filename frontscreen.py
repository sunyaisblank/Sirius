"""Isolate the compiler front end; this does not compile a usable shader."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
PHASE = WORK / "frontend"


def seal(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    PHASE.mkdir(exist_ok=False)
    parent = WORK / "compile.py"
    text = parent.read_text()
    assert hashlib.sha256(parent.read_bytes()).hexdigest() == "ad8bc38b743bc3a0cea61b28b126dc58f220278320231cc5a6870fa95ed831f0"
    assert text.count('"timeout_seconds": 300') == text.count('process.wait(timeout=300)') == 1
    # Same reviewed lifecycle wrapper, different diagnostic duration and depth.
    text = text.replace('"timeout_seconds": 300', '"timeout_seconds": 90').replace(
        'process.wait(timeout=300)', 'process.wait(timeout=90)').replace(
        'ROOT = WORK.parents[1]', 'ROOT = WORK.parents[2]')
    runner = PHASE / "phase_runner.py"
    runner.write_text(text)
    spec = importlib.util.spec_from_file_location("literal_frontend_owner", runner)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    commands = json.loads((WORK / "commands.json").read_text())
    assert len(commands) >= 2
    original = commands[1]["command"]
    assert Path(original[0]).name == "slangc" and original[-2] == "-o"
    command = [*original[:-1], str(PHASE / "unused.spv"), "-skip-codegen", "-report-detailed-perf-benchmark"]
    source = Path(command[1])
    refs = [source, parent, runner, Path(__file__), Path(command[0]),
            *sorted(path for path in source.parent.iterdir() if path != source)]
    before = [seal(path) for path in refs]
    facts = {"scope": "90-second complete-source compiler front-end phase only, skip-codegen; no usable shader, validation, device, numerical or benefit result",
             "owner": {"pid": __import__('os').getpid(), "birth": module.birth(__import__('os').getpid())},
             "command": command, "input_seals": before}
    (PHASE / "before.json").write_text(json.dumps(facts, indent=2) + "\n")
    started = time.monotonic()
    result = {"scope": facts["scope"], "frontend_completed": False}
    try:
        module.bounded_run(command)
        result["frontend_completed"] = True
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic() - started,
                  input_seals_unchanged=[seal(path) for path in refs] == before,
                  output_exists=(PHASE / "unused.spv").exists())
    if not result["input_seals_unchanged"]:
        result["frontend_completed"] = False
    (PHASE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
    return 0 if result["frontend_completed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
