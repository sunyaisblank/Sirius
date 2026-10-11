"""One bounded compiler/validator diagnostic, with no device invocation."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import struct
import sys
import time
import traceback

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
records = []


def seal(path):
    path = Path(path)
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def birth(pid):
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
        return text.rsplit(")", 1)[1].split()[19]
    except FileNotFoundError:
        return None


def persist():
    (WORK / "commands.json").write_text(json.dumps(records, indent=2) + "\n")


def group_members(group):
    members, unreadable = [], 0
    for directory in Path("/proc").iterdir():
        if not directory.name.isdecimal():
            continue
        try:
            fields = (directory / "stat").read_text().rsplit(")", 1)[1].split()
            if int(fields[2]) == group:
                members.append({"pid": int(directory.name), "birth": fields[19], "state": fields[0]})
        except (OSError, IndexError, ValueError):
            unreadable += 1
    return {"scope": "point-in-time accessible Linux /proc members of recorded new-session group; excludes escaped groups, hidden/scheduled work",
            "group": group, "members": members, "unreadable_or_raced": unreadable}


def signal_group(group, value):
    try:
        os.killpg(group, value)
    except ProcessLookupError:
        pass


def close_group(process):
    # The encompassing caller finally reaches here after every successful spawn.
    errors, actions = [], []

    def attempt(name, operation):
        try:
            return operation()
        except BaseException as error:
            errors.append({"operation": name, "type": type(error).__name__, "error": str(error)})
            return None

    def scan():
        return attempt("accessible group scan", lambda: group_members(process.pid))

    def send_group(value):
        actions.append({"group": process.pid, "signal": int(value)})
        attempt("signal recorded group", lambda: signal_group(process.pid, value))

    current = scan()
    running = attempt("poll known leader", process.poll) is None
    forced = running or current is None or bool(current["members"])
    if forced:
        send_group(signal.SIGTERM)
    # Reaping remains independent of scan and signal errors.
    attempt("wait known leader after TERM", lambda: process.wait(timeout=5))
    current = scan()
    if process.returncode is None or current is None or current["members"]:
        forced = True
        send_group(signal.SIGKILL)
        if process.returncode is None:
            attempt("kill known leader independently", process.kill)
        attempt("wait known leader after KILL", lambda: process.wait(timeout=5))
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            current = scan()
            if current is None or not current["members"]:
                break
            time.sleep(0.05)
    return {"forced_cleanup": forced, "signal_actions": actions,
            "errors": errors, "terminal_group_scan": scan(),
            "known_leader_reaped": process.returncode is not None}


def bounded_run(command, **kwargs):
    index = len(records)
    started = time.monotonic()
    stdout_path = WORK / f"tool-{index:02d}.stdout"
    stderr_path = WORK / f"tool-{index:02d}.stderr"
    record = {"command": command, "tool": seal(command[0]), "timeout_seconds": 300,
              "stdout": str(stdout_path), "stderr": str(stderr_path)}
    records.append(record)
    process = None
    code = None
    try:
        with stdout_path.open("wb") as out, stderr_path.open("wb") as err:
            process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True)
            record.update(pid=process.pid, birth=birth(process.pid), disposition="running")
            persist()
            code = process.wait(timeout=300)
        record["disposition"] = "completed"
    except BaseException as error:
        record.update(disposition="timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        try:
            if process is not None:
                cleanup = close_group(process)
                record["cleanup"] = cleanup
                record["returncode"] = process.returncode
                terminal = cleanup["terminal_group_scan"]
                if cleanup["errors"] or not cleanup["known_leader_reaped"] or terminal is None or terminal["members"]:
                    raise RuntimeError("recorded tool cleanup incomplete: " + json.dumps(cleanup))
                if cleanup["forced_cleanup"] and record["disposition"] == "completed":
                    record["disposition"] = "forced-cleanup"
                    raise RuntimeError("normal tool completion needed forced cleanup")
        finally:
            record["elapsed_seconds"] = time.monotonic() - started
            for name, path in (("stdout", stdout_path), ("stderr", stderr_path)):
                if path.exists():
                    record[name + "_seal"] = seal(path)
            persist()
    if code:
        raise subprocess.CalledProcessError(code, command)
    return subprocess.CompletedProcess(command, code)


def main():
    cache = (ROOT / "bin/linux-gcc/CMakeCache.txt").read_text()
    tools = {name: re.search(r"^SIRIUS_" + name + r":FILEPATH=(.*)$", cache, re.M)[1]
             for name in ("SLANGC", "SPIRV_AS", "SPIRV_DIS", "SPIRV_VAL", "SPIRV_OPT")}
    spec = importlib.util.spec_from_file_location("literal_screen_builder", ROOT / "scripts/build-retained-kernels.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    programs = json.loads((WORK / "programs.json").read_text())
    general, special = programs["General"], programs["Schwarzschild"]
    source = WORK / "kernels/retained_transport.slang"
    destination = WORK / "literal-transport.spv"
    stage_inputs = [seal(path) for path in sorted(source.parent.iterdir())]
    for path in sorted(source.parent.iterdir()):
        if path != source:
            assert path.read_bytes() == (ROOT / "src/sirius/kernels" / path.name).read_bytes()
    baseline_source = ROOT / "bin/linux-gcc/src/sirius/backend/retained/retained_transport_portable_normal_sum.spv"
    baseline = WORK / "baseline-transport.spv"
    baseline.write_bytes(baseline_source.read_bytes())
    header = ROOT / "bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h"
    table = re.search(r"kTransportPortableNormalSumShader\{\{(.*?)\}\};", header.read_text(), re.S)
    expected = [int(value, 16) if value.startswith("0x") else int(value)
                for value in re.findall(r"(0x[0-9a-fA-F]+|\d+)u?", table[1])]
    assert baseline.read_bytes() == struct.pack("<" + str(len(expected)) + "I", *expected)
    before = {"revision": json.loads((WORK / "emission.json").read_text())["revision"],
              "owner": {"pid": os.getpid(), "birth": birth(os.getpid())},
              "source": seal(source), "emitter": seal(WORK / "emit.py"),
              "runner": seal(Path(__file__)),
              "program_metadata": seal(WORK / "programs.json"),
              "emission_metadata": seal(WORK / "emission.json"),
              "schedule": seal(WORK / "schedule.json"), "header": seal(header),
              "stage_inputs": stage_inputs, "unchanged_imports_equal_current_source": True,
              "builder": seal(ROOT / "scripts/build-retained-kernels.py"),
              "baseline": seal(baseline_source), "baseline_copy": seal(baseline),
              "baseline_equals_current_embedded_array": True,
              "tools": {name: seal(path) for name, path in tools.items()},
              "scope": "only compilation, SPIR-V control inspection, assembly, validation, disassembly; no device workload"}
    (WORK / "compile-before.json").write_text(json.dumps(before, indent=2) + "\n")
    subprocess.run = bounded_run
    started = time.monotonic()
    result = {"scope": before["scope"], "pass": False}
    try:
        bounded_run([tools["SPIRV_DIS"], str(baseline), "-o", str(WORK / "baseline-transport.spvasm")])
        encoded = builder.compile_shader(source, destination, tools["SLANGC"], tools["SPIRV_AS"],
            tools["SPIRV_DIS"], tools["SPIRV_VAL"], general["registers"], 5,
            len(general["layer_offsets"]) - 1, portable=True, optimizer=tools["SPIRV_OPT"],
            normal_sum32=True, transport_specialized=(13096, 2871, len(special["layer_offsets"]) - 1))
        assert len(encoded) * 4 == destination.stat().st_size
        bounded_run([tools["SPIRV_DIS"], str(destination), "-o", str(WORK / "literal-transport.spvasm")])
        result.update({"pass": True, "module": seal(destination),
                       "assembly": seal(WORK / "literal-transport.spvasm")})
    except Exception as error:
        result.update({"pass": False, "error_type": type(error).__name__, "error": str(error),
                       "traceback": traceback.format_exc()})
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        result["source_unchanged"] = seal(source) == before["source"]
        result["stage_inputs_unchanged"] = [seal(path) for path in sorted(source.parent.iterdir())] == stage_inputs
        expected = [before[name] for name in ("source", "emitter", "runner", "program_metadata", "emission_metadata", "schedule", "header", "builder", "baseline", "baseline_copy")]
        expected += stage_inputs + list(before["tools"].values())
        changed = [record["path"] for record in expected if seal(record["path"]) != record]
        result["input_checks"] = {"reference_count": len(expected), "changed": changed}
        if changed or not result["source_unchanged"] or not result["stage_inputs_unchanged"]:
            result["pass"] = False
        result["commands"] = len(records)
        (WORK / "compile-result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result), flush=True)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
