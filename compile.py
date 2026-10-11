"""One bounded compiler/validator diagnostic, with no device invocation."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import struct
import sys
import time

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


def bounded_run(command, **kwargs):
    index = len(records)
    started = time.monotonic()
    stdout_path = WORK / f"tool-{index:02d}.stdout"
    stderr_path = WORK / f"tool-{index:02d}.stderr"
    record = {"command": command, "tool": seal(command[0]), "timeout_seconds": 300,
              "stdout": str(stdout_path), "stderr": str(stderr_path)}
    records.append(record)
    with stdout_path.open("wb") as out, stderr_path.open("wb") as err:
        process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True)
        record.update(pid=process.pid, birth=birth(process.pid), disposition="running")
        persist()
        try:
            code = process.wait(timeout=300)
        except subprocess.TimeoutExpired:
            import signal
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            record.update(disposition="timeout", elapsed_seconds=time.monotonic() - started,
                          returncode=process.returncode)
            persist()
            raise
    record.update(disposition="completed", elapsed_seconds=time.monotonic() - started,
                  returncode=code, stdout_seal=seal(stdout_path), stderr_seal=seal(stderr_path))
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
        result.update({"error_type": type(error).__name__, "error": str(error)})
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        result["source_unchanged"] = seal(source) == before["source"]
        result["stage_inputs_unchanged"] = [seal(path) for path in sorted(source.parent.iterdir())] == stage_inputs
        result["commands"] = len(records)
        (WORK / "compile-result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
