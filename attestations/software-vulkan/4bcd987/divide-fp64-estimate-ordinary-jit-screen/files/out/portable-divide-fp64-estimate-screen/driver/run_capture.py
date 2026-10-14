"""One whole-module canary under passive ordinary LLVM/JIT observation."""
from pathlib import Path
import hashlib
import json
import os
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
_, previous_file_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
file_limit = 64 * 1024**2
if previous_file_limit != resource.RLIM_INFINITY:
    file_limit = min(file_limit, previous_file_limit)
resource.setrlimit(resource.RLIMIT_FSIZE, (file_limit, file_limit))
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
folder = Path(sys.argv[1]).resolve()
module = Path(sys.argv[2]).resolve()
stage = root / "bin/linux-gcc/portable-divide-fp64-estimate-screen/driver"
executable = stage / "short_input_jit_probe"
observer = root / "attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so"
assert folder.is_relative_to(work) and folder.exists()
assert module.is_relative_to(work / "module-inputs")
context = json.loads((work / "runtime-context.json").read_text())
assert context["host"]["sha256"] == hashlib.sha256(executable.read_bytes()).hexdigest()
matched = [row for row in context["modules"].values() if Path(row["resolved_path"]) == module]
assert len(matched) == 1 and matched[0]["sha256"] == hashlib.sha256(module.read_bytes()).hexdigest()
assert hashlib.sha256(observer.read_bytes()).hexdigest() == "a1f8808cd8a3a793c59eeed10d78515a08d5aa09b1080ac47281c1a78c04fd68"
symbols = subprocess.check_output(["nm", "--defined-only", str(executable)], text=True).splitlines()
def symbol(part):
    matched = [line.split()[-1] for line in symbols if part in line and ".cold" not in line]
    assert len(matched) == 1, matched
    return matched[0]
dispatch = symbol("VulkanDevice8DispatchE")
load = symbol("VulkanDevice10LoadKernelE")
env = os.environ.copy()
env.update(SIRIUS_GEOMETRY_CAPTURE=str(folder), SIRIUS_DISPATCH_SYMBOL=dispatch,
           SIRIUS_LOAD_SYMBOL=load, SIRIUS_MODULE_INPUT=str(module),
           SIRIUS_FEATURE_LAYOUT=str(work / "feature-layout.json"),
           SIRIUS_FEATURE_LOADER="/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275")
script = folder / "capture.gdb"
assert not script.exists()
script.write_text("set pagination off\nset confirm off\nset breakpoint pending on\n"
    "set debuginfod enabled off\nset print thread-events off\nset disable-randomization off\n"
    "set environment LD_PRELOAD " + str(observer) + "\n"
    "set environment SIRIUS_OPT_IR_FOLDER " + str(folder) + "\n"
    "source " + str(work / "capture_jit.py") + "\n"
    "source " + str(work / "capture_features.py") + "\nrun\n")
inputs = {"module": str(module), "module_sha256": hashlib.sha256(module.read_bytes()).hexdigest(),
          "probe_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
          "observer_sha256": hashlib.sha256(observer.read_bytes()).hexdigest(),
          "dispatch_symbol": dispatch, "load_symbol": load,
          "inherited_core_limit_bytes": 0, "inherited_per_file_limit_bytes": file_limit,
          "scope": "Short-input whole-module compilation/code observation; no retained arithmetic, timing, coldness, or qualification verdict."}
(folder / "capture-inputs.json").write_text(json.dumps(inputs, indent=2) + "\n")
subprocess.run(["/usr/bin/gdb", "-q", "-batch", "-nx", "-x", str(script), "--args",
                str(executable), str(module)], env=env, check=True)
assert not (folder / "capture-error.txt").exists()
assert json.loads((folder / "inferior-exit.json").read_text())["exit_code"] == 0
loaded = json.loads((folder / "actual-loaded-module.json").read_text())
assert loaded["sha256"] == inputs["module_sha256"] and loaded["actual_span_matches_frozen_input"]
bindings = json.loads((folder / "bindings.json").read_text())
assert len(bindings["dispatches"]) == 1 and bindings["dispatches"][0]["completed"]
assert bindings["dispatches"][0]["kernel"] == 0 and bindings["dispatches"][0]["groups_x"] == 1
calls = [item for item in bindings["entry_calls"] if item["dispatch_id"] == 0]
assert len(calls) == 1 and calls[0]["kernel"] == 0
lookups = [item for item in bindings["pointer_lookups"] if calls[0]["file"] in item["matched_files"]]
generations = {item["module_generation"] for item in lookups}
assert len(generations) == 1
generation = generations.pop()
created = bindings["module_creations"][generation]
assert created["id"] == generation
assert all(item["module"] == created["module"] for item in lookups)
exports = [item for item in bindings["ir_exports"] if item["module_generation"] == generation]
assert len(exports) == 1 and exports[0]["optimized_module_observed"]
export = exports[0]
assert export["module"] == created["module"]
raw = (folder / export["path"]).read_bytes()
assert len(raw) == export["bytes"] and hashlib.sha256(raw).hexdigest() == export["sha256"]
passes = [bindings["optimizer_passes"][i] for i in export["successful_optimizer_pass_ids"]]
assert all(item["completed"] and item["error_pointer"] == 0 for item in passes)
assert all(item["id"] == i and item["module_generation"] == generation
           and item["module"] == created["module"]
           for i, item in zip(export["successful_optimizer_pass_ids"], passes))
assert [item["pipeline"] for item in passes] == ["default<O0>",
    "sroa,early-cse,simplifycfg,reassociate,mem2reg,instsimplify,instcombine<no-verify-fixpoint>"]
(folder / "capture-verification.json").write_text(json.dumps({"passed": True,
    "actual_loaded_module": loaded, "actual_entry": calls[0], "module_generation": generation,
    "ordinary_ir_export": export, "ordinary_pipelines": [item["pipeline"] for item in passes],
    "scope": inputs["scope"], "complete_linked_section_and_relocation_review_pending": True}, indent=2) + "\n")

# This local reader does not require the enclosing owner to be terminal.
# Original IR/JIT relocation analysis runs separately after root sees terminal success.
import importlib.util
spec = importlib.util.spec_from_file_location("raw_feature_arguments", work / "read_features.py")
feature_reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(feature_reader)
feature_report = feature_reader.inspect(folder)
(folder / "feature-verification.json").write_text(json.dumps(feature_report, indent=2) + "\n")
