"""One authorized Transport-only retry; original 1 GiB receipt stays unchanged."""
import difflib
import importlib.util
import json
from pathlib import Path
import re
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OWN = HERE / "six-stage/transport-2g"
ROOT = HERE.parents[4]
FLAGS = ["--target-env=vulkan1.2", "--ssa-rewrite", "--eliminate-dead-code-aggressive",
         "--preserve-bindings", "--preserve-interface"]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    OWN.mkdir(parents=True, exist_ok=True)
    h = load("binding_helpers", HERE / "reproduce.py")
    checks = load("module_surface_checks", HERE.parent / "scalar-replacement-review/review.py")
    stats = load("module_structure", HERE / "six-stage-review.py")
    prior_path = HERE / "six-stage/report.json"
    prior = json.loads(prior_path.read_text())
    source = HERE / "six-stage/transport-o1.spv"
    expected = json.loads((ROOT / "attestations/software-vulkan/850d9e7/ray-camera-first-and-repeat/identity.json").read_text())["artifacts"][
        "bin/linux-gcc/src/sirius/backend/retained/retained_transport_portable.spv"]
    assert h.sha(source) == expected["sha256"] and source.stat().st_size == expected["bytes"]
    tools = [Path("/usr/bin/spirv-opt"), Path("/usr/bin/spirv-val"), Path("/usr/bin/spirv-dis")]
    before = {str(path): h.binding(path) for path in tools}
    for name, identity in before.items():
        assert identity == prior["tools_before"][name]
    report = {"scope": "one offline Transport-only 2 GiB retry; no numerical/GPU/timing verdict",
        "status": "running", "prior_incomplete_receipt": h.binding(prior_path),
        "original_1g_error": prior["error"], "input": h.binding(source), "flags": FLAGS,
        "limits": {"per_tool_wall_seconds": 60, "per_tool_address_space_bytes": 2147483648},
        "tools_before": before, "commands": [], "script": h.binding(Path(__file__))}
    output = OWN / "transport-o1-ssa.spv"
    assembly = output.with_suffix(".spvasm")

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (2147483648, 2147483648))

    def run(name, command):
        stdout = OWN / (name + "-stdout.log")
        stderr = OWN / (name + "-stderr.log")
        entry = {"argv": [str(part) for part in command]}
        report["commands"].append(entry)
        with stdout.open("wb") as out, stderr.open("wb") as err:
            try:
                process = subprocess.run(command, stdout=out, stderr=err,
                    check=True, timeout=60, preexec_fn=limits)
                entry["exit_code"] = process.returncode
            except subprocess.CalledProcessError as error:
                entry["exit_code"] = error.returncode
                raise
            except subprocess.TimeoutExpired:
                entry["status"] = "timeout"
                raise
            finally:
                entry["stdout"] = h.binding(stdout)
                entry["stderr"] = h.binding(stderr)

    try:
        run("optimize", [str(tools[0]), *FLAGS, str(source), "-o", str(output)])
        run("validate", [str(tools[1]), "--target-env", "vulkan1.2", str(output)])
        run("disassemble", [str(tools[2]), str(output), "-o", str(assembly)])
        original = source.with_suffix(".spvasm").read_text()
        candidate = assembly.read_text()
        checks.constraints(candidate)
        assert checks.surface(original) == checks.surface(candidate)
        controls = lambda text: sorted(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)", text))
        assert controls(original) == controls(candidate)
        report["input_structure"] = stats.structure(source.with_suffix(".spvasm"))
        report["output_structure"] = stats.structure(assembly)
        report["output"] = h.binding(output)
        report["constraints_layout_interface_globals_loops_controls"] = "PASS"
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "incomplete"
        report["error"] = str(error)
        raise
    finally:
        assert h.binding(prior_path) == report["prior_incomplete_receipt"]
        assert h.binding(source) == report["input"]
        report["tools_after"] = {str(path): h.binding(path) for path in tools}
        assert report["tools_after"] == before
        report["source_after"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip()
        assert report["source_after"] == prior["source_before"]
        assert not subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip()
        (OWN / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    current = (ROOT / "scripts/build-retained-kernels.py").read_text()
    old = '    subprocess.run([compiler, str(source), *definitions, "-I", str(source.parent), "-O0",\n'
    new = ('    # Inline portable helpers before the existing SSA cleanup; preserve native -O0.\n'
           '    optimization = "-O1" if portable else "-O0"\n'
           '    subprocess.run([compiler, str(source), *definitions, "-I", str(source.parent), optimization,\n')
    assert current.count(old) == 1
    candidate_source = current.replace(old, new).replace(
        "        # Promote local state without the default inlining/unrolling expansion.\n",
        "        # Remove local temporaries before driver compilation.\n")
    patch = ''.join(difflib.unified_diff(current.splitlines(True), candidate_source.splitlines(True),
        fromfile="a/scripts/build-retained-kernels.py", tofile="b/scripts/build-retained-kernels.py"))
    patch_path = OWN / "product-o1-ssa-candidate-v2.patch"
    patch_path.write_text(patch)
    (OWN / "build-retained-kernels-candidate-v2.py").write_text(candidate_source)
    subprocess.run(["git", "apply", "--check", "--whitespace=error", str(patch_path)], cwd=ROOT, check=True)
    print(json.dumps({key: report[key] for key in ["status", "input", "output", "input_structure", "output_structure"]} |
                     {"candidate_patch": h.binding(patch_path)}, indent=2))


if __name__ == "__main__":
    main()
