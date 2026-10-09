"""One fixed offline cheap cleanup sequence; Transport first, RayCamera only on success."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OWN = HERE / "cheap-cleanup"
ROOT = HERE.parents[4]
FLAGS = ["--target-env=vulkan1.2", "--eliminate-local-single-block",
         "--eliminate-local-single-store", "--eliminate-dead-code-aggressive",
         "--compact-ids", "--preserve-bindings", "--preserve-interface"]


def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    OWN.mkdir(parents=True,exist_ok=True)
    h=load("bindings",HERE/"reproduce.py")
    checks=load("constraints",HERE.parent/"scalar-replacement-review/review.py")
    stats=load("structures",HERE/"six-stage-review.py")
    semantic=load("semantic_surface",HERE/"transport-clean-compact-review.py")
    identity_path=ROOT/"attestations/software-vulkan/850d9e7/ray-camera-first-and-repeat/identity.json"
    historical=json.loads(identity_path.read_text())
    tool_paths=[Path("/usr/bin/spirv-opt"),Path("/usr/bin/spirv-val"),Path("/usr/bin/spirv-dis")]
    tools_before={str(path):h.binding(path) for path in tool_paths}
    previous=json.loads((HERE/"six-stage/report.json").read_text())
    assert all(binding==previous["tools_before"][name] for name,binding in tools_before.items())
    kernels={name:h.binding(ROOT/name) for name in historical["artifacts"] if name.startswith("src/sirius/kernels/")}
    for name,binding in kernels.items():
        prior=historical["artifacts"][name]
        assert (binding["bytes"],binding["sha256"])==(prior["bytes"],prior["sha256"]),name
    preserved=[HERE/"six-stage/report.json",HERE/"six-stage/transport-2g/report.json",HERE/"six-stage/transport-clean-compact/report.json"]
    source_context=lambda:{"head":subprocess.run(["git","rev-parse","HEAD"],cwd=ROOT,check=True,capture_output=True,text=True).stdout.strip(),
        "status":subprocess.run(["git","status","--porcelain"],cwd=ROOT,check=True,capture_output=True,text=True).stdout}
    report={"scope":"one fixed offline cheap cleanup sequence, Transport then RayCamera only after success; no compilation/GPU/numerical/performance qualification",
        "status":"running","flags":FLAGS,"limits":{"per_tool_wall_seconds":60,"per_tool_address_space_bytes":2147483648},
        "script":h.binding(Path(__file__)),"original_850_identity":h.binding(identity_path),
        "preserved_prior_receipts":[h.binding(path) for path in preserved],"bound_kernel_inputs":kernels,
        "tools_before":tools_before,"source_context_before":source_context(),"commands":[],"stages":[]}

    def limits(): resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))

    def run(name,command):
        out_path=OWN/(name+"-stdout.log")
        err_path=OWN/(name+"-stderr.log")
        entry={"argv":[str(arg) for arg in command]}
        report["commands"].append(entry)
        start=time.monotonic()
        with out_path.open("wb") as out,err_path.open("wb") as err:
            try:
                process=subprocess.Popen(command,stdout=out,stderr=err,preexec_fn=limits)
                entry["pid"]=process.pid
                try: entry["exit_code"]=process.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    process.kill();entry["exit_code"]=process.wait()
                    entry["status"]="timeout-killed-reaped"
                    raise
                if entry["exit_code"]:raise subprocess.CalledProcessError(entry["exit_code"],command)
            finally:
                entry["wall_seconds"]=time.monotonic()-start
                entry["stdout"]=h.binding(out_path);entry["stderr"]=h.binding(err_path)

    def surface(text):
        # Global declarations precede OpFunction under the validated SPIR-V layout.
        # Avoid parsing irrelevant function bodies while comparing exact interfaces.
        first=re.search(r"^\s*%\S+\s*=\s*OpFunction ",text,re.M)
        assert first is not None
        return semantic.normalized_surface(text[:first.start()])

    try:
        for stage,source in [("transport",HERE/"six-stage/transport-o1.spv"),
                             ("ray_camera",HERE/"ray-camera-o1-reproduced.spv")]:
            expected=historical["artifacts"]["bin/linux-gcc/src/sirius/backend/retained/retained_"+stage+"_portable.spv"]
            source_binding=h.binding(source)
            assert (source_binding["bytes"],source_binding["sha256"])==(expected["bytes"],expected["sha256"])
            entry={"stage":stage,"input":source_binding,"input_matches_original_850":True,"status":"running"}
            report["stages"].append(entry)
            output=OWN/(stage+"-o1-cheap.spv")
            assembly=output.with_suffix(".spvasm")
            run(stage+"-optimize",[str(tool_paths[0]),*FLAGS,str(source),"-o",str(output)])
            run(stage+"-validate",[str(tool_paths[1]),"--target-env","vulkan1.2",str(output)])
            run(stage+"-disassemble",[str(tool_paths[2]),str(output),"-o",str(assembly)])
            original=source.with_suffix(".spvasm").read_text();candidate=assembly.read_text()
            checks.constraints(candidate)
            old_surface=surface(original);new_surface=surface(candidate)
            assert old_surface==new_surface,"normalized interface/layout/global declarations differ"
            controls=lambda text:sorted(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)",text))
            assert controls(original)==controls(candidate),"loop/control masks differ"
            entry.update({"status":"passed","output":h.binding(output),"input_structure":stats.structure(source.with_suffix(".spvasm")),
                "output_structure":stats.structure(assembly),"surface_sha256":hashlib.sha256(json.dumps(old_surface,sort_keys=True).encode()).hexdigest(),
                "constraints_normalized_layout_interface_globals_loops_controls":"PASS"})
            (OWN/"report.json").write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps(entry),flush=True)
            del original,candidate,old_surface,new_surface
        report["status"]="passed"
    except Exception as error:
        report["status"]="incomplete";report["error_type"]=type(error).__name__;report["error"]=str(error)
        if report["stages"]:report["stages"][-1]["status"]="incomplete"
    finally:
        assert [h.binding(path) for path in preserved]==report["preserved_prior_receipts"]
        for stage in report["stages"]:assert h.binding(Path(stage["input"]["path"]))==stage["input"]
        assert h.binding(identity_path)==report["original_850_identity"]
        for name,binding in kernels.items():assert h.binding(ROOT/name)==binding,name
        report["tools_after"]={str(path):h.binding(path) for path in tool_paths};assert report["tools_after"]==tools_before
        report["source_context_after"]=source_context()
        report["head_stability_requirement"]="Kernel/module/tool bytes are authoritative; unrelated wrapper-only logging changes may alter HEAD/status."
        (OWN/"report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({key:report.get(key) for key in ["status","error_type","error","stages"]},indent=2))
    return 0 if report["status"]=="passed" else 1


if __name__=="__main__":raise SystemExit(main())
