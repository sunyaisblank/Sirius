"""Read/disassemble exact frozen Camera/RayCamera artifacts only; no compiler/device use."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
sys.dont_write_bytecode=True
ROOT=Path.cwd()
BASE=ROOT/"out/retained-qualification/cold-render-diagnostic/shader-ssa/o1-structure-review"
OWN=Path(__file__).resolve().parent/"current125-structure"

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def main():
    OWN.mkdir(exist_ok=True)
    h=load("bindings",BASE/"reproduce.py")
    checks=load("constraints",BASE.parent/"scalar-replacement-review/review.py")
    semantic=load("semantic_surface",BASE/"transport-clean-compact-review.py")
    stats=load("structures",BASE/"six-stage-review.py")
    oldpath=ROOT/"attestations/software-vulkan/6cee359/ray-camera-first-and-repeat/identity.json"
    newpath=ROOT/"attestations/software-vulkan/125a890/ray-camera-first-and-repeat/identity.json"
    old=json.loads(oldpath.read_text());new=json.loads(newpath.read_text())
    assert old["source"]["revision"]=="6cee35906b092dedd520c99c6a76b752f1f9834b"
    assert new["source"]["revision"]=="125a8906639df0f09e461e6b0ee55bd76e366dae"
    head=subprocess.run(["git","rev-parse","HEAD"],check=True,capture_output=True,text=True).stdout.strip()
    assert head==new["source"]["revision"]
    assert not subprocess.run(["git","status","--porcelain"],check=True,capture_output=True,text=True).stdout.strip()
    modules={name:h.binding(ROOT/name) for name in new["artifacts"] if name.endswith(".spv")}
    assert len(modules)==24
    for name,current in modules.items():assert (current["bytes"],current["sha256"])==(new["artifacts"][name]["bytes"],new["artifacts"][name]["sha256"]),name
    native=[name for name in modules if "_portable" not in name]
    assert len(native)==12
    for name in native:assert modules[name]["sha256"]==old["artifacts"][name]["sha256"],name
    dis=Path("/usr/bin/spirv-dis");tool=h.binding(dis)
    report={"scope":"offline exact Camera/RayCamera structural comparison; no compiler/GPU/numerical/timing verdict", "source_revision":head,
        "baseline_identity":h.binding(oldpath),"current_identity":h.binding(newpath),"current24module_bindings":modules,
        "all12native_byteidentical_6cee":True,"disassembler":tool,"script":h.binding(Path(__file__)),"stages":[]}
    surface=lambda text:semantic.normalized_surface(text[:re.search(r"^\s*%\S+\s*=\s*OpFunction ",text,re.M).start()])
    for stage,previous in [("camera",BASE/"six-stage/camera-o1-ssa.spv"),("ray_camera",BASE/"ray-camera-o1-ssa-dead.spv")]:
        name="bin/linux-gcc/src/sirius/backend/retained/retained_"+stage+"_portable.spv"
        wide=name.replace("_portable.spv","_portable_fp64.spv")
        assert h.sha(previous)==old["artifacts"][name]["sha256"]
        assert old["artifacts"][name]["sha256"]==old["artifacts"][wide]["sha256"]
        assert modules[name]["sha256"]==modules[wide]["sha256"]
        assembly=OWN/(stage+"-125a890.spvasm")
        subprocess.run([str(dis),str(ROOT/name),"-o",str(assembly)],check=True,timeout=30)
        oldtext=previous.with_suffix(".spvasm").read_text();newtext=assembly.read_text()
        checks.constraints(oldtext);checks.constraints(newtext)
        assert surface(oldtext)==surface(newtext)
        before=h.structure(previous.with_suffix(".spvasm"));after=h.structure(assembly)
        before["bytes"]=previous.stat().st_size;after["bytes"]=modules[name]["bytes"]
        entry={"stage":stage,"baseline_module":h.binding(previous),"current_module":modules[name],"current_assembly":h.binding(assembly),
            "baseline_structure":before,"current_structure":after,
            "FindUMsb_before":len(re.findall(r"OpExtInst .*? FindUMsb ",oldtext)),"FindUMsb_after":len(re.findall(r"OpExtInst .*? FindUMsb ",newtext)),
            "FindSMsb_before":len(re.findall(r"OpExtInst .*? FindSMsb ",oldtext)),"FindSMsb_after":len(re.findall(r"OpExtInst .*? FindSMsb ",newtext)),
            "native_contracts_and_exact_interfaces_layout_globals":"PASS","both_public_portable_mode_bytes_equal":True}
        report["stages"].append(entry)
    for name,binding in modules.items():assert h.binding(ROOT/name)==binding
    assert h.binding(dis)==tool
    assert subprocess.run(["git","rev-parse","HEAD"],check=True,capture_output=True,text=True).stdout.strip()==head
    assert not subprocess.run(["git","status","--porcelain"],check=True,capture_output=True,text=True).stdout.strip()
    report["bindings_unchanged_after"]=True
    path=OWN/"report.json";assert not path.exists();path.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"report":h.binding(path),"native12_byteidentity":True,"stages":report["stages"]},indent=2))
if __name__=="__main__":main()
