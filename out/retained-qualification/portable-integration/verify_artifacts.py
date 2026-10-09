"""Fresh offline verification of exact canonical-header probe binaries."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]

def sha(data):
    return hashlib.sha256(data).hexdigest()

report = {"scope":"host/Slang artifact checks only; no GPU execution or retained-route admission",
          "verified_utc":datetime.now(timezone.utc).isoformat(),"modules":{},
          "ubsan_existing_output_equal":(ROOT/"host-ubsan-output.txt").read_bytes()==(ROOT/"expected.txt").read_bytes(),
          "header_sha256":sha((REPO/"src/sirius/kernels/portable_binary32.h").read_bytes()),
          "gcc":subprocess.check_output(["g++-14","--version"],text=True).splitlines()[0],
          "slang":subprocess.check_output(["/opt/slang/bin/slangc","-version"],text=True,stderr=subprocess.STDOUT).strip()}
assert report["ubsan_existing_output_equal"]
for stem in ("probe","probe_o3"):
    binary = (ROOT/(stem+".spv")).read_bytes()
    with tempfile.TemporaryDirectory(prefix="verify-",dir=ROOT) as directory:
        snapshot = Path(directory)/"probe.spv"
        snapshot.write_bytes(binary)
        validation = subprocess.run(["spirv-val","--target-env","vulkan1.2",str(snapshot)],
                                    capture_output=True,text=True,check=True)
        assembly = subprocess.check_output(["spirv-dis",str(snapshot)])
    text = assembly.decode("utf-8")
    capabilities = re.findall(r"OpCapability (\w+)",text)
    widths = re.findall(r"OpTypeInt (\d+) [01]",text)
    modes = re.findall(r"OpExecutionMode .+",text)
    assert capabilities==["Shader"] and widths and set(widths)=={"32"}
    assert "OpTypeFloat" not in text and "OpExecutionModeId" not in text
    assert len(modes)==1 and modes[0].endswith("LocalSize 64 1 1")
    (ROOT/(stem+".spvasm")).write_bytes(assembly)
    report["modules"][stem]={"bytes":len(binary),"sha256":sha(binary),
        "validation_exit":validation.returncode,"fresh_disassembly_sha256":sha(assembly),
        "capabilities":capabilities,"integer_widths":widths,"float_types":0,"execution_modes":modes}
report["source_sha256"]={name:sha((ROOT/name).read_bytes()) for name in
                         ("host_probe.cpp","probe.slang","exact_product_reference.py","verify_artifacts.py")}
report["input_sha256"]=sha((ROOT/"shader-input.bin").read_bytes())
report["expected_sha256"]=sha((ROOT/"shader-expected.bin").read_bytes())
(ROOT/"artifact-report.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
