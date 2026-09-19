"""Compile recovered numerical prototypes with their required float controls.

This checks emitted SPIR-V, not device arithmetic or production integration.
All outputs stay under the checkout's ignored out/ directory.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path(__file__).resolve().parent
OUTPUT = ROOT / "out" / "retained-camera-build"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    global OUTPUT
    probes = ("pair_operations_probe", "pair_camera_probe", "full_camera_probe", "program_camera_probe")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="append", choices=probes)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--compiler")
    args = parser.parse_args()
    OUTPUT = args.output
    compiler = args.compiler or shutil.which("slangc") or "/opt/slang/bin/slangc"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    receipt = OUTPUT / "checks.json"
    receipt.unlink(missing_ok=True)
    sources = {path.name: digest(path) for path in sorted(SOURCE.glob("*.slang"))}
    results = []
    for probe in args.probe or probes:
        for compensated in (False, True):
            name = probe + ("-fp32comp" if compensated else "-fp32")
            raw = OUTPUT / (name + ".compiler.spv")
            assembly = OUTPUT / (name + ".spvasm")
            product = OUTPUT / (name + ".spv")
            flags = ["-DSIRIUS_FP32_COMP"] if compensated else []
            if probe == "program_camera_probe":
                flags.append("-O0")
            command = [compiler, str(SOURCE / (probe + ".slang")), "-I", str(SOURCE),
                       *flags, "-target", "spirv", "-profile", "spirv_1_5", "-entry",
                       "ComputeMain", "-stage", "compute", "-denorm-mode-fp32", "preserve",
                       "-o", str(raw)]
            subprocess.run(command, check=True, timeout=600)
            subprocess.run(["spirv-dis", str(raw), "-o", str(assembly)], check=True)
            text = assembly.read_text()
            assert "OpCapability Float64" not in text
            assert "OpTypeFloat 64" not in text and "OpTypeInt 64" not in text
            assert " Fma " not in text
            entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text).group(1)
            assert re.search(r"OpExecutionMode " + re.escape(entry) + r" DenormPreserve 32", text)
            operations = re.findall(r"(%\S+) = OpF(?:Add|Sub|Mul) ", text)
            decorated = set(re.findall(r"OpDecorate (%\S+) NoContraction", text))
            assert operations and all(operation in decorated for operation in operations)
            # The compiler exposes denormal control but not this RTE32 entry
            # mode. Preserve the original prototype's explicit assembly edit.
            text = text.replace("OpCapability Shader", "OpCapability Shader\n"
                                "               OpCapability RoundingModeRTE", 1)
            local = re.search(r"^.*OpExecutionMode " + re.escape(entry) + r" LocalSize.*$",
                              text, re.M).group(0)
            text = text.replace(local, local + "\n               OpExecutionMode " + entry +
                                " RoundingModeRTE 32", 1)
            assembly.write_text(text)
            subprocess.run(["spirv-as", "--target-env", "spv1.5", str(assembly), "-o",
                            str(product)], check=True)
            subprocess.run(["spirv-val", "--target-env", "vulkan1.2", str(product)], check=True)
            results.append({"probe": name, "command": command, "sha256": digest(product),
                            "ordered_arithmetic_instructions": len(operations)})
            print(name + ": compiled and validated", flush=True)
    assert sources == {path.name: digest(path) for path in sorted(SOURCE.glob("*.slang"))}
    receipt.write_text(json.dumps({"scope": "prototype compilation only; no runtime admission",
                                   "sources": sources, "products": results}, indent=2) + "\n")


if __name__ == "__main__":
    main()
