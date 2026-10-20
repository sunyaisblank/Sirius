from pathlib import Path
import hashlib, json, re, subprocess, sys
sys.dont_write_bytecode = True
if sys.flags.optimize:
    raise RuntimeError("assertions must remain enabled")
WORK = Path(__file__).resolve().parent

def native_controls(text, fma=False):
    if "OpCapability Float64" in text or "OpTypeFloat 64" in text or "OpTypeInt 64" in text or " Fma " in text:
        raise ValueError("retained stage introduced wide arithmetic or contraction")
    if re.search(r"\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision)\b", text):
        raise ValueError("retained stage introduced relaxed arithmetic")
    operations = re.findall(r"(%\S+) = OpF(?:Add|Sub|Mul) ", text)
    fmas = re.findall(r"(%\S+) = OpFmaKHR ", text)
    decorated = set(re.findall(r"OpDecorate (%\S+) NoContraction", text))
    if (not operations and not fmas) or any(operation not in decorated for operation in operations + fmas):
        raise ValueError("retained arithmetic lost NoContraction")
    if bool(fmas) != fma:
        raise ValueError("retained FMA instructions disagree with the selected variant")
    entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
    if f"OpExecutionMode {entry} DenormPreserve 32" not in text:
        raise ValueError("retained arithmetic lost subnormal preservation")
    capabilities = "OpCapability Shader\n               OpCapability RoundingModeRTE"
    if fma:
        capabilities += "\n               OpCapability FMAKHR\n               OpCapability SignedZeroInfNanPreserve"
    text = re.sub(r"(?m)^([ \t]*)OpCapability Shader$",
                  lambda match: match[1] + capabilities, text, count=1)
    if fma:
        end = list(re.finditer(r"^.*OpCapability .*?$", text, re.M))[-1].end()
        text = text[:end] + '\n               OpExtension "SPV_KHR_fma"' + text[end:]
    local_size = re.search(r"^.*OpExecutionMode " + re.escape(entry) + r" LocalSize.*$", text, re.M)[0]
    modes = "\n               OpExecutionMode " + entry + " RoundingModeRTE 32"
    if fma:
        modes += "\n               OpExecutionMode " + entry + " SignedZeroInfNanPreserve 32"
    return text.replace(local_size, local_size + modes, 1)

def identity(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    assert len(sys.argv) == 2
    output = Path(sys.argv[1]).resolve(strict=True)
    assert output == WORK / "candidate" and output.is_dir()
    recipe = json.loads((WORK / "recipes.json").read_text())
    assert recipe["schema"] == "sirius-native-coefficient-sharing-one-compiler-recipes-v1"
    assert recipe["optimizer_invoked"] is False
    names = ["retained_transport"]
    assert [row["name"] for row in recipe["profiles"]] == names
    source = Path(recipe["source_directory"])
    assert source == WORK / "source/src/sirius/kernels"
    tools = recipe["tools"]
    records = []
    products = {}
    assert not (output / "command-results.json").exists()
    assert not (output / "products.json").exists()
    for name in names:
        for suffix in recipe["expected_products_per_profile"]:
            assert not (output / (name + suffix)).exists()

    def invoke(label, argv):
        stdout = output / (label + ".stdout.log")
        stderr = output / (label + ".stderr.log")
        with stdout.open("xb") as out, stderr.open("xb") as err:
            result = subprocess.run(argv, stdout=out, stderr=err, check=False)
        records.append({"label": label, "argv": argv, "exit_code": result.returncode,
                        "stdout": identity(stdout), "stderr": identity(stderr)})
        (output / "command-results.json").write_text(json.dumps(records, indent=2) + "\n")
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, argv)

    for row in recipe["profiles"]:
        name = row["name"]
        raw = output / (name + ".compiler.spv")
        raw_asm = output / (name + ".raw.spvasm")
        assembly = output / (name + ".spvasm")
        final = output / (name + ".spv")
        final_asm = output / (name + ".final.spvasm")
        invoke(name + "-compile", [tools["compiler"], str(source / (name + ".slang")),
            *row["definitions"], "-I", str(source), "-O0", "-target", "spirv", "-profile", "spirv_1_5",
            "-entry", "ComputeMain", "-stage", "compute", "-denorm-mode-fp32", "preserve", "-o", str(raw)])
        invoke(name + "-raw-disassemble", [tools["disassembler"], str(raw), "-o", str(raw_asm)])
        assembly.write_text(native_controls(raw_asm.read_text(), fma=False))
        invoke(name + "-assemble", [tools["assembler"], "--target-env", "spv1.5", str(assembly), "-o", str(final)])
        invoke(name + "-validate", [tools["validator"], "--target-env", "vulkan1.2", str(final)])
        invoke(name + "-final-disassemble", [tools["disassembler"], str(final), "-o", str(final_asm)])
        products[name] = {suffix: identity(output / (name + suffix))
                          for suffix in recipe["expected_products_per_profile"]}
    assert len(records) == 5 and all(row["exit_code"] == 0 for row in records)
    result = {"schema": "sirius-native-coefficient-sharing-one-compiler-products-v1", "passed": True,
              "scope": "One ordinary native FP32 Transport module compiled/annotated/assembled/validated only",
              "runtime_admitted": False, "emitted_semantics_admitted": False,
              "optimizer_invoked": False, "actual_tool_calls": len(records),
              "recipe": identity(WORK / "recipes.json"), "products": products}
    (output / "products.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
