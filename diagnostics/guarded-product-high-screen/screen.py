"""Bounded host-neighbor and complete-SPIR-V screen; no device/timing verdict."""
import argparse
import collections
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import resource
import shutil
import struct
import subprocess
import sys
import time
import types

sys.dont_write_bytecode = True
REPO = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
PRODUCTS = REPO / "bin/linux-gcc/guarded-product-high-screen"


def record(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def write(name, data):
    (WORK / name).write_text(json.dumps(data, indent=2) + "\n")


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def execute(argv, **kwargs):
    index = len(list((WORK / "actions").glob("*.json")))
    prefix = WORK / "actions" / f"{index:02}"
    started = time.monotonic()
    with prefix.with_suffix(".stdout").open("wb") as stdout, prefix.with_suffix(".stderr").open("wb") as stderr:
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
        child = subprocess.Popen(argv, cwd=REPO, stdout=stdout, stderr=stderr,
                                 preexec_fn=limits)
        try:
            birth = Path(f"/proc/{child.pid}/stat").read_text().rsplit(")", 1)[1].split()[19]
        except FileNotFoundError:
            birth = None
        timed_out = False
        try:
            code = child.wait(timeout=300)
        except subprocess.TimeoutExpired:
            timed_out = True
            child.kill()
            code = child.wait()
    owner = {"argv": list(map(str, argv)), "cwd": str(REPO), "pid": child.pid,
             "saved_start_ticks": birth, "elapsed_seconds": time.monotonic() - started,
             "return_code": code, "timed_out": timed_out,
             "stdout": record(prefix.with_suffix(".stdout")),
             "stderr": record(prefix.with_suffix(".stderr"))}
    prefix.with_suffix(".json").write_text(json.dumps(owner, indent=2) + "\n")
    print(f"action {index}: {Path(argv[0]).name}, exit {code}, {owner['elapsed_seconds']:.3f}s", flush=True)
    if code != 0:
        raise RuntimeError(prefix.with_suffix(".json"))
    return subprocess.CompletedProcess(argv, code)


def prepare():
    assert not (WORK / "stock").exists() and not PRODUCTS.exists()
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
    assert not subprocess.check_output(["git", "status", "--porcelain=v1"], text=True).strip()
    (WORK / "actions").mkdir()
    PRODUCTS.mkdir()
    original = REPO / "src/sirius/kernels"
    for owner in ("stock", "candidate"):
        target = WORK / owner / "src/sirius/kernels"
        target.mkdir(parents=True)
        for path in original.iterdir():
            if path.is_file() and path.suffix in (".h", ".slang", ".py"):
                shutil.copyfile(path, target / path.name)
    header = WORK / "candidate/src/sirius/kernels/portable_binary32.h"
    text = header.read_text()
    needle = "PB32Product PB32ProductResidual(uint a, uint b) {"
    assert text.count(needle) == 1
    text = text.replace(needle, "PB32Product PB32ProductResidualFallback(uint a, uint b) {", 1)
    text = text.replace("#undef SIRIUS_PB32_INLINE", (WORK / "guarded_product_body.inc").read_text() + "\n#undef SIRIUS_PB32_INLINE", 1)
    header.write_text(text)
    pins = json.loads((REPO / "attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/protected-before.json").read_text())
    for path, expected in pins.items():
        assert record(REPO / path) == expected
    write("protected-before.json", pins)
    write("source-freeze.json", {str(p.relative_to(WORK)): record(p)
          for p in sorted(WORK.rglob("*")) if p.is_file() and "actions" not in p.parts})


def host():
    reference = module("pb32_reference", REPO / "tests/support/portable_binary32/reference.py")
    reference.selfcheck()
    counts = collections.Counter()
    rows = []
    corpus = list(reference.residual_cases())
    # Product below an exact binade with an upward hardware estimate.
    for sign_a in (0, 0x80000000):
        for sign_b in (0, 0x80000000):
            corpus.append((0x3f7ff400 | sign_a, 0x3f800600 | sign_b, "downward-binade-witness"))
    for a, b, category in corpus:
        high, low, valid = reference.product_expected(a, b)
        am, bm = a & reference.MASK, b & reference.MASK
        guarded = (0x00800000 <= am <= reference.BOUND and
                   0x00800000 <= bm <= reference.BOUND and
                   (am >> 23) + (bm >> 23) - 254 >= -126)
        if guarded:
            exact = abs(reference.exact(a) * reference.exact(b))
            magnitude = high & reference.MASK
            rounded = reference.exact(magnitude)
            neighbors = ([magnitude] if rounded == exact else
                         [magnitude, magnitude + 1] if rounded < exact else
                         [magnitude - 1, magnitude])
            assert all(0x00800000 <= value < reference.INF for value in neighbors)
        else:
            neighbors = [0xffffffff]  # Poison: fallback must not read the estimate.
        for value in neighbors:
            rows.append((a, b, value, high, low, valid, int(guarded)))
            counts[("guarded/" if guarded else "fallback/") + category] += 1
    fixture = WORK / "neighbor-cases.bin"
    with fixture.open("wb") as output:
        for row in rows:
            output.write(struct.pack("<7I", *row))
    execute(["/usr/bin/g++", "-std=c++20", "-O2", "-Wall", "-Wextra", "-Werror", "-MMD",
             "-MF", str(PRODUCTS / "neighbor_probe.d"), "-I", str(WORK / "candidate/src/sirius/kernels"),
             str(WORK / "neighbor_probe.cpp"), "-o", str(PRODUCTS / "neighbor_probe")])
    execute([str(PRODUCTS / "neighbor_probe"), str(fixture)])
    observed = json.loads((WORK / "actions/01.stdout").read_text())
    assert observed["records"] == len(rows) and observed["corrections"] > 0
    assert observed["downward_binades"] >= 4
    write("host-result.json", {"passed": True, "scope": "Exact C++ correction under both permitted neighboring estimates; no shader/provider execution or all-domain proof.",
          "original_residual_cases": len(corpus) - 4, "additional_binade_cases": 4,
          "categories": dict(counts), "fixture": record(fixture), "observed": observed,
          "reference": record(REPO / "tests/support/portable_binary32/reference.py")})


def compile_modules():
    builder = module("retained_builder", REPO / "scripts/build-retained-kernels.py")
    builder.subprocess = types.SimpleNamespace(run=execute)
    programs = module("retained_program", WORK / "stock/src/sirius/kernels/retained_program.py")
    cache = (REPO / "bin/linux-gcc/CMakeCache.txt").read_text()
    def tool(name):
        return re.search(r"^" + name + r":FILEPATH=(.*)$", cache, re.M)[1]
    compiler, assembler, disassembler, validator, optimizer = map(tool,
        ("SIRIUS_SLANGC", "SIRIUS_SPIRV_AS", "SIRIUS_SPIRV_DIS", "SIRIUS_SPIRV_VAL", "SIRIUS_SPIRV_OPT"))
    write("tools.json", {p: record(Path(p)) for p in (compiler, assembler, disassembler, validator, optimizer, "/usr/bin/g++")})
    modules = WORK / "modules"
    modules.mkdir()
    results = []
    for stage, normal in (("endpoint", False), ("endpoint", True), ("transport", True), ("dopri_phase", True), ("dense", False)):
        program = getattr(programs, "build_" + stage + "_program")(parallel=True)
        kwargs = {}
        if stage in ("endpoint", "transport"):
            special = getattr(programs, "build_schwarzschild_" + stage + "_program")(parallel=True)
            encode = lambda p: [p["instructions"], p["registers"], len(p["outputs"])] + p["outputs"] + p["operations"] + p["layer_offsets"]
            dimensions = (len(encode(program)), len(encode(special)), len(special["layer_offsets"]) - 1)
            if stage == "endpoint": dimensions += (special["prefix_instructions"],)
            kwargs[stage + "_specialized"] = dimensions
        variant = stage + ("-normal" if normal else "-core")
        texts = {}
        for owner in ("stock", "candidate"):
            source = WORK / owner / "src/sirius/kernels" / ("retained_" + stage + ".slang")
            destination = modules / (owner + "-" + variant + ".spv")
            try:
                builder.compile_shader(source, destination, compiler, assembler, disassembler, validator,
                    program["registers"], 5, len(program["layer_offsets"]) - 1,
                    program.get("prefix_instructions", 0), portable=True, optimizer=optimizer,
                    normal_sum32=normal, **kwargs)
                assert owner == "stock"
            except ValueError as error:
                assert owner == "candidate" and str(error) in (
                    "portable retained stage depends on native floating arithmetic",
                    "portable normal sum changed its ordered three-operation body")
                (modules / (owner + "-" + variant + "-production-policy-rejection.txt")).write_text(str(error) + "\n")
                raw = destination.with_suffix(".compiler.spv")
                assembly = destination.with_suffix(".spvasm")
                text = assembly.read_text()
                floating = re.findall(r"(%\S+) = (OpF(?!unction)\S+) ", text)
                multiplications = [(value, op) for value, op in floating if op == "OpFMul"]
                assert len(multiplications) == 1
                value = multiplications[0][0]
                assert "OpDecorate " + value + " NoContraction" in text
                bodies = re.findall(r"%\S+ = OpFunction .*?OpFunctionEnd", text, re.S)
                product_bodies = [body for body in bodies if "= OpFMul " in body]
                assert len(product_bodies) == 1 and re.match(r"%PB32ProductResidual(?:_\d+)? = OpFunction ", product_bodies[0])
                assert re.findall(r"OpCapability (\S+)", text) == ["Shader"]
                assert set(re.findall(r"OpTypeInt (\d+) [01]", text)) == {"32"}
                assert re.findall(r"OpTypeFloat (\d+)", text) == ["32"]
                assert not re.search(r"\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFmaKHR)\b", text)
                without_multiply = re.sub(r"(?m)^.* = OpFMul .*\n", "", text)
                if normal:
                    # Reuse original structural checks for every old floating body.
                    # This filtered text is never assembled or executed.
                    builder.portable_normal_sum_controls(without_multiply)
                    entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
                    text = text.replace("OpCapability Shader", "OpCapability Shader\n               OpCapability RoundingModeRTE", 1)
                    local = re.search(r"^.*OpExecutionMode " + re.escape(entry) + r" LocalSize.*$", text, re.M)[0]
                    text = text.replace(local, local + "\n               OpExecutionMode " + entry + " RoundingModeRTE 32", 1)
                else:
                    assert len(floating) == 1
                assembly.write_text(text)
                execute([assembler, "--target-env", "spv1.5", str(assembly), "-o", str(destination)])
                execute([validator, "--target-env", "vulkan1.2", str(destination)])
            assembly = modules / (owner + "-" + variant + "-final.spvasm")
            execute([disassembler, str(destination), "-o", str(assembly)])
            texts[owner] = assembly.read_text()
        canonical = REPO / "bin/linux-gcc/src/sirius/backend/retained" / ("retained_" + stage + ("_portable_normal_sum.spv" if normal else "_portable.spv"))
        assert (modules / ("stock-" + variant + ".spv")).read_bytes() == canonical.read_bytes()
        for opcode in ("OpControlBarrier", "OpMemoryBarrier"):
            assert len(re.findall(r"\b" + opcode + r"\b", texts["stock"])) == len(re.findall(r"\b" + opcode + r"\b", texts["candidate"]))
        layouts = lambda text: re.findall(r"OpExecutionMode %\S+ LocalSize (\d+ \d+ \d+)", text)
        assert layouts(texts["stock"]) == layouts(texts["candidate"])
        results.append({"variant": variant, "stock": record(modules / ("stock-" + variant + ".spv")),
            "candidate": record(modules / ("candidate-" + variant + ".spv")),
            "baseline_byte_exact_to_current_product": True,
            "local_size": layouts(texts["stock"]),
            "instruction_counts": {owner: dict(collections.Counter(re.findall(r"\b(Op\w+)\b", text))) for owner, text in texts.items()}})
    write("compiler-result.json", {"passed": True, "scope": "Private complete SPIR-V screen, original production policy rejects added FMul; no adopted policy, JIT footprint, shader execution or timing verdict.", "modules": results})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("prepare", "host", "compile"))
    phase = parser.parse_args().phase
    if phase != "prepare":
        for name, expected in json.loads((WORK / "source-freeze.json").read_text()).items():
            assert record(WORK / name) == expected, name
    {"prepare": prepare, "host": host, "compile": compile_modules}[phase]()
