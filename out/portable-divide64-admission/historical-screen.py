"""Prepare a clean context, compile one complete Endpoint pair, then read its terminal owner.

No mode creates a Vulkan instance, device, pipeline or dispatch. Root operates
the unchanged copied owner; runtime feature enablement and JIT capture are deferred.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import resource
import stat
import struct
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
PRODUCTS = ROOT / "bin/linux-gcc/portable-divide64-screen"
CANONICAL = ROOT / "bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv"
PRODUCT_REVISION = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"


def record(path):
    path = Path(path)
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"bytes": path.stat().st_size, "sha256": digest,
            "mode": stat.S_IMODE(path.stat().st_mode),
            "resolved_path": str(path.resolve(strict=True)),
            "symlink_target": os.readlink(path) if path.is_symlink() else None}


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def verify_frozen():
    frozen = json.loads((WORK / "source-freeze.json").read_text())
    for name, expected in frozen["private_inputs"].items():
        assert record(WORK / name) == expected, name
    for name, expected in frozen["external_inputs"].items():
        assert record(name) == expected, name
    for path in (WORK / "stock/src/sirius/kernels").iterdir():
        assert path.read_bytes() == (ROOT / "src/sirius/kernels" / path.name).read_bytes(), path
    assert (WORK / "stock/scripts/build-retained-kernels.py").read_bytes() == (
        ROOT / "scripts/build-retained-kernels.py").read_bytes()


def pin_context():
    verify_frozen()
    assert Path.cwd().resolve() == ROOT
    assert not PRODUCTS.exists(), PRODUCTS
    assert not (WORK / "compile-endpoint-normal").exists()
    assert not (WORK / "operator-context.json").exists(), "pin once before consumption"
    assert not subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    assert (WORK / "actions.json").read_bytes() == (WORK / "actions.unpinned.json").read_bytes()
    actions = json.loads((WORK / "actions.json").read_text())
    assert list(actions) == ["compile-endpoint-normal"]
    assert actions["compile-endpoint-normal"]["revision"] == "PENDING_OPERATOR_CLEAN_CONTEXT"
    actions["compile-endpoint-normal"]["revision"] = revision
    write("actions.json", actions)
    write("operator-context.json", {"revision": revision, "status": "",
          "product_revision": PRODUCT_REVISION,
          "actions": record(WORK / "actions.json"),
          "source_freeze": record(WORK / "source-freeze.json"),
          "preparation": record(WORK / "preparation.json"),
          "scope": "Clean source pin only; no compiler or device was invoked by pin-context."})
    print("Pinned clean operator context", revision, flush=True)


def function_body(text, name):
    match = re.search(r'OpName (%\S+) "' + re.escape(name) + '"', text)
    assert match, name
    found = re.search(r"^\s*" + re.escape(match[1]) + r" = OpFunction\b.*?OpFunctionEnd",
                      text, re.M | re.S)
    assert found, name
    return found.group(0) + "\n"


def interface_signature(text):
    """Normalize IDs only for the complete descriptor/shared type closure."""
    definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
    names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))
    decorations = {}
    for owner, suffix in re.findall(r"Op(?:Member)?Decorate (%\S+) (.+)$", text, re.M):
        decorations.setdefault(owner, []).append(suffix)

    def signature(identifier, active=()):
        assert identifier not in active, "cyclic interface type"
        line = definitions[identifier]
        assert line.startswith(("OpType", "OpConstant ")), line
        values = [signature(word, active + (identifier,)) if word.startswith("%") else word
                  for word in line.split()]
        return [values, sorted(decorations.get(identifier, []))]

    answer = []
    for identifier, line in definitions.items():
        words = line.split()
        if words[0] == "OpVariable" and words[2] in (
                "StorageBuffer", "Uniform", "Workgroup", "PushConstant"):
            answer.append([names.get(identifier, ""), words[2], signature(words[1]),
                           sorted(decorations.get(identifier, []))])
    assert answer
    return sorted(answer, key=lambda row: json.dumps(row))


def compile_pair():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    _, prior = resource.getrlimit(resource.RLIMIT_FSIZE)
    limit = 64 * 1024**2 if prior == resource.RLIM_INFINITY else min(prior, 64 * 1024**2)
    resource.setrlimit(resource.RLIMIT_FSIZE, (limit, limit))
    verify_frozen()
    context = json.loads((WORK / "operator-context.json").read_text())
    assert context["actions"] == record(WORK / "actions.json")
    assert context["source_freeze"] == record(WORK / "source-freeze.json")
    assert context["preparation"] == record(WORK / "preparation.json")
    assert not PRODUCTS.exists(), PRODUCTS
    (PRODUCTS / "modules").mkdir(parents=True)
    recipe = json.loads((WORK / "generated/endpoint-arguments.json").read_text())
    tools = recipe["tools"]
    texts, modules = {}, {}
    for owner in ("stock", "candidate"):
        builder = load_module("divide64_" + owner + "_builder",
                              WORK / owner / "scripts/build-retained-kernels.py")
        output = PRODUCTS / "modules" / (owner + "-endpoint-normal.spv")
        keywords = {"portable": True, "optimizer": tools["optimizer"],
                    "normal_sum32": True,
                    "endpoint_specialized": tuple(recipe["endpoint_specialized"])}
        if owner == "candidate":
            keywords["diagnostic_divide64"] = True
        builder.compile_shader(WORK / owner / "src/sirius/kernels/retained_endpoint.slang",
            output, tools["compiler"], tools["assembler"], tools["disassembler"], tools["validator"],
            recipe["registers"], 5, recipe["layers"], recipe["prefix"], **keywords)
        assembly = output.with_suffix(".final.spvasm")
        subprocess.run([tools["disassembler"], str(output), "-o", str(assembly)], check=True)
        text = assembly.read_text()
        body = function_body(text, "PB32Divide")
        body_path = output.with_suffix(".PB32Divide.spvasm")
        body_path.write_text(body)
        widths = dict(re.findall(r"(%\S+) = OpTypeInt (\d+) [01]", text))
        div_widths = [int(widths[kind]) for kind in re.findall(r"= OpUDiv (%\S+) ", body)]
        if owner == "stock":
            assert div_widths and set(div_widths) == {32}, div_widths
        else:
            assert div_widths == [64], div_widths
            assert body.count("OpLoopMerge") == 0
            assert re.findall(r"OpTypeInt 64 ([01])", text) == ["0"]
            wide = next(identifier for identifier, width in widths.items() if width == "64")
            assert re.search(r"= OpIMul " + re.escape(wide) + r" ", body)
            assert re.search(r"= OpISub " + re.escape(wide) + r" ", body)
        texts[owner] = text
        modules[owner] = {"module": record(output), "assembly": record(assembly),
                          "divide_body": record(body_path), "divide_udiv_widths": div_widths,
                          "divide_loop_count": body.count("OpLoopMerge"),
                          "capabilities": re.findall(r"OpCapability (\S+)", text),
                          "integer_widths": sorted(set(widths.values()))}
    stock = PRODUCTS / "modules/stock-endpoint-normal.spv"
    candidate = PRODUCTS / "modules/candidate-endpoint-normal.spv"
    assert stock.read_bytes() == CANONICAL.read_bytes(), "stock differs from canonical complete module"
    assert candidate.read_bytes() != stock.read_bytes()
    assert interface_signature(texts["stock"]) == interface_signature(texts["candidate"])
    modes = lambda text: sorted(re.findall(r"OpExecutionMode %\S+ (.*)$", text, re.M))
    assert modes(texts["stock"]) == modes(texts["candidate"])
    for opcode in ("OpControlBarrier", "OpMemoryBarrier"):
        assert texts["stock"].count(opcode) == texts["candidate"].count(opcode)
    expected = {"stock": {"Shader", "RoundingModeRTE"},
                "candidate": {"Shader", "Int64", "RoundingModeRTE"}}
    for owner, result in modules.items():
        assert len(result["capabilities"]) == len(expected[owner])
        assert set(result["capabilities"]) == expected[owner]
    verify_frozen()
    write("compiler-result.json", {"schema": "sirius-divide64-compiler-screen/v1",
          "stock_equals_current_canonical": True, "modules": modules,
          "complete_descriptor_and_workgroup_type_decorations_exact": True,
          "execution_modes_exact": True, "barrier_counts_exact": True,
          "stock_builder_policy_unchanged": True,
          "candidate_only_allowance": "Int64 and unsigned internal 64-bit integers; original float/order/layout/barrier checks retained",
          "runtime_or_numerical_admission": False, "benefit_established": False,
          "next_review": "Inspect complete wide division/remainder lowering and code/frame costs on the actual explicitly enabled provider; compiler output alone does not establish a native divide instruction or benefit."})
    print("Complete Endpoint compiler pair accepted structurally; runtime and benefit remain unmeasured.", flush=True)


def review_terminal():
    verify_frozen()
    folder = WORK / "compile-endpoint-normal"
    owner = json.loads((folder / "owner.json").read_text())
    assert owner["status"] == "terminal" and owner["passed"] is True
    assert owner["child_exit"] == 0 and owner["stop_reason"] is None
    assert owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["observed_births"] and all(row["absent"] for row in owner["observed_births"])
    assert (folder / "inputs-before.json").read_bytes() == (folder / "inputs-after.json").read_bytes()
    result = json.loads((WORK / "compiler-result.json").read_text())
    for row in result["modules"].values():
        for key in ("module", "assembly", "divide_body"):
            assert record(row[key]["resolved_path"]) == row[key]
    write("compiler-acceptance.json", {"compiler_owner": record(folder / "owner.json"),
          "compiler_result": record(WORK / "compiler-result.json"), "passed": True,
          "scope": "Terminal compiler ownership and structural screen only; no Vulkan execution or numerical/performance verdict."})
    print("Terminal compiler owner and unchanged inputs accepted.", flush=True)


if __name__ == "__main__":
    {"pin-context": pin_context, "compile": compile_pair, "review": review_terminal}[sys.argv[1]]()
