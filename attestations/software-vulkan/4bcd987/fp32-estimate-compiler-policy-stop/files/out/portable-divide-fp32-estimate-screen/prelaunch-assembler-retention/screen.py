"""One complete FP32-proposal compiler producer; graph acceptance is separate.

Root alone operates the unchanged owner. This source creates no Vulkan objects
and issues no numerical or cost acceptance. Completed outputs receive a producer
seal before any later separately owned structural reader can inspect them.
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
PRODUCTS = ROOT / "bin/linux-gcc/portable-divide-fp32-estimate-screen/compiler"
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
    verify_stock_reuse()
    assert Path.cwd().resolve() == ROOT
    assert not PRODUCTS.exists(), PRODUCTS
    assert not (WORK / "compile-endpoint-normal").exists()
    assert not (WORK / "operator-context.json").exists(), "pin once before consumption"
    assert not any((WORK / name).exists() for name in ("compiler-products.json", "compiler-producer-acceptance.json"))
    assert not subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    assert len(sys.argv) == 3 and revision == sys.argv[2] and re.fullmatch(r"[0-9a-f]{40}", revision)
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


def content(path):
    value = record(path)
    return {key: value[key] for key in ("bytes", "sha256")}


def verify_stock_reuse():
    reuse = json.loads((WORK / "stock-reuse.json").read_text())
    archive = Path(reuse["archive"])
    assert content(archive / "manifest.json")["sha256"] == reuse["manifest_sha256"]
    manifest = json.loads((archive / "manifest.json").read_text())
    for name, expected in reuse["archive_pins"].items():
        assert content(name) == expected, name
        if Path(name) != archive / "manifest.json":
            assert expected == manifest["files"][str(Path(name).relative_to(archive))]
    history = Path(reuse["historical_work"])
    owner = json.loads((history / "compile-endpoint-normal/owner.json").read_text())
    assert owner["status"] == "terminal" and owner["passed"] is True and owner["child_exit"] == 0
    assert owner["stop_reason"] is None and owner["cleanup_errors"] == []
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["remaining_owned_processes"] == [] and owner["observed_births"]
    assert all(row["absent"] for row in owner["observed_births"])
    before_path = history / "compile-endpoint-normal/inputs-before.json"
    assert before_path.read_bytes() == (history / "compile-endpoint-normal/inputs-after.json").read_bytes()
    before = json.loads(before_path.read_text())
    for row in reuse["source_joins"]:
        expected = {key: before[row["original_input_key"]][key] for key in ("bytes", "sha256")}
        assert all(content(row[key]) == expected for key in ("historical", "private", "canonical")), row
    for row in reuse["tool_joins"]:
        expected = before[row["original_input_key"]]
        assert content(row["path"]) == {key: expected[key] for key in ("bytes", "sha256")}
        assert str(Path(row["path"]).resolve(strict=True)) == expected["resolved_path"]
    result = json.loads((history / "compiler-result.json").read_text())
    acceptance = json.loads((history / "compiler-acceptance.json").read_text())
    assert acceptance["passed"] and acceptance["compiler_owner"]["sha256"] == content(history / "compile-endpoint-normal/owner.json")["sha256"]
    assert acceptance["compiler_result"]["sha256"] == content(history / "compiler-result.json")["sha256"]
    assert result["stock_equals_current_canonical"] and result["stock_builder_policy_unchanged"]
    for key, name in (("module", "stock_module"), ("assembly", "stock_assembly")):
        assert content(reuse[name]) == {field: result["modules"]["stock"][key][field] for field in ("bytes", "sha256")}
    assert Path(reuse["stock_module"]).read_bytes() == CANONICAL.read_bytes()
    controls = json.loads(Path(reuse["guard_controls"]).read_text())
    assert controls["passed"] and len(controls["controls"]) == 12
    assert all(row["actual_acceptance"] == row["expected_acceptance"] for row in controls["controls"])
    assert sum(row["expected_acceptance"] for row in controls["controls"]) == 2
    header = (WORK / "candidate/src/sirius/kernels/portable_binary32.h").read_text()
    insertion = (WORK / "candidate-insertion.txt").read_text()
    assert header.count(insertion) == 1
    reversed_header = header.replace(insertion, "", 1)
    revised = "// Exact arithmetic uses uint32, including explicit two-word exact products.\n"
    original = "// All arithmetic uses uint32, including explicit two-word exact products.\n"
    caveat = ("// The same source is used by C++ and Slang. The default path uses no native\n"
              "// floating operation. Optional FP32 quotient proposals are accepted only after\n"
              "// exact integer verification; their accuracy is not a correctness premise.\n")
    assert reversed_header.count(revised) == reversed_header.count(caveat) == 1
    reversed_header = reversed_header.replace(revised, original, 1).replace(caveat,
        "// The same source is used by C++ and Slang; no native floating operation is used.\n", 1)
    assert reversed_header == (WORK / "stock/src/sirius/kernels/portable_binary32.h").read_text()
    assert reversed_header == (ROOT / "src/sirius/kernels/portable_binary32.h").read_text()
    return {"archive_manifest": record(archive / "manifest.json"), "original_owner": record(history / "compile-endpoint-normal/owner.json"),
            "original_acceptance": record(history / "compiler-acceptance.json"), "guard_controls": record(reuse["guard_controls"]),
            "source_joins": len(reuse["source_joins"]), "tool_joins": len(reuse["tool_joins"]),
            "scope": "Exact accepted full stock/compiler-policy reuse; no stock compiler or original Int64-policy execution repeated. New FP32-proposal policy remains separately scoped."}


def normalized_function(text, name):
    body = function_body(text, name)
    definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
    names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))
    local = {identifier: "local" + str(i) for i, identifier in enumerate(re.findall(r"^\s*(%\S+) =", body, re.M))}
    def value(identifier):
        if identifier in local:
            return local[identifier]
        expression = definitions[identifier]
        if expression.startswith("OpFunction "):
            return ["function", names[identifier]]
        assert expression.startswith(("OpType", "OpConstant", "OpSpecConstant")), expression
        return [value(word) if word.startswith("%") else word for word in expression.split()]
    return [[value(word) if word.startswith("%") else word for word in line.split()]
            for line in body.splitlines() if line.strip()]



def compile_candidate():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    _, prior = resource.getrlimit(resource.RLIMIT_FSIZE)
    limit = 64 * 1024**2 if prior == resource.RLIM_INFINITY else min(prior, 64 * 1024**2)
    resource.setrlimit(resource.RLIMIT_FSIZE, (limit, limit))
    verify_frozen()
    reused = verify_stock_reuse()
    context = json.loads((WORK / "operator-context.json").read_text())
    assert context["actions"] == record(WORK / "actions.json")
    assert context["source_freeze"] == record(WORK / "source-freeze.json")
    assert context["preparation"] == record(WORK / "preparation.json")
    assert not PRODUCTS.exists()
    assert not (WORK / "compiler-products.json").exists()
    folder = WORK / "compile-endpoint-normal"
    assert (folder / "inputs-before.json").is_file()
    modules_dir = PRODUCTS / "modules"
    modules_dir.mkdir(parents=True)
    reuse = json.loads((WORK / "stock-reuse.json").read_text())
    stock = Path(reuse["stock_module"])
    stock_assembly = Path(reuse["stock_assembly"])
    recipe = json.loads((WORK / "generated/endpoint-arguments.json").read_text())
    tools = recipe["tools"]
    builder = load_module("fp32_estimate_candidate_builder", WORK / "candidate/scripts/build-retained-kernels.py")
    candidate = modules_dir / "candidate-endpoint-normal.spv"
    builder.compile_shader(WORK / "candidate/src/sirius/kernels/retained_endpoint.slang",
        candidate, tools["compiler"], tools["assembler"], tools["disassembler"], tools["validator"],
        recipe["registers"], 5, recipe["layers"], recipe["prefix"], portable=True,
        optimizer=tools["optimizer"], normal_sum32=True,
        endpoint_specialized=tuple(recipe["endpoint_specialized"]), diagnostic_divide_fp32_estimate=True)
    candidate_assembly = candidate.with_suffix(".final.spvasm")
    subprocess.run([tools["disassembler"], str(candidate), "-o", str(candidate_assembly)], check=True)
    texts = {"stock": stock_assembly.read_text(), "candidate": candidate_assembly.read_text()}
    assert stock.read_bytes() == CANONICAL.read_bytes() and stock.read_bytes() != candidate.read_bytes()
    assert interface_signature(texts["stock"]) == interface_signature(texts["candidate"])
    modes = lambda text: sorted(re.findall(r"OpExecutionMode %\S+ (.*)$", text, re.M))
    assert modes(texts["stock"]) == modes(texts["candidate"])
    def barriers(text):
        definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
        return [[opcode, *[definitions[word] for word in operands.split()]]
                for opcode, operands in re.findall(r"^\s*(OpControlBarrier|OpMemoryBarrier) (.+)$", text, re.M)]
    assert barriers(texts["stock"]) == barriers(texts["candidate"])
    modules = {}
    for mode, output, assembly in (("stock", stock, stock_assembly), ("candidate", candidate, candidate_assembly)):
        text = texts[mode]
        capabilities = re.findall(r"OpCapability (\S+)", text)
        assert sorted(capabilities) == sorted(["Shader", "RoundingModeRTE"])
        widths = re.findall(r"OpTypeInt (\d+) [01]", text)
        assert widths and set(widths) == {"32"}
        assert re.findall(r"OpTypeFloat (\d+)", text) == ["32"]
        assert "OpConvertFToU" not in text and "OpConvertFToS" not in text
        modules[mode] = {"module": record(output), "assembly": record(assembly),
                         "capabilities": capabilities, "integer_widths": ["32"], "float_widths": ["32"]}
    verify_frozen()
    result_path = WORK / "compiler-products.json"
    write("compiler-products.json", {"schema": "sirius-fp32-estimate-compiler-products-v1",
          "passed": True, "original_compile_ssa_assemble_validate_complete": True,
          "stock_reuse": reused, "stock_equals_current_canonical": True, "modules": modules,
          "complete_descriptor_and_workgroup_type_decorations_equal": True,
          "execution_modes_equal": True, "ordered_barrier_operands_equal": True,
          "private_single_FDiv32_policy_and_original_normal_policy_passed": True,
          "structural_graph_accepted": False, "runtime_or_numerical_admission": False,
          "benefit_established": False,
          "scope": "Successful complete normal Endpoint compiler producer and type/interface/mode/barrier policy only. Exact correction/dominator/fallback graph remains for a separate reader against these immutable outputs."})
    write("compile-endpoint-normal/output-seal.json", {
          "schema": "sirius-fp32-estimate-compiler-output-seal-v1", "label": "compile-endpoint-normal",
          "source_freeze": record(WORK / "source-freeze.json"),
          "operator_context": record(WORK / "operator-context.json"),
          "inputs_before": record(folder / "inputs-before.json"),
          "recipe": record(WORK / "compiler-recipes.json"),
          "outputs": {str(path): record(path) for path in (candidate, candidate_assembly, result_path)},
          "structural_graph_accepted": False})
    print("Complete FP32-proposal Endpoint compiler outputs sealed; structural reader and device work unexecuted.", flush=True)


def review_terminal():
    verify_frozen(); verify_stock_reuse()
    context = json.loads((WORK / "operator-context.json").read_text())
    folder = WORK / "compile-endpoint-normal"
    owner = json.loads((folder / "owner.json").read_text())
    action = json.loads((WORK / "actions.json").read_text())["compile-endpoint-normal"]
    assert owner["status"] == "terminal" and owner["passed"] is True and owner["child_exit"] == 0
    assert owner["stop_reason"] is None and owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["source"] == {"revision": context["revision"], "status": ""}
    assert owner["argv"] == action["argv"] and owner["limits"] == action["limits"]
    assert owner["observed_births"] and all(row["absent"] for row in owner["observed_births"])
    before = json.loads((folder / "inputs-before.json").read_text())
    assert (folder / "inputs-before.json").read_bytes() == (folder / "inputs-after.json").read_bytes()
    for name, pin in before.items():
        path = Path(name); path = ROOT / path if not path.is_absolute() else path
        assert content(path) == {key: pin[key] for key in ("bytes", "sha256")}
        assert str(path.resolve(strict=True)) == pin["resolved_path"]
    assert context["actions"] == record(WORK / "actions.json")
    assert context["source_freeze"] == record(WORK / "source-freeze.json")
    assert context["preparation"] == record(WORK / "preparation.json")
    seal = json.loads((folder / "output-seal.json").read_text())
    assert seal["label"] == "compile-endpoint-normal" and seal["structural_graph_accepted"] is False
    for key, path in (("source_freeze", WORK / "source-freeze.json"),
                      ("operator_context", WORK / "operator-context.json"),
                      ("inputs_before", folder / "inputs-before.json"),
                      ("recipe", WORK / "compiler-recipes.json")):
        assert seal[key] == record(path)
    expected = [PRODUCTS / "modules/candidate-endpoint-normal.spv",
                PRODUCTS / "modules/candidate-endpoint-normal.final.spvasm", WORK / "compiler-products.json"]
    assert seal["outputs"] == {str(path): record(path) for path in expected}
    result = json.loads((WORK / "compiler-products.json").read_text())
    assert result["passed"] and result["structural_graph_accepted"] is False
    assert not result["runtime_or_numerical_admission"] and not result["benefit_established"]
    assert not (WORK / "compiler-producer-acceptance.json").exists()
    write("compiler-producer-acceptance.json", {
          "schema": "sirius-fp32-estimate-compiler-producer-acceptance-v1", "passed": True,
          "compiler_owner": record(folder / "owner.json"), "output_seal": record(folder / "output-seal.json"),
          "compiler_products": record(WORK / "compiler-products.json"),
          "structural_graph_accepted": False, "runtime_or_numerical_admission": False,
          "scope": "One successful original tool producer with exact consumed-input and output-seal joins. A separately admitted owned structural reader is still required; no graph/cost/numerical conclusion."})
    print("Terminal compiler producer accepted; separate structural graph admission remains pending.", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "pin-context": pin_context()
    elif sys.argv[1] == "compile": compile_candidate()
    elif sys.argv[1] == "review": review_terminal()
    else: raise SystemExit("unknown compiler-producer phase")
