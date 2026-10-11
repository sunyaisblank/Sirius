"""Inspect the exact validated diagnostic, without loading a Vulkan device."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

WORK = Path(__file__).resolve().parent


def read(path):
    text = path.read_text()
    functions = {match[1]: match[0] for match in re.finditer(
        r"^\s*(%\S+) = OpFunction .*?OpFunctionEnd", text, re.M | re.S)}
    constants = {match[1]: int(match[2]) for match in re.finditer(
        r"^\s*(%\S+) = OpConstant %(?:uint|int) (-?\d+)\s*$", text, re.M)}
    names = {match[1]: match[2] for match in re.finditer(r'OpName (%\S+) "([^"]+)"', text)}
    return text, functions, constants, names


def seal(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def facts(path, text, functions, names):
    return {"module": seal(path), "assembly_bytes": len(text.encode()),
            "functions": len(functions),
            "opcodes": dict(Counter(re.findall(r"\b(Op\w+)", text))),
            "capabilities": re.findall(r"OpCapability (\S+)", text),
            "execution_modes": re.findall(r"^\s*OpExecutionMode .*?$", text, re.M),
            "workgroup_variables": [line.strip() for line in text.splitlines()
                                    if re.search(r"= OpVariable .* Workgroup\b", line)],
            "private_variables": [line.strip() for line in text.splitlines()
                                  if re.search(r"= OpVariable .* Private\b", line)]}


def body_by_name(functions, names, name):
    matches = [(key, body) for key, body in functions.items() if names.get(key) == name]
    assert len(matches) == 1, (name, len(matches))
    return matches[0]


def storage_facts(text, constants):
    types = {m[1]: (m[2], m[3].split()) for m in re.finditer(
        r"^\s*(%\S+) = (OpType\S+) ([^\n]+)", text, re.M)}

    def size(key):
        operation, args = types[key]
        if operation in ("OpTypeInt", "OpTypeFloat"):
            return int(args[0]) // 8
        if operation in ("OpTypeVector", "OpTypeMatrix"):
            return size(args[0]) * int(args[1])
        if operation == "OpTypeArray":
            return size(args[0]) * constants[args[1]]
        if operation == "OpTypeStruct":
            return sum(size(arg) for arg in args)
        raise ValueError(key)

    variables = []
    for match in re.finditer(r"^\s*(%\S+) = OpVariable (%\S+) (Workgroup|Private|Function)(?: ([^\n]+))?", text, re.M):
        operation, args = types[match[2]]
        assert operation == "OpTypePointer" and args[0] == match[3]
        try:
            byte_count = size(args[1])
        except (ValueError, KeyError):
            byte_count = None
        variables.append({"id": match[1], "storage": match[3], "type": args[1],
                          "logical_member_bytes": byte_count, "initializer": match[4]})
    return {"scope": "SPIR-V declared logical scalar/aggregate member bytes; excludes ABI padding, driver lowering and physical residency",
            "variables": variables,
            "declared_workgroup_member_bytes": sum(v["logical_member_bytes"] for v in variables if v["storage"] == "Workgroup"),
            "aggregate_variables_at_least_10000_bytes": [v for v in variables if (v["logical_member_bytes"] or 0) >= 10000]}


def main():
    compile_result = json.loads((WORK / "compile-result.json").read_text())
    compile_before = json.loads((WORK / "compile-before.json").read_text())
    assert compile_result["pass"] and compile_result["source_unchanged"] and compile_result["stage_inputs_unchanged"]
    assert seal(WORK / "literal-transport.spv") == compile_result["module"]
    assert seal(WORK / "literal-transport.spvasm") == compile_result["assembly"]
    assert seal(WORK / "baseline-transport.spv") == compile_before["baseline_copy"]
    commands = json.loads((WORK / "commands.json").read_text())
    assert len(commands) == 7 and all(command["disposition"] == "completed" and command["returncode"] == 0 for command in commands)
    assert commands[0]["command"][1:] == [str(WORK / "baseline-transport.spv"), "-o", str(WORK / "baseline-transport.spvasm")]
    assert commands[-1]["command"][1:] == [str(WORK / "literal-transport.spv"), "-o", str(WORK / "literal-transport.spvasm")]
    baseline, bf, bc, bn = read(WORK / "baseline-transport.spvasm")
    candidate, cf, cc, cn = read(WORK / "literal-transport.spvasm")
    report = {"scope": "structural compiler-only facts; no numerical or timing result",
              "compile_result": seal(WORK / "compile-result.json"),
              "exact_compile_result_module_and_assembly_binding": True,
              "baseline": facts(WORK / "baseline-transport.spv", baseline, bf, bn),
              "diagnostic": facts(WORK / "literal-transport.spv", candidate, cf, cn),
              "literal_functions": {}}
    layer_names = sorted({name for name in cn.values() if re.fullmatch(r"Literal(?:General|Schwarzschild)Layer\d{3}", name)})
    assert len(layer_names) == 159
    report["literal_layer_function_count"] = len(layer_names)
    for name in [*layer_names, "LiteralGeneralProgram", "LiteralSchwarzschildProgram"]:
        key, body = body_by_name(cf, cn, name)
        accesses = re.findall(r"(%\S+) = Op(?:InBounds)?AccessChain %\S+ %retainedScratch ([^\n]+)", body)
        assert accesses, name
        assert all(len(indices.split()) == 1 and indices.strip() in cc for _, indices in accesses)
        input_accesses = re.findall(r"= Op(?:InBounds)?AccessChain %\S+ %inputs\b", body)
        assert not input_accesses, name
        calls = re.findall(r"= OpFunctionCall %\S+ (%\S+)", body)
        call_names = Counter(cn.get(callee, callee) for callee in calls)
        assert not any(name in call_names for name in ("EvaluateOne", "ReadValue", "WriteValue", "ScratchRead", "ScratchWrite"))
        report["literal_functions"][name] = {
            "id": key, "body_bytes": len(body.encode()),
            "literal_shared_access_chains": len(accesses),
            "literal_shared_indices": sorted({cc[indices.strip()] for _, indices in accesses}),
            "direct_input_buffer_access_chains": len(input_accesses),
            "calls": dict(call_names),
            "barriers": len(re.findall(r"\bOpControlBarrier\b", body)),
            "instructions": dict(Counter(re.findall(r"\b(Op\w+)", body)))}
    _, recognition = body_by_name(cf, cn, "RecognizeLiteralTransport")
    report["recognition"] = {"body_bytes": len(recognition.encode()),
        "barriers": len(re.findall(r"\bOpControlBarrier\b", recognition)),
        "atomic_and": len(re.findall(r"\bOpAtomicAnd\b", recognition)),
        "direct_input_buffer_access_chains": len(re.findall(r"= Op(?:InBounds)?AccessChain %\S+ %inputs\b", recognition)),
        "array_access_chains": [line.strip() for line in recognition.splitlines()
                                if "AccessChain" in line],
        "function_variables": [line.strip() for line in recognition.splitlines()
                               if "OpVariable" in line]}
    assert report["recognition"]["barriers"] == 2
    assert report["recognition"]["atomic_and"] == 2
    _, main_body = body_by_name(cf, cn, "ComputeMain")
    blocks = re.split(r"\n\s*(%\S+) = OpLabel\n", main_body)
    routing = []
    for i in range(1, len(blocks), 2):
        for line in blocks[i + 1].splitlines():
            match = re.search(r"= OpFunctionCall %\S+ (%\S+)", line)
            if match and cn.get(match[1]) in ("LiteralGeneralProgram", "LiteralSchwarzschildProgram", "EvaluateProgramGeneric"):
                routing.append({"block": blocks[i], "callee": cn[match[1]], "line": line.strip()})
    assert len(routing) == 3 and len({item["block"] for item in routing}) == 3
    report["routing_call_blocks"] = routing
    report["main_cfg_requires_independent_review"] = True
    report["generic_fallback_retained"] = all(any(cn.get(key) == name for key in cf)
                                               for name in ("EvaluateOne", "EvaluateProgramGeneric", "ReadValue", "WriteValue"))
    assert report["generic_fallback_retained"]
    call_graph = {key: re.findall(r"= OpFunctionCall %\S+ (%\S+)", body) for key, body in cf.items()}
    report["canonical_program_reachable_functions"] = {}
    for family in ("General", "Schwarzschild"):
        entry, _ = body_by_name(cf, cn, "Literal" + family + "Program")
        visited, pending = set(), [entry]
        while pending:
            key = pending.pop()
            if key in visited:
                continue
            visited.add(key)
            pending.extend(call_graph[key])
        named = sorted(cn.get(key, key) for key in visited)
        assert not set(named).intersection(("EvaluateOne", "EvaluateProgramGeneric", "ReadValue", "WriteValue", "ScratchRead", "ScratchWrite", "ScratchAddress", "RecognizeLiteralTransport"))
        report["canonical_program_reachable_functions"][family] = named
    report["declared_storage"] = storage_facts(candidate, cc)
    report["baseline_declared_storage"] = storage_facts(baseline, bc)
    assert report["declared_storage"]["declared_workgroup_member_bytes"] == 9688
    assert report["baseline_declared_storage"]["declared_workgroup_member_bytes"] == 9684
    programs = json.loads((WORK / "programs.json").read_text())
    composites = [(m[1], m[2], m[3].split()) for m in re.finditer(r"^\s*(%\S+) = OpConstantComposite (%\S+) ([^\n]+)", candidate, re.M)]
    report["complete_literal_table_constant_joins"] = {}
    for family, program in programs.items():
        expected = [program["instructions"], program["registers"], len(program["outputs"]),
                    *program["outputs"], *program["operations"], *program["layer_offsets"]]
        matches = [key for key, _, values in composites if len(values) == len(expected) and
                   all(value in cc for value in values) and [cc[value] for value in values] == expected]
        assert len(matches) == 1, (family, len(matches))
        report["complete_literal_table_constant_joins"][family] = {"words": len(expected), "constant_composite": matches[0], "exact": True}
    report["module_size_factor"] = report["diagnostic"]["module"]["bytes"] / report["baseline"]["module"]["bytes"]
    report["pass"] = True
    (WORK / "compiled-inspection.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"pass": True, "module_size_factor": report["module_size_factor"],
                      "baseline_bytes": report["baseline"]["module"]["bytes"],
                      "diagnostic_bytes": report["diagnostic"]["module"]["bytes"],
                      "literal_layer_functions": len(layer_names),
                      "literal_shared_accesses": sum(value["literal_shared_access_chains"] for value in report["literal_functions"].values()),
                      "declared_workgroup_member_bytes": report["declared_storage"]["declared_workgroup_member_bytes"],
                      "large_declared_aggregates": report["declared_storage"]["aggregate_variables_at_least_10000_bytes"]}))


if __name__ == "__main__":
    main()
