"""One owned pure reader of completed immutable Endpoint compiler products."""
from pathlib import Path
import hashlib, importlib.util, json, os, re, subprocess, sys
sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[2]
P = json.loads((WORK / "inputs.json").read_text())


def record(path):
    path = Path(path)
    return {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "resolved_path": str(path.resolve(strict=True))}


def content(path):
    return {key: value for key, value in record(path).items() if key != "resolved_path"}


def load(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    assert not path.exists(), path
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def verify_frozen():
    for name, pin in load(WORK / "source-freeze.json")["pins"].items():
        assert record(name) == pin, name


def original_path(path):
    path = Path(path)
    if not path.is_absolute(): path = ROOT / path
    archived = Path(P["archive"]) / "files" / path.relative_to(ROOT) if path.is_relative_to(ROOT) else None
    if archived is not None and archived.is_file(): return archived
    assert not path.is_relative_to(ROOT / "out/portable-divide-fp64-estimate-screen/compiler")
    assert not path.is_relative_to(ROOT / "bin/linux-gcc/portable-divide-fp64-estimate-screen/compiler")
    return path


def original_completed():
    archive = Path(P["archive"])
    assert content(archive / "manifest.json")["sha256"] == P["manifest_sha256"]
    manifest = load(archive / "manifest.json")
    for name, expected in manifest["files"].items():
        assert content(archive / name) == expected
    old = archive / "files/out/portable-divide-fp64-estimate-screen/compiler"
    owner = load(old / "compile-endpoint-normal/owner.json")
    assert owner["status"] == "terminal" and owner["passed"] is False and owner["child_exit"] == 1
    assert owner["stop_reason"] is None and owner["cleanup_errors"] == []
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["remaining_owned_processes"] == [] and len(owner["observed_births"]) == 2
    assert all(row["absent"] is True for row in owner["observed_births"])
    assert owner["source"] == {"revision": "4498cc41168ad29fc4c111819a4cb7ccaa54a8bf", "status": ""}
    assert owner["argv"] == ["/usr/bin/python3", "-B", str(ROOT / "out/portable-divide-fp64-estimate-screen/compiler/screen.py"), "compile"]
    assert content(old / "compile-endpoint-normal/stderr.log") == P["original_stderr"]
    assert not (old / "compiler-result.json").exists() and not (old / "compiler-acceptance.json").exists()
    before_path = old / "compile-endpoint-normal/inputs-before.json"
    assert before_path.read_bytes() == (old / "compile-endpoint-normal/inputs-after.json").read_bytes()
    before = load(before_path)
    frozen = load(old / "source-freeze.json")
    for kind in ("private_inputs", "external_inputs"):
        for name, pin in frozen[kind].items():
            original = ROOT / "out/portable-divide-fp64-estimate-screen/compiler" / name if kind == "private_inputs" else Path(name)
            assert content(original_path(original)) == {key: pin[key] for key in ("bytes", "sha256")}
    action = load(old / "actions.json")["compile-endpoint-normal"]
    for name in action["frozen_inputs"]:
        original = Path(name); key = str(original.relative_to(ROOT)) if original.is_relative_to(ROOT) else str(original)
        actual = original_path(original)
        assert content(actual) == {field: before[key][field] for field in ("bytes", "sha256")}
        if actual == original:
            assert record(actual) == before[key]
        else:
            # These immutable relocated entries were regular private files,
            # not the five lexical tool aliases retained at their real paths.
            assert before[key]["resolved_path"] == str(original)
    header = old / "candidate/src/sirius/kernels/portable_binary32.h"
    insertion = (old / "candidate-insertion.txt").read_text()
    assert header.read_text().count(insertion) == 1
    assert header.read_text().replace(insertion, "", 1) == (old / "stock/src/sirius/kernels/portable_binary32.h").read_text()
    assert (old / "stock/src/sirius/kernels/portable_binary32.h").read_bytes() == (ROOT / "src/sirius/kernels/portable_binary32.h").read_bytes()
    return old, {"original_failed_wrapper": record(old / "compile-endpoint-normal/owner.json"),
                 "original_passed": False, "original_child_exit": 1,
                 "exact_post_validation_reader_failure": True,
                 "archive_manifest": record(archive / "manifest.json"),
                 "original_input_byte_joins": len(action["frozen_inputs"]),
                 "no_compiler_validator_or_policy_rerun": True}


def pin_context():
    verify_frozen(); original_completed()
    assert len(sys.argv) == 3 and re.fullmatch(r"[0-9a-f]{40}", sys.argv[2])
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip() == sys.argv[2]
    assert not subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT)
    assert not (WORK / "operator-context.json").exists() and not (WORK / "read-completed-endpoint").exists()
    assert (WORK / "actions.json").read_bytes() == (WORK / "actions.unpinned.json").read_bytes()
    actions = load(WORK / "actions.json")
    assert list(actions) == ["read-completed-endpoint"]
    assert actions["read-completed-endpoint"]["revision"] == "PENDING_OPERATOR_CLEAN_CONTEXT"
    actions["read-completed-endpoint"]["revision"] = sys.argv[2]
    (WORK / "actions.json").write_text(json.dumps(actions, sort_keys=True, indent=2) + "\n")
    write(WORK / "operator-context.json", {"revision": sys.argv[2], "actions": record(WORK / "actions.json"),
          "source_freeze": record(WORK / "source-freeze.json"), "preparation": record(WORK / "preparation.json")})


def read_products():
    verify_frozen(); old, producer = original_completed()
    context = load(WORK / "operator-context.json")
    assert context["actions"] == record(WORK / "actions.json") and context["source_freeze"] == record(WORK / "source-freeze.json")
    assert context["preparation"] == record(WORK / "preparation.json")
    original = module("preserved_compiler_pure_functions", old / "screen.py")
    graph = module("corrected_estimate_graph", WORK / "estimate_graph.py")
    controls = module("complete_module_parser_controls", WORK / "parser_controls.py")
    products = Path(P["products"])
    texts = {mode: (products / (mode + "-endpoint-normal.final.spvasm")).read_text() for mode in ("stock", "candidate")}
    shape = graph.inspect(texts["stock"], texts["candidate"], original.function_body)
    checked = controls.controls(texts["stock"], texts["candidate"], graph, original.function_body)
    comparisons = {}
    for name in ("PB32Unpack", "PB32Round", "PB32Multiply", "PB32WideProduct", "PB32ProductResidual"):
        present = [bool(re.search(r'OpName %\S+ "' + name + '"', texts[mode])) for mode in ("stock", "candidate")]
        assert present[0] == present[1]
        if not present[0]:
            assert name == "PB32WideProduct"
            comparisons[name] = "No named function in either module; unchanged source intrinsic and actual OpUMulExtended graph checked."
        else:
            assert original.normalized_function(texts["stock"], name) == original.normalized_function(texts["candidate"], name)
            comparisons[name] = "Complete normalized function equal."
    assert original.interface_signature(texts["stock"]) == original.interface_signature(texts["candidate"])
    modes = lambda text: sorted(re.findall(r"OpExecutionMode %\S+ (.*)$", text, re.M))
    assert modes(texts["stock"]) == modes(texts["candidate"])
    def barriers(text):
        definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
        return [[opcode, *[definitions[word] for word in operands.split()]]
                for opcode, operands in re.findall(r"^\s*(OpControlBarrier|OpMemoryBarrier) (.+)$", text, re.M)]
    assert barriers(texts["stock"]) == barriers(texts["candidate"])
    assert not (WORK / "products").exists(); (WORK / "products").mkdir()
    modules = {}
    for mode in ("stock", "candidate"):
        text = texts[mode]
        assert sorted(re.findall(r"OpCapability (\S+)", text)) == sorted(["Shader", "RoundingModeRTE"] + (["Float64"] if mode == "candidate" else []))
        widths = dict(re.findall(r"(%\S+) = OpTypeInt (\d+) [01]", text)); assert set(widths.values()) == {"32"}
        assert sorted(re.findall(r"OpTypeFloat (\d+)", text)) == (["32", "64"] if mode == "candidate" else ["32"])
        assert "OpConvertFToU" not in text and "OpConvertFToS" not in text
        assert original.function_body(text, "PB32Divide").count("OpLoopMerge") == 1
        assert all(widths[kind] == "32" for kind in re.findall(r"= OpUDiv (%\S+) ", text))
        bodies = {}
        for name in ("PB32Multiply", "PB32Divide"):
            path = WORK / "products" / (mode + "." + name + ".spvasm")
            path.write_text(original.function_body(text, name)); bodies[name] = record(path)
        modules[mode] = {"module": record(products / (mode + "-endpoint-normal.spv")),
                         "assembly": record(products / (mode + "-endpoint-normal.final.spvasm")), "bodies": bodies}
    assert Path(modules["stock"]["module"]["resolved_path"]).read_bytes() == (ROOT / "bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv").read_bytes()
    write(WORK / "parser-control-result.json", checked)
    write(WORK / "compiler-result.json", {"passed": True, "producer": producer, "modules": modules,
          "candidate_estimate_graph": shape, "unchanged_function_comparisons": comparisons,
          "original_multiply_residual_unpack_round_normalized_bodies_equal": True,
          "original_fallback_complete_suffix_normalized_equal": True, "candidate_has_no_Int64_or_FToU": True,
          "complete_descriptor_and_workgroup_type_decorations_exact": True, "execution_modes_exact": True,
          "ordered_barrier_operands_exact": True, "runtime_or_numerical_admission": False,
          "parser_controls": record(WORK / "parser-control-result.json"),
          "scope": "Fresh parse-only acceptance of immutable completed compiler products. Original failed wrapper remains false; no recompilation, shader validation rerun, arithmetic or cost verdict."})


def review():
    verify_frozen(); original_completed()
    folder = WORK / "read-completed-endpoint"; owner = load(folder / "owner.json")
    assert owner["status"] == "terminal" and owner["passed"] is True and owner["child_exit"] == 0
    assert owner["stop_reason"] is None and owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["observed_births"] and all(row["absent"] is True for row in owner["observed_births"])
    assert (folder / "inputs-before.json").read_bytes() == (folder / "inputs-after.json").read_bytes()
    context = load(WORK / "operator-context.json"); action = load(WORK / "actions.json")["read-completed-endpoint"]
    assert owner["source"] == {"revision": context["revision"], "status": ""}
    assert owner["argv"] == action["argv"] and owner["limits"] == action["limits"]
    assert context["actions"] == record(WORK / "actions.json") and context["source_freeze"] == record(WORK / "source-freeze.json")
    before = load(folder / "inputs-before.json")
    for name in action["frozen_inputs"]:
        path = Path(name); key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        assert before[key] == record(path)
    result = load(WORK / "compiler-result.json"); assert result["passed"] and result["producer"]["original_passed"] is False
    for values in result["modules"].values():
        for pin in (values["module"], values["assembly"], *values["bodies"].values()): assert record(pin["resolved_path"]) == pin
    assert record(WORK / "parser-control-result.json") == result["parser_controls"]
    write(WORK / "compiler-acceptance.json", {"passed": True, "reader_owner": record(folder / "owner.json"),
          "compiler_result": record(WORK / "compiler-result.json"), "original_failed_wrapper": result["producer"]["original_failed_wrapper"],
          "original_wrapper_passed": False, "original_wrapper_child_exit": 1,
          "scope": "Separate terminal pure-reader acceptance; original post-validation wrapper failure preserved. No compiler/provider/numerical execution or benefit claim."})


if __name__ == "__main__":
    {"pin-context": pin_context, "read": read_products, "review": review}[sys.argv[1]]()
