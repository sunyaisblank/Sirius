"""Compose diagnostic inputs only; ROOT executes after emitted admission.

Usage: python3 -B compose_inputs.py EMITTED_ADMISSION.json
The single admission has schema sirius-native-constant-emitted-admission-v1,
passed/source_revision, compiler_owner/compiler_products/command_results pins,
modules {SPV filename: identity}, and true outside_constructor_semantics_preserved_with_explicit_pure_CSE,
public_ABI_equal, bounded_native_constructor_refinement_validated. Each identity has
path/bytes/sha256. It grants no provider, runtime, benefit or adoption authority.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import struct
import sys

ROOT = Path(__file__).resolve().parents[2]
BASIS = ROOT / "out/native-fp32-packed-scratch-producer-inputs"
COMPILE = ROOT / "out/native-constant-compiler-source/candidate"
OUTPUT = ROOT / "out/native-constant-producer-inputs"
CANONICAL = ROOT / "bin/linux-gcc/src/sirius/backend/retained"
REVISION = "5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e"
TARGETS = {"retained_transport.spv": "kTransportShader",
           "retained_endpoint.spv": "kEndpointShader",
           "retained_dense.spv": "kDenseShader",
           "retained_initialize.spv": "kInitializeShader",
           "retained_dopri_phase.spv": "kDopriPhaseShader"}
ARRAY = re.compile(rb"inline constexpr std::array<std::uint32_t, ([0-9]+)> "
                   rb"([A-Za-z0-9_]+)\{\{\n(.*?)\n\}\};\n", re.S)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(data):
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def read(path, expected=None):
    require(path.is_file() and not path.is_symlink(), f"not a regular file: {path}")
    data = path.read_bytes()
    if expected is not None:
        require(identity(data) == {k: expected[k] for k in ("bytes", "sha256")},
                f"identity mismatch: {path}")
        if "path" in expected:
            reference = Path(expected["path"])
            require((reference if reference.is_absolute() else ROOT / reference).resolve()
                    == path.resolve(), f"identity path mismatch: {path}")
    return data


def arrays(header, count):
    result = {}
    for match in ARRAY.finditer(header):
        name = match[2].decode("ascii")
        require(name not in result, f"duplicate array: {name}")
        require(re.fullmatch(rb"(?:[0-9]+u,\s*)+", match[3]) is not None,
                f"unexpected array syntax: {name}")
        words = [int(token) for token in re.findall(rb"([0-9]+)u", match[3])]
        require(len(words) == int(match[1]) and all(0 <= w <= 0xffffffff for w in words),
                f"invalid array words: {name}")
        result[name] = (match, struct.pack("<" + "I" * len(words), *words))
    require(len(result) == count, f"expected {count} arrays")
    return result


def format_array(name, data):
    # Exact scripts/build-retained-kernels.py array() decimal-u/16-word format.
    words = struct.unpack("<" + "I" * (len(data) // 4), data)
    lines = [f"inline constexpr std::array<std::uint32_t, {len(words)}> {name}{{{{"]
    lines += [",".join(str(v) + "u" for v in words[start:start + 16]) + ","
              for start in range(0, len(words), 16)]
    return ("\n".join(lines) + "\n}};\n").encode("ascii")


def main():
    require(len(sys.argv) == 2, "expected one emitted-admission JSON path")
    require(not OUTPUT.exists(), "output already exists; do not overwrite")
    admission_path = Path(sys.argv[1]).resolve(strict=True)
    require(admission_path.is_relative_to(ROOT), "admission must be inside this checkout")
    admission_bytes = read(admission_path)
    admission = json.loads(admission_bytes)
    require(admission["schema"] == "sirius-native-constant-emitted-admission-v1"
            and admission["passed"] is True and admission["source_revision"] == REVISION,
            "emitted admission refused")
    for field in ("outside_constructor_semantics_preserved_with_explicit_pure_CSE", "public_ABI_equal",
                  "bounded_native_constructor_refinement_validated"):
        require(admission[field] is True, f"missing emitted check: {field}")
    owner_bytes = read(COMPILE / "owner.json", admission["compiler_owner"])
    product_bytes = read(COMPILE / "products.json", admission["compiler_products"])
    command_bytes = read(COMPILE / "command-results.json", admission["command_results"])
    owner, products, commands = map(json.loads, (owner_bytes, product_bytes, command_bytes))
    require(owner["passed"] is True and owner["status"] == "terminal"
            and type(owner["child_exit"]) is int and owner["child_exit"] == 0
            and owner["stop_reason"] is None and owner["cleanup_errors"] == []
            and owner["remaining_owned_processes"] == []
            and owner["source_unchanged"] is True and owner["whole_inputs_unchanged"] is True
            and owner["source"] == {"revision": REVISION, "status": ""}, "compiler owner refused")
    require(products["schema"] == "sirius-native-constant-five-compiler-products-v1"
            and products["passed"] is True and products["actual_tool_calls"] == 25
            and products["optimizer_invoked"] is False and products["runtime_admitted"] is False,
            "compiler products refused")
    require(set(products["products"]) == {n.removesuffix(".spv") for n in TARGETS},
            "compiler product inventory differs")
    require(len(commands) == 25 and all(type(r["exit_code"]) is int and r["exit_code"] == 0
                                       for r in commands), "compiler calls incomplete")
    read(COMPILE.parent / "source/src/sirius/kernels/retained_triple.slang", {
        "bytes": 15525, "sha256": "c035d93912824a7321a30700433aa20911408a4734de50b44cad15971db10426"})
    read(COMPILE.parent / "source/src/sirius/kernels/retained_scratch.slang", {
        "bytes": 892, "sha256": "cd4a8c27cd2936d78f893c1d825eb988fb9c7f993433a5756d24e1b423626ed8"})
    read(COMPILE.parent / "recipes.json", products["recipe"])

    prior_bytes = read(BASIS / "root-composed-inputs.json", {
        "bytes": 30848, "sha256": "5d4062652bcad0e10807c8195b6279d05f9fd5a2a9dbad592424d1f977da5922"})
    prior = json.loads(prior_bytes)
    require(prior["passed"] is True and prior["source_revision"] == REVISION,
            "baseline composition authority refused")
    files = {name: read(BASIS / name, row) for name, row in prior["files"].items()
             if name.startswith("baseline/") or name == "observer/audited-observer.cpp"}
    require(len(files) == 39, "expected38 baseline files and observer")
    for name, data in files.items():
        if name.startswith("baseline/"):
            require(read(CANONICAL / name.removeprefix("baseline/").removeprefix("shaders/"))
                    == data, f"current canonical baseline differs: {name}")
    require(identity(files["observer/audited-observer.cpp"]) == {
        "bytes": 342470, "sha256": "c9485d03ec55b5c1307a34097f15efbc58c4f84749e9f3736ad063a68676b601"},
        "observer differs from inherited27cases/1492assertions authority")
    header = files["baseline/retained_kernels.h"]
    main_arrays = arrays(header, 40)
    owner_arrays = arrays(files["baseline/retained_owner_kernels.h"], 2)
    kinds = {"camera": "Camera", "transport": "Transport", "endpoint": "Endpoint",
             "dense": "Dense", "initialize": "Initialize", "ray_camera": "RayCamera",
             "dopri_phase": "DopriPhase"}
    shader_map = {}
    for stem, kind in kinds.items():
        for suffix, variant in (("", ""), ("_fp64", "Fp64"), ("_portable", "Portable"),
                                ("_portable_fp64", "PortableFp64"), ("_fma", "Fma"),
                                ("_portable_normal_sum", "PortableNormalSum")):
            filename = f"retained_{stem}{suffix}.spv"
            if "baseline/shaders/" + filename in files:
                shader_map[filename] = main_arrays[f"k{kind}{variant}Shader"][1]
    for stem, kind in (("transport", "Transport"), ("endpoint", "Endpoint")):
        shader_map[f"retained_{stem}_portable_stage7_owner.spv"] = owner_arrays[f"k{kind}Shader"][1]
    require(len(shader_map) == 35 and {f"baseline/shaders/{n}" for n in shader_map}
            == {n for n in files if n.startswith("baseline/shaders/")}, "shader inventory differs")
    require(all(files["baseline/shaders/" + name] == data for name, data in shader_map.items()),
            "baseline header/SPV payload mismatch")
    require({n for n in main_arrays if n.endswith("Program")}
            == {f"k{k}Program" for k in kinds.values()}, "seven table inventory differs")

    require(set(admission["modules"]) == set(TARGETS), "emitted module inventory differs")
    replacements = {}
    for filename, array_name in TARGETS.items():
        path = COMPILE / filename
        data = read(path, admission["modules"][filename])
        read(path, products["products"][filename.removesuffix(".spv")][".spv"])
        require(len(data) >= 20 and len(data) % 4 == 0 and data[:4] == b"\x03\x02\x23\x07",
                f"invalid SPIR-V framing: {filename}")
        require(data != shader_map[filename], f"selected module unchanged: {filename}")
        replacements[array_name] = format_array(array_name, data)
    candidate = header
    for name in sorted(replacements, key=lambda n: main_arrays[n][0].start(), reverse=True):
        match = main_arrays[name][0]
        candidate = candidate[:match.start()] + replacements[name] + candidate[match.end():]
    candidate_arrays = arrays(candidate, 40)
    require(set(candidate_arrays) == set(main_arrays), "candidate array names differ")
    restored = candidate
    for name in sorted(replacements, key=lambda n: candidate_arrays[n][0].start(), reverse=True):
        match = candidate_arrays[name][0]
        restored = restored[:match.start()] + main_arrays[name][0][0] + restored[match.end():]
    require(restored == header, "text outside five array blocks changed")
    require(sum(main_arrays[n][0][0] == candidate_arrays[n][0][0] for n in main_arrays) == 35,
            "expected35 byte-exact arrays including seven tables")
    for name, data in list(files.items()):
        if name.startswith("baseline/"):
            files[name.replace("baseline/", "candidate/", 1)] = data
    files["candidate/retained_kernels.h"] = candidate
    for filename, array_name in TARGETS.items():
        files["candidate/shaders/" + filename] = candidate_arrays[array_name][1]
    require(len(files) == 77, "expected77 payload files")
    require(sum(files["candidate/shaders/" + n] == data for n, data in shader_map.items()) == 30,
            "excluded30 shader payloads changed")
    OUTPUT.mkdir()
    for name, data in files.items():
        path = OUTPUT / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(data)
        path.chmod(0o644)
    report = {"schema": "sirius-native-constant-diagnostic-input-composition-v1", "passed": True,
              "created_utc": datetime.now(timezone.utc).isoformat(), "source_revision": REVISION,
              "prior_baseline_authority": identity(prior_bytes), "emitted_admission": {
                  "path": str(admission_path.relative_to(ROOT)), **identity(admission_bytes)},
              "compiler_owner": identity(owner_bytes), "compiler_products": identity(product_bytes),
              "command_results": identity(command_bytes), "target_array_names": TARGETS,
              "main_arrays": {"changed": 5, "exact": 35, "seven_tables_and_other_text_exact": True},
              "excluded_SPVs_exact": 30, "both_owner_headers_exact": True,
              "observer_exact": True, "observer_original_cases": 27, "observer_ordered_assertions": 1492,
              "files": {n: {**identity(d), "mode": "0o644"} for n, d in sorted(files.items())},
              "cpp_producer_compiled": False, "runtime_released": False, "adoption_admitted": False,
              "scope": "Five new compiler/emitted-reviewed ordinary nativeFP32 RTConstant SPVs composed with current5b2 canonical unpacked baseline. Diagnostic input composition only; no provider/science/benefit/fullqualification authority."}
    with (OUTPUT / "root-composed-inputs.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
