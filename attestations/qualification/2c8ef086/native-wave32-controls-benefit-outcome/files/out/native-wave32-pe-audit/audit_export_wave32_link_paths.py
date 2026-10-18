"""Static reader for the one matched Windows export; never execute its products.

Header decoding is retained from f851 read_postproof.py::arrays. Windows source/
artifact joins follow the selected-record approach of 9d56406 reconstruction,
without its project imports, staging, subprocesses or native-gate authority.
Wave32 adds two private Vulkan objects and identical shader payloads. LINK input
paths are resolved against the recorded Windows build directory; caret headers
are excluded. COFF/MAP/VERBOSE records establish selection, not textual order. This source
has not been executed. ROOT must review it before use.
"""
from pathlib import Path
import argparse
import hashlib
import json
import ntpath
import re
import struct
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "out/native-wave32-producer-inputs"
CANONICAL_CXX_REVISION = "2c8ef08693d46f57d5dea14efdbddd4874dac684"
NS = {"m": "http://schemas.microsoft.com/developer/msbuild/2003"}


def require(value, message):
    if not value:
        raise ValueError(message)


def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def document(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def arrays(path):
    result = {}
    pattern = r"inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};"
    for count, name, body in re.findall(pattern, path.read_text(), re.S):
        words = [int(word) for word in re.findall(r"(\d+)u", body)]
        require(len(words) == int(count) and name not in result, "Malformed/duplicate header array")
        result[name] = struct.pack("<" + str(len(words)) + "I", *words)
    return result


def norm(value):
    return str(value).replace("\\", "/").lower()


def text(path):
    data = path.read_bytes()
    return data.decode("utf-16" if data.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8-sig")


def logs(base, prefix):
    return sorted(p for p in base.rglob("*.tlog") if p.name.lower().startswith(prefix.lower()))


def check_flags(command, label):
    for option in ("/std:c++latest", "/W4", "/WX", "/permissive-", "/utf-8"):
        require(option.lower() in command.lower(), "Missing strict flag: " + label + " " + option)
    require(re.search(r"SIRIUS_CONTRACT_MODE=2(?![0-9])", command), "Missing exact contract mode2: " + label)


def source_records(paths, suffix):
    result = []
    for path in paths:
        selected = False
        for line in text(path).splitlines():
            if line.startswith("^"):
                selected = any(norm(v).endswith(norm(suffix)) for v in line.split("^") if v)
            if selected:
                result.append(line)
    require(result, "No actual tlog record for " + suffix)
    return "\n".join(result)


def archive_member(archive, object_path):
    data = archive.read_bytes()
    expected = object_path.read_bytes()
    require(data.startswith(b"!<arch>\n"), "Not COFF archive")
    cursor, matches = 8, []
    while cursor < len(data):
        require(cursor + 60 <= len(data) and data[cursor + 58:cursor + 60] == b"`\n", "Bad archive member")
        length = int(data[cursor + 48:cursor + 58])
        require(length >= 0, "Negative archive member length")
        begin = cursor + 60
        require(begin + length <= len(data), "Truncated archive member")
        if data[begin:begin + length] == expected:
            matches.append(data[cursor:cursor + 16].decode("ascii").strip())
        cursor = begin + length + (length & 1)
    require(cursor == len(data) and len(matches) == 1, "Device object not uniquely present in backend archive")
    return matches[0]


def selected_row(rows, suffix, expected):
    found = [v for k, v in rows.items() if norm(k).endswith(norm(suffix))]
    require(found, "Missing sealed input: " + suffix)
    require(all({k: v[k] for k in ("bytes", "sha256")} == expected for v in found),
            "Sealed input differs: " + suffix)


def pe_sections(path):
    data = path.read_bytes()
    require(len(data) >= 64 and data[:2] == b"MZ", "Not DOS/PE")
    nt = struct.unpack_from("<I", data, 60)[0]
    require(nt + 24 <= len(data) and data[nt:nt + 4] == b"PE\0\0", "Bad PE header")
    machine, count = struct.unpack_from("<HH", data, nt + 4)
    opt_size = struct.unpack_from("<H", data, nt + 20)[0]
    opt = nt + 24
    require(machine == 0x8664 and opt_size >= 112 and opt + opt_size <= len(data), "Not x64 PE32+")
    require(struct.unpack_from("<H", data, opt)[0] == 0x20b, "Not PE32+")
    image_base = struct.unpack_from("<Q", data, opt + 24)[0]
    begin = opt + opt_size
    require(begin + count * 40 <= len(data), "Truncated section table")
    sections = []
    for index in range(count):
        name, virtual_size, rva, size, offset, _, _, _, _, flags = struct.unpack_from(
            "<8sIIIIIIHHI", data, begin + index * 40)
        require(offset + size <= len(data), "Truncated PE section")
        if flags & 0x40 and flags & 0x40000000 and not flags & 0x20000000:
            initialized = min(size, virtual_size) if virtual_size else size
            sections.append((name.rstrip(b"\0").decode("ascii"), rva, offset,
                             data[offset:offset + initialized]))
    require(sections, "No initialized non-code data section")
    return image_base, sections


def verify_arrays(exe, map_text, namespace, expected, image_base, sections):
    result = {}
    for name, payload in expected.items():
        require(payload, "Empty expected array: " + name)
        locations = []
        for section, rva, raw_offset, data in sections:
            start = 0
            while True:
                offset = data.find(payload, start)
                if offset < 0:
                    break
                locations.append({"section": section, "file_offset": raw_offset + offset,
                                  "rva": rva + offset})
                start = offset + 1
        require(locations, "Expected full PE array missing: " + name)
        symbol_prefix = "?" + name + "@" + namespace + "@backend@sirius@@"
        addresses = []
        for line in map_text.splitlines():
            match = re.match(r"\s*[0-9a-fA-F]{4}:[0-9a-fA-F]+\s+(\S+)\s+([0-9a-fA-F]{8,16})\b", line)
            if match and match[1].startswith(symbol_prefix):
                addresses.append(int(match[2], 16) - image_base)
        require(addresses and set(addresses) <= {v["rva"] for v in locations},
                "Missing/ambiguous map-to-PE array join: " + name)
        result[name] = {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(),
                        "locations": locations, "map_rvas": sorted(set(addresses))}
    return result


def coff_functions(path):
    """Read ordinary AMD64 COFF external function definitions; fail on other formats."""
    data = path.read_bytes()
    require(len(data) >= 20, "Truncated COFF object")
    machine, sections, _, pointer, count, optional, _ = struct.unpack_from("<HHIIIHH", data)
    require(machine == 0x8664 and sections > 0 and optional == 0,
            "Expected ordinary x64 COFF object; BigObj/LTCG needs separate review")
    strings = pointer + count * 18
    require(pointer >= 20 + sections * 40 and strings + 4 <= len(data), "COFF symbol bounds")
    length = struct.unpack_from("<I", data, strings)[0]
    require(length >= 4 and strings + length <= len(data), "COFF string bounds")
    result, index = {}, 0
    while index < count:
        entry = pointer + index * 18
        raw = data[entry:entry + 8]
        if raw[:4] == b"\0" * 4:
            offset = struct.unpack_from("<I", raw, 4)[0]
            require(4 <= offset < length, "COFF name offset")
            end = data.find(b"\0", strings + offset, strings + length)
            require(end >= 0, "Unterminated COFF symbol")
            name = data[strings + offset:end].decode("ascii")
        else:
            name = raw.rstrip(b"\0").decode("ascii")
        value, section, kind, storage, auxiliary = struct.unpack_from("<IhHBB", data, entry + 8)
        require(index + auxiliary < count, "COFF auxiliary symbol bounds")
        if storage == 2 and kind & 0x20 and section > 0:
            require(section <= sections and name not in result, "COFF duplicate/invalid function")
            result[name] = {"section": section, "value": value}
        index += auxiliary + 1
    require(result, "No external COFF function definitions")
    return result


def selected_command_semantics(command):
    """Compare actual non-path compiler switches; retain includes for manual audit."""
    body = " ".join(line for line in command.splitlines() if not line.startswith("^"))
    body = re.sub(r'(?i)(?<!\S)/(?:external:)?I\s*(?:"[^"\r\n]*"|\S+)', " ", body)
    body = re.sub(r'(?i)(?<!\S)/F[od]\s*(?:"[^"\r\n]*"|\S+)', " ", body)
    body = re.sub(r'(?i)(?:[A-Z]:[^\s"]*|"[A-Z]:[^"\r\n]*")\.cpp(?=\s|$)', " ", body)
    # The remaining switches include all /D feature defines and arithmetic/ABI flags.
    return " ".join(body.lower().split())


def map_device_functions(mapping, obj, arm):
    definitions = coff_functions(obj)
    result = {}
    for prefix in ("?CreateVulkanDevice@backend@sirius@@", "?EnumerateVulkanDevices@backend@sirius@@",
                   "?LoadKernel@VulkanDevice@backend@sirius@@",
                   "?GetOrCreatePipeline@VulkanDevice@backend@sirius@@"):
        rows = []
        for line in mapping.splitlines():
            match = re.match(r"\s*[0-9a-fA-F]{4}:[0-9a-fA-F]+\s+(\S+)\s+([0-9a-fA-F]{8,16})\b(.*)", line)
            if match and match[1].startswith(prefix):
                owner = match[3].split()[-1]
                require(norm(owner) == "vulkan_device.obj", "Device definition is not a direct private object: " + arm)
                require(match[1] in definitions, "MAP definition absent from selected private COFF object")
                rows.append({"symbol": match[1], "address": match[2], "map_object": owner,
                             "coff_definition": definitions[match[1]]})
        require(len(rows) == 1, "Missing/ambiguous private device method MAP join: " + prefix)
        result[prefix] = rows[0]
    return result


def windows_file(value, build_root):
    """Resolve a selected explicit file operand in the recorded MSBuild directory."""
    require(ntpath.isabs(build_root) and ntpath.splitdrive(build_root)[0], "Not an absolute Windows build root")
    return ntpath.normcase(ntpath.normpath(value if ntpath.isabs(value) else ntpath.join(build_root, value)))


def link_command_files(records, build_root):
    # A caret header lists tlog tracking keys, not LINK command arguments.
    command = "\n".join(line for line in records.splitlines() if not line.startswith("^"))
    require(not re.search(r"(?i)(?<!\S)/(?:force|wholearchive)(?=[:\s]|$)", command),
            "Unexpected FORCE or WHOLEARCHIVE")
    paths = set()
    # Current selected object/archive operands are standalone bare or quoted
    # path arguments. Output switches (e.g. /IMPLIB:"...") are not operands.
    for quoted, bare in re.findall(r'(?<!\S)(?:"([^"\r\n]+)"|([^\s"]+))(?=\s|$)', command):
        value = quoted or bare
        if value.startswith(("/", "-")) or not value.lower().endswith((".obj", ".lib")):
            continue
        # Bare library names use LINK search rules. Do not invent their resolved
        # paths; only explicit selected archive/object paths are joined below.
        if value.lower().endswith(".lib") and not ntpath.dirname(value):
            continue
        paths.add(windows_file(value, build_root))
    return command, paths


def link_read_files(records, build_root):
    return {windows_file(line.strip(), build_root) for line in records.splitlines()
            if not line.startswith("^") and line.strip().lower().endswith((".obj", ".lib"))}


def main(args):
    require(re.fullmatch(r"[0-9a-f]{40}", args.revision), "Expected exact private revision")
    require(re.fullmatch(r"[0-9]+", args.run_id), "Expected exact run id")
    export = args.export.resolve()
    require(not args.output.resolve().is_relative_to(export), "Output must be outside immutable export")
    manifest = document(export / "manifest.json")
    require(manifest["schema"] == "sirius-native-packed-ci-export-manifest-v1" and manifest["passed"] is True,
            "No successful producer manifest")
    require(manifest["source_revision"] == args.revision, "Manifest revision mismatch")
    names = set()
    for relative, row in manifest["files"].items():
        p = Path(relative)
        require(not p.is_absolute() and ".." not in p.parts and "\\" not in relative and ":" not in relative,
                "Unsafe export member")
        require(relative.lower() not in names, "Case-colliding export member")
        names.add(relative.lower())
        path = export / p
        require(path.is_file() and not path.is_symlink(), "Missing/aliased member")
        require(identity(path) == {k: row[k] for k in ("bytes", "sha256")}, "Member differs: " + relative)
    actual = {str(p.relative_to(export)).replace("\\", "/").lower()
              for p in export.rglob("*") if p.is_file() and p != export / "manifest.json"}
    require(actual == names, "Unlisted/missing export files")
    result = document(export / "producer-result.json")
    for key in ("passed", "configure_passed", "build_passed", "source_unchanged", "frozen_build_inputs_unchanged"):
        require(result[key] is True, "Producer did not pass: " + key)
    require(result["source_revision"] == result["github_sha"] == args.revision, "Producer revision mismatch")
    require(result["repository"] == "sunyaisblank/Sirius" and str(result["run_id"]) == args.run_id,
            "Producer run/repository mismatch")
    require(str(result["run_attempt"]) == "1", "Unexpected run attempt")
    require(result["schema"] == "sirius-native-packed-matched-ci-producer-v1" and
            result["producer_profile"] == "native-wave32-transport-device-only" and
            result["shader_bytes_identical_between_arms"] is True and
            type(result["expected_private_compile_objects"]) is int and
            result["expected_private_compile_objects"] == 7 and
            result["private_compile_object_count"] == 7, "Wrong wave32 producer/object profile")
    before = document(export / "build-inputs-before.json")
    require(before == document(export / "build-inputs-after.json"), "Build seals differ")
    source = document(export / "source-before-configure.json")
    for path in sorted(p for p in INPUTS.rglob("*") if p.is_file()):
        require(not path.is_symlink(), "Aliased local producer input")
        suffix = str(path.relative_to(ROOT)).replace("\\", "/")
        selected_row(source, suffix, identity(path))
        selected_row(before, suffix, identity(path))
    require(identity(INPUTS / "observer/audited-observer.cpp")["sha256"] ==
            "c9485d03ec55b5c1307a34097f15efbc58c4f84749e9f3736ad063a68676b601", "Observer changed")
    require((INPUTS / "baseline/vulkan_device.cpp").read_bytes() ==
            (ROOT / "src/sirius/backend/vulkan/vulkan_device.cpp").read_bytes(), "Baseline Vulkan source is not canonical")
    require((INPUTS / "candidate/vulkan_device.cpp").read_bytes() ==
            (ROOT / "out/native-wave32-source/vulkan_device.cpp").read_bytes(), "Candidate Vulkan source mismatch")
    for header in ("retained_kernels.h", "retained_owner_kernels.h", "retained_owner_layout.h"):
        require((INPUTS / "baseline" / header).read_bytes() == (INPUTS / "candidate" / header).read_bytes(),
                "Full retained header differs between arms")
    shaders = {arm: {p.name: identity(p) for p in (INPUTS / arm / "shaders").glob("*.spv")}
               for arm in ("baseline", "candidate")}
    require(len(shaders["baseline"]) == 35 and shaders["baseline"] == shaders["candidate"], "Expected identical35 SPVs")
    for arm in ("baseline", "candidate"):
        for header in ("retained_kernels.h", "retained_owner_kernels.h", "retained_owner_layout.h"):
            suffix = "out/native-wave32-producer-inputs/" + arm + "/" + header
            selected_row(source, suffix, identity(INPUTS / arm / header))
            selected_row(before, suffix, identity(INPUTS / arm / header))
    for relative in ("src/sirius/backend/retained_compute.cpp", "tests/backend/retained_dopri_test.cpp",
                     "src/sirius/backend/vulkan/vulkan_device.cpp", "cmake/sirius_toolchain.cmake"):
        selected_row(source, relative, identity(ROOT / relative))
    selected_row(source, "out/native-wave32-producer-inputs/observer/audited-observer.cpp",
                 identity(INPUTS / "observer/audited-observer.cpp"))
    expected = {arm: arrays(INPUTS / arm / "retained_kernels.h") for arm in ("baseline", "candidate")}
    require(all(len(v) == 40 for v in expected.values()), "Expected40 main arrays")
    changed = [name for name in expected["baseline"] if expected["baseline"][name] != expected["candidate"][name]]
    require(not changed and expected["baseline"] == expected["candidate"], "Wave32 may not change any main array")
    owner_arrays = arrays(INPUTS / "baseline/retained_owner_kernels.h")
    require(len(owner_arrays) == 2 and owner_arrays == arrays(INPUTS / "candidate/retained_owner_kernels.h"),
            "Owner payload delta")
    products = {}
    for arm in ("baseline", "candidate"):
        exe = export / "payload" / arm / ("original_" + arm + "_controls.exe")
        mapping = text(exe.with_suffix(".map"))
        image_base, sections = pe_sections(exe)
        require(not re.search(r"sirius_backend(?:\.lib)?[:(][^\s)]*retained_compute\.obj", mapping, re.I),
                "Canonical compute object extracted")
        require(not re.search(r"sirius_backend(?:\.lib)?[:(][^\s)]*vulkan_device\.obj", mapping, re.I),
                "Canonical device object extracted instead of private implementation")
        products[arm] = {"PE": identity(exe), "map": identity(exe.with_suffix(".map")),
                         "main_arrays": verify_arrays(exe, mapping, "retained_program", expected[arm], image_base, sections),
                         "owner_arrays": verify_arrays(exe, mapping, "retained_owner", owner_arrays, image_base, sections),
                         "canonical_device_map_members": []}
        for name in changed:
            foreign = expected["candidate" if arm == "baseline" else "baseline"][name]
            require(all(data.find(foreign) < 0 for _, _, _, data in sections),
                    "Foreign changed-arm shader remains in PE: " + arm + " " + name)
    tlogs = {}
    targets = ["native_wave32_dopri"] + ["native_wave32_" + arm + "_" + kind
               for arm in ("baseline", "candidate") for kind in ("compute", "controls")]
    all_objects = []
    for target in targets:
        base = export / "metadata/private-objects" / (target + ".dir") / "Release"
        objects = sorted(base.rglob("*.obj"))
        require(len(objects) == (2 if target.endswith("compute") else 1), "Private object count: " + target)
        all_objects += objects
        commands, reads = logs(base, "CL.command."), logs(base, "CL.read.")
        require(commands and reads, "Missing CL command/read tlogs: " + target)
        project = ET.parse(export / "metadata/projects" / (target + ".vcxproj")).getroot()
        sources = [v.attrib["Include"] for v in project.findall(".//m:ClCompile[@Include]", NS)]
        arm = "baseline" if "baseline" in target else "candidate"
        expected_sources = (["tests/backend/retained_dopri_test.cpp"] if target == "native_wave32_dopri" else
                            ["src/sirius/backend/retained_compute.cpp", "out/native-wave32-producer-inputs/" + arm + "/vulkan_device.cpp"]
                            if target.endswith("compute") else ["out/native-wave32-producer-inputs/observer/audited-observer.cpp"])
        require(len(sources) == len(expected_sources) and
                all(sum(norm(v).endswith(suffix) for v in sources) == 1 for suffix in expected_sources),
                "Unexpected private target sources: " + target)
        records = {}
        for suffix in expected_sources:
            command = source_records(commands, suffix)
            dependencies = norm(source_records(reads, suffix))
            check_flags(command, target + " " + suffix)
            basename = Path(suffix).stem.lower()
            selected_objects = [p for p in objects if p.stem.lower() in (basename, basename + ".cpp")]
            require(len(selected_objects) == 1, "Cannot uniquely join private object by source basename")
            if target == "native_wave32_dopri":
                require("retained_kernels.h" not in dependencies, "Shared DP consumes arm header")
            else:
                other = "candidate" if arm == "baseline" else "baseline"
                headers = (("retained_kernels.h", "retained_owner_layout.h", "retained_owner_kernels.h")
                           if basename == "retained_compute" else ("retained_kernels.h", "retained_owner_layout.h")
                           if not target.endswith("compute") else ("retained_kernels.h",) if arm == "candidate" else ())
                for header in headers:
                    require("out/native-wave32-producer-inputs/" + arm + "/" + header in dependencies,
                            "Expected private header dependency missing")
                require("/" + other + "/retained_kernels.h" not in dependencies, "Other arm header consumed")
            if not target.endswith("compute"):
                require("fail_on_skip.h" in command.lower(), "Missing required-runtime skip refusal")
            if target.endswith("compute"):
                require("SIRIUS_HAS_RETAINED_COMPUTE=1" in command and "SIRIUS_RELEASE_RESOURCE_LOCKED=1" in command,
                        "Missing canonical retained/resource feature define")
            records[suffix] = {"object": identity(selected_objects[0]),
                               "object_path": str(selected_objects[0].relative_to(export)),
                               "command_semantics": selected_command_semantics(command),
                               "source_record": suffix}
        tlogs[target] = {"sources": sources, "records": records,
                         "command_tlogs": [str(p.relative_to(export)) for p in commands],
                         "read_tlogs": [str(p.relative_to(export)) for p in reads]}
    require(len(all_objects) == 7, "Expected exactly seven exported private objects")
    for arm in ("baseline", "candidate"):
        base = export / "metadata/private-objects" / ("native_wave32_" + arm + "_controls.dir") / "Release"
        link_reads, link_commands = logs(base, "link.read."), logs(base, "link.command.")
        require(link_reads and link_commands, "Missing LINK records: " + arm)
        read_records = "\n".join(text(p) for p in link_reads)
        command_records = "\n".join(text(p) for p in link_commands)
        command, command_paths = link_command_files(command_records, result["build_root"])
        read_paths = link_read_files(read_records, result["build_root"])
        dependencies = norm("\n".join(line for line in read_records.splitlines() if not line.startswith("^")))
        require("native_wave32_" + ("candidate" if arm == "baseline" else "baseline") + "_compute.dir" not in dependencies,
                "Other arm object linked")
        archive_path = windows_file("src/sirius/backend/Release/sirius_backend.lib", result["build_root"])
        require(archive_path in read_paths and archive_path in command_paths,
                "Selected backend archive missing from LINK operands/read records")
        for target in ("native_wave32_" + arm + "_compute", "native_wave32_" + arm + "_controls", "native_wave32_dopri"):
            folder = export / "metadata/private-objects" / (target + ".dir") / "Release"
            for obj in folder.rglob("*.obj"):
                suffix = norm(obj.relative_to(export / "metadata/private-objects"))
                selected_path = windows_file(suffix, result["build_root"])
                require(selected_path in read_paths and selected_path in command_paths,
                        "Private object missing from actual LINK operands/read records: " + suffix)
        device_obj = export / tlogs["native_wave32_" + arm + "_compute"]["records"][
            "out/native-wave32-producer-inputs/" + arm + "/vulkan_device.cpp"]["object_path"]
        mapping = text(export / "payload" / arm / ("original_" + arm + "_controls.map"))
        products[arm]["selected_LINK_archive"] = archive_path
        products[arm]["LINK_operands_and_read_records_joined"] = True
        products[arm]["link_textual_order_is_selection_authority"] = False
        products[arm]["private_device_object"] = identity(device_obj)
        products[arm]["private_device_COFF_MAP_joins"] = map_device_functions(mapping, device_obj, arm)
        products[arm]["link_read_tlogs"] = [str(p.relative_to(export)) for p in link_reads]
        products[arm]["link_command_tlogs"] = [str(p.relative_to(export)) for p in link_commands]
    build_log = text(export / "build.stdout.log")
    require(not re.search(r"Loaded[^\r\n]*sirius_backend\.lib\([^)]*vulkan_device\.obj\)", build_log, re.I),
            "VERBOSE:LIB reports canonical Vulkan object extraction")
    backend = export / "metadata/canonical-backend/src/sirius/backend"
    device_root = backend / "sirius_backend.dir/Release"
    devices = [p for p in device_root.rglob("*.obj") if p.name.lower() in ("vulkan_device.obj", "vulkan_device.cpp.obj")]
    require(len(devices) == 1, "Missing/ambiguous canonical device object")
    device_source = "src/sirius/backend/vulkan/vulkan_device.cpp"
    device_commands = source_records(logs(device_root, "CL.command."), device_source)
    device_reads = source_records(logs(device_root, "CL.read."), device_source)
    check_flags(device_commands, "canonical VulkanDevice")
    require(norm(device_source) in norm(device_reads), "Device source missing from actual dependencies")
    archive = backend / "Release/sirius_backend.lib"
    archive_label = archive_member(archive, devices[0])
    canonical_semantics = selected_command_semantics(device_commands)
    for arm in ("baseline", "candidate"):
        records = tlogs["native_wave32_" + arm + "_compute"]["records"]
        for row in records.values():
            require(row["command_semantics"] == canonical_semantics,
                    "Private compute/device command switches differ from canonical Vulkan usage")
    require(tlogs["native_wave32_baseline_compute"]["records"][
        "out/native-wave32-producer-inputs/baseline/vulkan_device.cpp"]["object"]["sha256"] !=
        tlogs["native_wave32_candidate_compute"]["records"][
        "out/native-wave32-producer-inputs/candidate/vulkan_device.cpp"]["object"]["sha256"],
        "Private candidate Vulkan object is identical to baseline")
    return {"schema": "sirius-native-packed-static-PE-audit-v1", "audit_completed": True,
            "array_and_private_compile_checks_passed": True, "full_producer_eligibility_passed": False,
            "native_wave32_device_object_and_map_joins_passed": True,
            "canonical_CXX_source_revision": CANONICAL_CXX_REVISION,
            "zero_changed_SPV_payloads": True, "private_compile_object_count": 7,
            "manifest_identity": identity(export / "manifest.json"),
            "producer_identity": identity(export / "producer-result.json"),
            "source_revision": args.revision, "run_id": args.run_id, "products": products,
            "private_compiles": tlogs, "changed_main_arrays": changed,
            "canonical_device": {"object": identity(devices[0]), "archive": identity(archive),
                                 "exact_archive_member_label": archive_label,
                                 "actual_CL_source_records_joined": True,
                                 "scope": "Canonical archive object retained for provenance only; not selected by either arm MAP."},
            "remaining_manual_gates": [
                "Independent read of actual source/tool/API provenance, complete MSVC commands and VERBOSE:LIB records; this narrow reader is not that whole admission.",
                "Actual selected DLL/CRT/SDK imports, command/source/API/ZIP provenance and full current producer admission; no PE/provider execution is performed here."],
            "runtime_admitted": False, "scope": "Static PE arrays and selected private build records only; not numerical/performance qualification."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # This guard must precede the reporting catch: never write a failure report
    # into the immutable artifact merely because main rejected that destination.
    require(not args.output.resolve().is_relative_to(args.export.resolve()), "Output must be outside export")
    try:
        report = main(args)
    except (ValueError, KeyError, IndexError, TypeError, OSError, struct.error, UnicodeError, ET.ParseError) as error:
        report = {"schema": "sirius-native-packed-static-PE-audit-v1", "audit_completed": False,
                  "array_and_private_compile_checks_passed": False, "full_producer_eligibility_passed": False,
                  "native_wave32_device_object_and_map_joins_passed": False,
                  "canonical_CXX_source_revision": CANONICAL_CXX_REVISION,
                  "failure": str(error), "runtime_admitted": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    raise SystemExit(0 if report["array_and_private_compile_checks_passed"] else 1)
