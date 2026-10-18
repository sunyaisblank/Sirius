"""Static reader for the one matched Windows export; never execute its products.

Header decoding is retained from f851 read_postproof.py::arrays. Windows source/
artifact joins follow the selected-record approach of 9d56406 reconstruction,
without its project imports, staging, subprocesses or native-gate authority.
This source has not been executed. ROOT must review it before use.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import struct
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "out/native-fp32-packed-scratch-producer-inputs"
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
    for option in ("/std:c++latest", "/W4", "/WX", "/permissive-", "/utf-8", "SIRIUS_CONTRACT_MODE=2"):
        require(option.lower() in command.lower(), "Missing strict flag: " + label + " " + option)


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
    before = document(export / "build-inputs-before.json")
    require(before == document(export / "build-inputs-after.json"), "Build seals differ")
    source = document(export / "source-before-configure.json")
    for arm in ("baseline", "candidate"):
        for header in ("retained_kernels.h", "retained_owner_kernels.h", "retained_owner_layout.h"):
            suffix = "out/native-fp32-packed-scratch-producer-inputs/" + arm + "/" + header
            selected_row(source, suffix, identity(INPUTS / arm / header))
            selected_row(before, suffix, identity(INPUTS / arm / header))
    for relative in ("src/sirius/backend/retained_compute.cpp", "tests/backend/retained_dopri_test.cpp",
                     "src/sirius/backend/vulkan/vulkan_device.cpp", "cmake/sirius_toolchain.cmake"):
        selected_row(source, relative, identity(ROOT / relative))
    selected_row(source, "out/native-fp32-packed-scratch-producer-inputs/observer/audited-observer.cpp",
                 identity(INPUTS / "observer/audited-observer.cpp"))
    expected = {arm: arrays(INPUTS / arm / "retained_kernels.h") for arm in ("baseline", "candidate")}
    require(all(len(v) == 40 for v in expected.values()), "Expected40 main arrays")
    changed = [name for name in expected["baseline"] if expected["baseline"][name] != expected["candidate"][name]]
    require(set(changed) == {"kTransportShader", "kEndpointShader", "kDenseShader", "kInitializeShader", "kDopriPhaseShader"},
            "Unexpected arm payload delta")
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
        device_members = set(re.findall(r"sirius_backend(?:\.lib)?[:(][^\s)]*vulkan_device\.obj", mapping, re.I))
        require(len({v.lower() for v in device_members}) == 1, "Canonical device map member missing/ambiguous")
        products[arm] = {"PE": identity(exe), "map": identity(exe.with_suffix(".map")),
                         "main_arrays": verify_arrays(exe, mapping, "retained_program", expected[arm], image_base, sections),
                         "owner_arrays": verify_arrays(exe, mapping, "retained_owner", owner_arrays, image_base, sections),
                         "canonical_device_map_members": sorted(device_members)}
        for name in changed:
            foreign = expected["candidate" if arm == "baseline" else "baseline"][name]
            require(all(data.find(foreign) < 0 for _, _, _, data in sections),
                    "Foreign changed-arm shader remains in PE: " + arm + " " + name)
    tlogs = {}
    targets = ["native_packed_dopri"] + ["native_packed_" + arm + "_" + kind
               for arm in ("baseline", "candidate") for kind in ("compute", "controls")]
    for target in targets:
        base = export / "metadata/private-objects" / (target + ".dir") / "Release"
        objects = list(base.rglob("*.obj"))
        require(len(objects) == 1, "Missing/extra private object: " + target)
        commands = logs(base, "CL.command.")
        reads = logs(base, "CL.read.")
        require(commands and reads, "Missing CL command/read tlogs: " + target)
        command = "\n".join(text(p) for p in commands)
        dependencies = norm("\n".join(text(p) for p in reads))
        check_flags(command, target)
        if target == "native_packed_dopri":
            require("retained_kernels.h" not in dependencies, "Shared DP consumes an arm header")
        else:
            arm = "baseline" if "baseline" in target else "candidate"
            other = "candidate" if arm == "baseline" else "baseline"
            for header in (("retained_kernels.h", "retained_owner_layout.h", "retained_owner_kernels.h")
                           if target.endswith("compute") else ("retained_kernels.h", "retained_owner_layout.h")):
                suffix = "out/native-fp32-packed-scratch-producer-inputs/" + arm + "/" + header
                require(suffix in dependencies, "Private dependency missing: " + suffix)
            require("/" + other + "/retained_kernels.h" not in dependencies, "Other arm header consumed")
        if not target.endswith("compute"):
            require("fail_on_skip.h" in command.lower(), "Missing required-runtime skip refusal")
        project = ET.parse(export / "metadata/projects" / (target + ".vcxproj")).getroot()
        sources = [v.attrib["Include"] for v in project.findall(".//m:ClCompile[@Include]", NS)]
        require(len(sources) == 1, "Private target source count differs")
        expected_source = ("tests/backend/retained_dopri_test.cpp" if target == "native_packed_dopri"
                           else "src/sirius/backend/retained_compute.cpp" if target.endswith("compute")
                           else "out/native-fp32-packed-scratch-producer-inputs/observer/audited-observer.cpp")
        require(norm(sources[0]).endswith(expected_source), "Wrong private source: " + target)
        tlogs[target] = {"object": identity(objects[0]), "source": sources[0],
                         "command_tlogs": [str(p.relative_to(export)) for p in commands],
                         "read_tlogs": [str(p.relative_to(export)) for p in reads]}
    for arm in ("baseline", "candidate"):
        base = export / "metadata/private-objects" / ("native_packed_" + arm + "_controls.dir") / "Release"
        link_reads = logs(base, "link.read.")
        link_commands = logs(base, "link.command.")
        require(link_reads and link_commands, "Missing LINK records: " + arm)
        dependencies = norm("\n".join(text(p) for p in link_reads))
        command = "\n".join(text(p) for p in link_commands)
        require("/wholearchive" not in command.lower(), "Unexpected WHOLEARCHIVE")
        for target in ("native_packed_" + arm + "_compute", "native_packed_" + arm + "_controls", "native_packed_dopri"):
            folder = export / "metadata/private-objects" / (target + ".dir") / "Release"
            obj = next(folder.rglob("*.obj"))
            suffix = norm(obj.relative_to(export / "metadata/private-objects"))
            require(suffix in dependencies, "Private object absent from LINK read: " + suffix)
        require("/src/sirius/backend/release/sirius_backend.lib" in dependencies, "Backend archive not linked")
        products[arm]["link_read_tlogs"] = [str(p.relative_to(export)) for p in link_reads]
        products[arm]["link_command_tlogs"] = [str(p.relative_to(export)) for p in link_commands]
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
    return {"schema": "sirius-native-packed-static-PE-audit-v1", "audit_completed": True,
            "array_and_private_compile_checks_passed": True, "full_producer_eligibility_passed": False,
            "source_revision": args.revision, "run_id": args.run_id, "products": products,
            "private_compiles": tlogs, "changed_main_arrays": changed,
            "canonical_device": {"object": identity(devices[0]), "archive": identity(archive),
                                 "exact_archive_member_label": archive_label,
                                 "actual_CL_source_records_joined": True},
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
                  "failure": str(error), "runtime_admitted": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    raise SystemExit(0 if report["array_and_private_compile_checks_passed"] else 1)
