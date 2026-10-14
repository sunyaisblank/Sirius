#!/usr/bin/env python3
"""Prepare the finite issue-48 source/translation experiment; never build or execute it."""
import difflib
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import stat
import tarfile
import urllib.request

ROOT = Path(__file__).absolute().parent
REPO = ROOT.parents[1]
PRODUCT_SOURCE = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
RUNNER_SOURCE = "5733117ba8ba4430d723d62294f991db423fbf4c"
PINS = [
    ("MoltenVK", "KhronosGroup", "db445ff2042d9ce348c439ad8451112f354b8d2a", 1251145,
     "54e1881e8f10c7dcc4828c14d3e0fc699a6f6c563a78ba8a631c626c836ae429"),
    ("SPIRV-Cross", "KhronosGroup", "adec7acbf41a988713cdb85f93f26c8ca5ea863e", 1952650,
     "1bed78fb2c43c877058c21bf1d9bfa2332df935960c6901bb141f4ef34cd9f03"),
    ("Vulkan-Headers", "KhronosGroup", "6aefb8eb95c8e170d0805fd0f2d02832ec1e099a", 2855412,
     "f7e54b9e95a3da54802d62ce7b4995cbad03ae2e4325c26377ea4ec04d0429f1"),
    ("SPIRV-Headers", "KhronosGroup", "b824a462d4256d720bebb40e78b9eb8f78bbb305", 558965,
     "c693867f10a7760ef1bcf85419d51783586768cc2c601d03841bc6a8b2554b9c"),
    ("SPIRV-Tools", "KhronosGroup", "262bdab48146c937467f826699a40da0fdfc0f1a", 3428734,
     "8f830936b057003be77e90d37068ff6066eb87a1bab0033f8f2f5feaa966ce56"),
    ("cereal", "USCiLab", "a56bad8bbb770ee266e930c95d37fff2a5be7fea", 386830,
     "12e10fb541270632f6683286ce8c7e8c321fe4addb53fa91ed498b95c1758f93"),
]
CPM_URL = "https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake"
CPM_PIN = (42313, "78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791")
CAPTURE = Path("attestations/native-vulkan/4bcd987/macos-retained-diagnostic-a1d6e75/files/out/macos-retained-diagnostic/observation/metal-source")
CAPTURE_PINS = {
    "initialize.spv": (77444, "d22a8770b3f52a49d6075cdc9d1884115eb0c69b76f72b6cba8b1b793ad9861c"),
    "initialize.metal": (85091, "b2d7daec9928356de4ca24a724d135adc59685a330b043e4fa55f2fe2f9b5550"),
}
DERIVED_PIN = (82487, "80b33059934cc77a685c5738b9413c94379e45f1b18d0e6635a6432ff0407921")
FORCED = b"static inline __attribute__((always_inline))"


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def encoded(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def write_json(path, value):
    path.write_bytes(encoded(value))


def pin_bytes(data):
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def expect(data, pin, label):
    require(pin_bytes(data) == {"bytes": pin[0], "sha256": pin[1]}, "different bytes: " + label)


def fetch(url, destination, expected):
    require(not destination.exists() and not destination.is_symlink(), "existing download: " + str(destination))
    partial = destination.with_name(destination.name + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "Sirius-issue48-source-preparation"})
    digest, count = hashlib.sha256(), 0
    with urllib.request.urlopen(request, timeout=45) as response, partial.open("xb") as output:
        require(response.geturl().startswith("https://"), "non-HTTPS redirect")
        while chunk := response.read(1024 * 1024):
            count += len(chunk)
            require(count <= expected[0], "download exceeds exact byte count")
            digest.update(chunk)
            output.write(chunk)
    require((count, digest.hexdigest()) == expected, "archive/script byte seal failed: " + url)
    partial.rename(destination)


def extract(archive, destination, prefix):
    """No traversal, special nodes, hardlinks, link ancestors or escaping symlinks."""
    require(not destination.exists(), "existing source destination")
    with tarfile.open(archive, "r:gz") as bundle:
        members = bundle.getmembers()
        require(len(members) < 50000 and sum(m.size for m in members) < 256 * 1024**2,
                "source archive exceeds preparation bounds")
        plan, links = {}, {}
        for member in members:
            raw = member.name.rstrip("/")
            require("\\" not in raw and "\0" not in raw and not raw.startswith("/"), "unsafe archive name")
            parts = PurePosixPath(raw).parts
            require(parts and parts[0] == prefix and all(p not in ("..", ".") for p in parts),
                    "archive prefix/traversal differs")
            if len(parts) == 1:
                require(member.isdir(), "archive root is not a directory")
                continue
            name = "/".join(parts[1:])
            require(name not in plan and not member.mode & 0o7000, "duplicate or privileged archive node")
            require(member.isfile() or member.isdir() or member.issym(), "unsupported archive node")
            require(member.size < 32 * 1024**2, "oversized archive member")
            plan[name] = member
            if member.issym():
                target = member.linkname
                require(target and "\\" not in target and "\0" not in target and not target.startswith("/"), "unsafe symlink")
                joined = posixpath.normpath(posixpath.join(posixpath.dirname(name), target))
                require(joined != ".." and not joined.startswith("../"), "symlink escapes source root")
                links[name] = target
        for name in plan:
            require(not any(str(parent) in links for parent in PurePosixPath(name).parents),
                    "archive member descends through a symlink")
        destination.mkdir()
        for name, member in sorted(plan.items(), key=lambda item: (len(PurePosixPath(item[0]).parts), item[0])):
            path = destination / name
            if member.isdir():
                path.mkdir(exist_ok=True)
            elif member.isfile():
                path.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(member)
                require(source is not None, "unreadable archive member")
                data = source.read()
                require(len(data) == member.size, "short archive member")
                with path.open("xb") as output:
                    output.write(data)
                path.chmod(member.mode & 0o777)
        for name, target in links.items():
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.symlink_to(target)
        for name in links:
            require((destination / name).resolve().is_relative_to(destination.resolve()), "resolved source symlink escapes")
        for name, member in sorted(plan.items(), reverse=True):
            if member.isdir():
                (destination / name).chmod(member.mode & 0o777)
    return {"members": len(plan), "symlinks": len(links)}


def replace_once(text, original, replacement, label):
    require(text.count(original) == 1, "source guard differs: " + label)
    return text.replace(original, replacement, 1)


def source_patch(path, transform, relative):
    original = path.read_text()
    changed = transform(original)
    require(changed != original, "empty source patch")
    path.write_text(changed)
    return "".join(difflib.unified_diff(original.splitlines(True), changed.splitlines(True),
                                      fromfile="a/" + relative, tofile="b/" + relative))


def metadata_cmake(original):
    begin = "# Record the MoltenVK GIT revision as a derived header file suitable for including in a build\n"
    end = '\n\ntarget_compile_definitions(MoltenVK PUBLIC\n\t"MVK_FRAMEWORK_VERSION=${PROJECT_VERSION}"\n)'
    require(original.count(begin) == original.count(end) == 1, "CMake revision block boundary differs")
    start, finish = original.index(begin), original.index(end)
    old = original[start:finish]
    require('COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD' in old and 'OUTPUT_VARIABLE MVK_GIT_REV' in old,
            "CMake revision block body differs")
    replacement = '''# Sirius issue-48 diagnostic: this is a source/configuration content discriminator.
# It is NOT an upstream or local Git revision, despite the legacy template name.
if(NOT DEFINED MVK_DIAGNOSTIC_CONTENT_ID)
    message(FATAL_ERROR "The diagnostic provider requires its declared content discriminator")
endif()
string(LENGTH "${MVK_DIAGNOSTIC_CONTENT_ID}" MVK_DIAGNOSTIC_CONTENT_ID_LENGTH)
if(NOT MVK_DIAGNOSTIC_CONTENT_ID_LENGTH EQUAL 64 OR
   NOT MVK_DIAGNOSTIC_CONTENT_ID MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR "The diagnostic content discriminator must be 64 lowercase hex digits")
endif()
string(SUBSTRING "${MVK_DIAGNOSTIC_CONTENT_ID}" 0 8 MVK_DIAGNOSTIC_CACHE_PREFIX)
if(MVK_DIAGNOSTIC_CACHE_PREFIX STREQUAL "00000000" OR
   MVK_DIAGNOSTIC_CACHE_PREFIX STREQUAL "db445ff2")
    message(FATAL_ERROR "The diagnostic cache prefix is zero or aliases the upstream base")
endif()
set(MVK_GIT_REV "${MVK_DIAGNOSTIC_CONTENT_ID}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/../Templates/cmake/mvkGitRevDerived.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/mvkGitRevDerived.h"
)
target_compile_definitions(MoltenVK PRIVATE
    "MVK_DIAGNOSTIC_CONTENT_ID=\\"${MVK_DIAGNOSTIC_CONTENT_ID}\\""
)
'''
    return original[:start] + replacement + original[finish:]


def metadata_device(original):
    result = replace_once(original,
        '\tstrcpy(supportedProps12.driverInfo, MVK_VERSION_STRING);',
        '\tstrcpy(supportedProps12.driverInfo, MVK_VERSION_STRING "+sirius-issue48-unforced." MVK_DIAGNOSTIC_CONTENT_ID);',
        "driverInfo")
    result = replace_once(result,
        '\t// First 4 bytes contains MoltenVK revision.\n\t// This is captured either as the MoltenVK Git revision, or if that\'s not available, as the MoltenVK version.',
        '\t// Diagnostic: first four bytes contain the declared source/config content token.\n\t// This discriminator is not a Git revision; full provenance is sealed separately.', "UUID description")
    begin = '// Retrieve the SPIRV-Cross Git revision hash from a derived header file,\n'
    end = 'uint32_t MVKPhysicalDevice::getMoltenVKGitRevision() {'
    require(result.count(begin) == result.count(end) == 1, "revision reader comment differs")
    start, finish = result.index(begin), result.index(end)
    return result[:start] + ('// Read the diagnostic content token from the legacy generated-header field.\n'
                              '// The legacy function name does not confer Git provenance.\n') + result[finish:]


def main():
    require(ROOT == REPO / "out/macos-moltenvk-unforced" and ROOT.resolve() == ROOT,
            "preparation must use its exact non-symlink owned root")
    for folder in ("archives", "sources", "captured", "patches", "cpm-cache", "observations"):
        require(not (ROOT / folder).exists(), "existing preparation products; preserve failure before a new attempt")
        (ROOT / folder).mkdir()
    (ROOT / "cpm-cache/cpm").mkdir()
    extraction, archive_pins = {}, {}
    for name, org, revision, count, digest in PINS:
        url = f"https://codeload.github.com/{org}/{name}/tar.gz/{revision}"
        archive = ROOT / "archives" / (name + ".tar.gz")
        fetch(url, archive, (count, digest))
        archive_pins[name] = {"url": url, "commit": revision, "bytes": count, "sha256": digest}
        destinations = ["cross-stock", "cross-patched"] if name == "SPIRV-Cross" else ["moltenvk-patched" if name == "MoltenVK" else name]
        for destination in destinations:
            extraction[destination] = extract(archive, ROOT / "sources" / destination, name + "-" + revision)
    fetch(CPM_URL, ROOT / "cpm-cache/cpm/CPM_0.40.8.cmake", CPM_PIN)
    revisions = ROOT / "sources/moltenvk-patched/ExternalRevisions"
    for name, _, revision, _, _ in PINS[1:]:
        require((revisions / (name + "_repo_revision")).read_text().strip() == revision,
                "MoltenVK dependency revision differs: " + name)

    patch = source_patch(ROOT / "sources/cross-patched/spirv_msl.cpp",
        lambda text: replace_once(text,
            '\tif (!processing_entry_point)\n\t\tstatement(force_inline);',
            '\tif (!processing_entry_point)\n\t\tstatement("static inline");', "non-entry prototype emitter"),
        "spirv_msl.cpp")
    (ROOT / "patches/SPIRV-Cross-unforced.patch").write_text(patch)
    metadata = source_patch(ROOT / "sources/moltenvk-patched/MoltenVK/CMakeLists.txt", metadata_cmake,
                            "MoltenVK/CMakeLists.txt")
    metadata += source_patch(ROOT / "sources/moltenvk-patched/MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm",
                             metadata_device, "MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm")
    (ROOT / "patches/MoltenVK-diagnostic-identity.patch").write_text(metadata)

    captured_origins = {}
    for name, pin in CAPTURE_PINS.items():
        origin = REPO / CAPTURE / name
        data = origin.read_bytes()
        expect(data, pin, str(origin))
        (ROOT / "captured" / name).write_bytes(data)
        captured_origins[name] = {"origin": str(CAPTURE / name), **pin_bytes(data)}
    original = (ROOT / "captured/initialize.metal").read_bytes()
    require(original.count(FORCED) == original.count(b"__attribute__((always_inline))") == 84,
            "original annotation inventory differs")
    derived = original.replace(FORCED, b"static inline")
    expect(derived, DERIVED_PIN, "annotation-only control")
    (ROOT / "captured/initialize-without-forced-inline.metal").write_bytes(derived)

    profile = {
        "schema": 1, "native_target": "MoltenVK", "native_arch": "arm64", "configuration": "Release",
        "deployment_target": "11.0", "public_metal_api_only": True, "shader_converter_tool": False,
        "cereal_enabled": True, "spirv_tools_enabled": True, "ccache_enabled": False,
        "build_parallelism": 1, "library_optimization": "unchanged MoltenVK Default",
        "native_toolchain": "Root must separately seal the selected Xcode, macOS SDK, clang, CMake and Ninja before build",
        "source_overrides": [name for name, *_ in PINS[1:]],
        "diagnostic_provider": True, "upstream_git_attestation": False,
    }
    write_json(ROOT / "build-profile.json", profile)
    content_inputs = {
        "archives": archive_pins,
        "patches": {p.name: pin_bytes(p.read_bytes()) for p in sorted((ROOT / "patches").iterdir())},
        "build_profile": profile,
        "meaning": "Source/configuration cache discriminator, not a Git commit or native-build identity. The ID itself is excluded from its preimage.",
    }
    content_id = hashlib.sha256(encoded(content_inputs)).hexdigest()
    require(content_id[:8] not in ("00000000", "db445ff2"), "invalid/aliased cache discriminator")
    write_json(ROOT / "content-id.json", {"content_id": content_id, "preimage": content_inputs,
        "driver_info_expected": "1.4.1+sirius-issue48-unforced." + content_id,
        "cache_uuid_expected_first_four_bytes_hex": content_id[:8],
        "required_native_check": "Query full stock and diagnostic UUIDs on the same host; require distinct values and the declared diagnostic prefix before original numerical controls."})
    write_json(ROOT / "preparation.json", {
        "status": "prepared_without_build_or_execution", "product_source": PRODUCT_SOURCE,
        "prepared_for_runner": RUNNER_SOURCE, "root": str(ROOT), "archives": archive_pins,
        "cpm": {"url": CPM_URL, "bytes": CPM_PIN[0], "sha256": CPM_PIN[1]},
        "extraction": extraction, "captured_origins": captured_origins, "content_id": content_id,
        "translation_configuration": "Source-backed reconstruction, not original full option readback; exact stock MSL equality is mandatory before candidate acceptance.",
        "native_build_products": "bin/macos/macos-moltenvk-unforced",
        "linux_build_products": "bin/linux-gcc/macos-moltenvk-unforced",
        "compiler_executed": False, "provider_executed": False, "qualification_claimed": False,
    })
    files, links = {}, {}
    for base, dirs, names in os.walk(ROOT, followlinks=False):
        directory = Path(base)
        for name in [*dirs, *names]:
            path = directory / name
            relative = str(path.relative_to(ROOT))
            if path.is_symlink():
                require(path.resolve().is_relative_to(ROOT), "source link escapes preparation")
                links[relative] = os.readlink(path)
            elif path.is_file():
                files[relative] = {**pin_bytes(path.read_bytes()), "mode": stat.S_IMODE(path.stat().st_mode)}
    write_json(ROOT / "input-freeze.json", {"schema": 1, "files": files, "symlinks": links,
        "scope": "Prepared inputs only; future build products and observations are separate outputs."})
    print(json.dumps({"status": "prepared_without_execution", "files": len(files), "symlinks": len(links),
        "content_id": content_id, "input_freeze": pin_bytes((ROOT / "input-freeze.json").read_bytes())}))


if __name__ == "__main__":
    main()
