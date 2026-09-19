"""First-party ownership, dependency layers and repository-policy orchestration."""

from __future__ import annotations

import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

from .build_policy import (
    GIT_ATTRIBUTES,
    STRICT_TEST_VOLUME_TARGETS,
    attestation_preflight_errors,
    attestation_source_authority_errors,
    authority_checkout_errors,
    immutable_input_errors,
    integration_boundary_errors,
    strict_test_volume_errors,
    verify_attestation_source_authority_policy,
    verify_authority_checkout_policy,
    verify_immutable_input_policy,
    verify_integration_boundary_policy,
    verify_strict_test_volume_policy,
)
from .colour_policy import (
    LIVE_VIEWER_SHADERS,
    aces_contract_errors,
    blackbody_laws_authority_errors,
    cie1931_observer_authority_errors,
    srgb_transfer_authority_errors,
    verify_aces_contract_policy,
    verify_blackbody_laws_authority_policy,
    verify_cie1931_observer_authority_policy,
    verify_srgb_transfer_authority_policy,
    verify_xyz_srgb_authority_policy,
    xyz_srgb_authority_errors,
)
from .common import (
    CPP_NON_CODE,
    OPERATING_MODEL,
    ROOT,
    SOURCE_ROOT,
    relative,
)
from .metric_policy import (
    KERR_SCHILD_CURVATURE_ORACLE,
    KERR_SCHILD_FIELD_EQUATION_GATE,
    KERR_SCHILD_OPERATING_MODEL,
    KERR_SCHILD_TEST_CMAKE,
    alcubierre_authority_errors,
    kerr_schild_field_equation_errors,
    morris_thorne_authority_errors,
    verify_alcubierre_authority_policy,
    verify_kerr_schild_field_equation_policy,
    verify_morris_thorne_authority_policy,
)
from .transport_policy import (
    kerr_zamo_transfer_authority_errors,
    page_thorne_edge_authority_errors,
    thin_lens_authority_errors,
    verify_kerr_zamo_transfer_authority_policy,
    verify_page_thorne_edge_authority_policy,
    verify_thin_lens_authority_policy,
    verify_volumetric_transfer_authority_policy,
    volumetric_transfer_authority_errors,
)


LAYER_DEPENDENCIES = {
    "base": {"base"},
    "core": {"base", "core"},
    "oracle": {"base", "oracle"},
    "backend": {"base", "core", "backend"},
    "render": {"base", "core", "backend", "render"},
    "app": {"base", "core", "backend", "render", "app"},
}
EXPECTED_VENDOR_DIRECTORIES = {"glad", "glfw", "stb", "tinyexr"}
RETIRED_PATHS = (
    ROOT / "lib" / "glm",
    ROOT / "lib" / "imgui",
    ROOT / "lib" / "glfw" / "examples",
    ROOT / "lib" / "glfw" / "tests",
    ROOT / "lib" / "glfw" / "deps" / "glad",
    ROOT / "lib" / "glfw" / "deps" / "getopt.c",
    ROOT / "lib" / "glfw" / "deps" / "getopt.h",
    ROOT / "lib" / "glfw" / "deps" / "linmath.h",
    ROOT / "lib" / "glfw" / "deps" / "nuklear.h",
    ROOT / "lib" / "glfw" / "deps" / "nuklear_glfw_gl2.h",
    ROOT / "lib" / "glfw" / "deps" / "stb_image_write.h",
    ROOT / "lib" / "glfw" / "deps" / "tinycthread.c",
    ROOT / "lib" / "glfw" / "deps" / "tinycthread.h",
    SOURCE_ROOT / "render" / "render_config.h",
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD004A.frag",
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD004A.vert",
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD005A.frag",
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD005A.vert",
)
CPP_TOKEN = re.compile(r'(?<![A-Za-z0-9_${])([A-Za-z0-9_./-]+\.cpp)')
PROJECT_INCLUDE = re.compile(r'^\s*#\s*include\s+"sirius/([^/]+)/[^\"]+"', re.MULTILINE)
CAMEL_IDENTIFIER = re.compile(r"\b[a-z][a-z0-9_]*[A-Z][A-Za-z0-9_]*\b")
BOUNDARY_HEADERS = (
    SOURCE_ROOT / "app" / "config" / "config_schema.h",
    SOURCE_ROOT / "app" / "config" / "session_config_adapter.h",
    SOURCE_ROOT / "app" / "viewer" / "interactive_viewer.h",
    SOURCE_ROOT / "backend" / "cpu" / "geodesic_tracer.h",
    SOURCE_ROOT / "render" / "film_config.h",
    SOURCE_ROOT / "render" / "exr_writer.h",
    SOURCE_ROOT / "render" / "session" / "render_session.h",
)


def cmake_source_owners() -> Counter[Path]:
    owners: Counter[Path] = Counter()
    for cmake_file in SOURCE_ROOT.rglob("CMakeLists.txt"):
        cmake = re.sub(r"#[^\n]*", "", cmake_file.read_text(encoding="utf-8"))
        cmake = re.sub(r"set_source_files_properties\s*\([^)]*\)", "", cmake)
        for token in CPP_TOKEN.findall(cmake):
            candidate = (cmake_file.parent / token).resolve()
            if candidate.is_file() and candidate.is_relative_to(SOURCE_ROOT):
                owners[candidate] += 1
    return owners


def verify() -> list[str]:
    errors: list[str] = []

    verify_immutable_input_policy()
    verify_strict_test_volume_policy()
    verify_integration_boundary_policy()
    verify_authority_checkout_policy()
    verify_attestation_source_authority_policy()
    verify_srgb_transfer_authority_policy()
    verify_xyz_srgb_authority_policy()
    verify_cie1931_observer_authority_policy()
    verify_aces_contract_policy()
    verify_blackbody_laws_authority_policy()
    verify_page_thorne_edge_authority_policy()
    verify_thin_lens_authority_policy()
    verify_kerr_zamo_transfer_authority_policy()
    verify_volumetric_transfer_authority_policy()
    verify_morris_thorne_authority_policy()
    verify_alcubierre_authority_policy()
    verify_kerr_schild_field_equation_policy()
    try:
        attribute_source = GIT_ATTRIBUTES.read_text(encoding="utf-8")
        attributes = subprocess.run(
            ["git", "check-attr", "text", "eol", "--", relative(OPERATING_MODEL)],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError, UnicodeDecodeError) as error:
        errors.append(f"could not inspect operating-model checkout attributes: {error}")
    else:
        errors.extend(
            authority_checkout_errors(
                OPERATING_MODEL.read_bytes(), attributes, attribute_source
            )
        )
    errors.extend(
        attestation_source_authority_errors(
            (ROOT / "scripts" / "verify-attestation.py").read_text(encoding="utf-8")
        )
    )
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    dependencies = (ROOT / "cmake" / "sirius_dependencies.cmake").read_text(
        encoding="utf-8"
    )
    errors.extend(immutable_input_errors(workflow, dependencies))
    errors.extend(integration_boundary_errors(workflow))
    strict_test_documents = {
        relative(path): (target, path.read_text(encoding="utf-8"))
        for path, target in STRICT_TEST_VOLUME_TARGETS.items()
    }
    errors.extend(strict_test_volume_errors(strict_test_documents))
    governed_sources = {
        path: path.read_text(encoding="utf-8")
        for path in SOURCE_ROOT.rglob("*")
        if path.suffix in {".h", ".cpp", ".slang", ".frag", ".vert"}
    }
    errors.extend(srgb_transfer_authority_errors(governed_sources))
    errors.extend(xyz_srgb_authority_errors(governed_sources))
    errors.extend(cie1931_observer_authority_errors(governed_sources))
    errors.extend(aces_contract_errors(governed_sources))
    errors.extend(blackbody_laws_authority_errors(governed_sources))
    errors.extend(page_thorne_edge_authority_errors(governed_sources))
    errors.extend(thin_lens_authority_errors(governed_sources))
    errors.extend(kerr_zamo_transfer_authority_errors(governed_sources))
    errors.extend(volumetric_transfer_authority_errors(governed_sources))
    errors.extend(morris_thorne_authority_errors(governed_sources))
    errors.extend(alcubierre_authority_errors(governed_sources))
    kerr_schild_documents = dict(governed_sources)
    for path in (
        KERR_SCHILD_CURVATURE_ORACLE,
        KERR_SCHILD_FIELD_EQUATION_GATE,
        KERR_SCHILD_TEST_CMAKE,
        KERR_SCHILD_OPERATING_MODEL,
    ):
        if path.is_file():
            kerr_schild_documents[path] = path.read_text(encoding="utf-8")
    errors.extend(kerr_schild_field_equation_errors(kerr_schild_documents))
    shader_root = SOURCE_ROOT / "app" / "viewer" / "shaders"
    actual_viewer_shaders = {
        path for path in shader_root.iterdir() if path.suffix in {".frag", ".vert"}
    }
    if actual_viewer_shaders != LIVE_VIEWER_SHADERS:
        errors.append(
            "live viewer shader inventory differs from the exact RDSD003A transfer pair"
        )
    preflight_path = ROOT / "scripts" / "attestation_preflight.py"
    reuse_path = ROOT / "scripts" / "reuse_qualification_evidence.py"
    if not preflight_path.is_file():
        errors.append("attestation preflight producer is missing")
    elif not reuse_path.is_file():
        errors.append("qualification evidence reuse producer is missing")
    else:
        errors.extend(attestation_preflight_errors(
            preflight_path.read_text(encoding="utf-8"),
            reuse_path.read_text(encoding="utf-8"),
            (ROOT / "scripts" / "validate-native-runtime.py").read_text(
                encoding="utf-8"
            ),
            (ROOT / "scripts" / "validate-viewer.sh").read_text(encoding="utf-8"),
            (ROOT / "tests" / "CMakeLists.txt").read_text(encoding="utf-8"),
            (ROOT / "docs" / "ATTESTATION.md").read_text(encoding="utf-8"),
        ))

    for path in RETIRED_PATHS:
        if path.exists():
            errors.append(f"retired repository content returned: {relative(path)}")

    vendor_root = ROOT / "lib"
    actual_vendors = {path.name for path in vendor_root.iterdir() if path.is_dir()}
    unexpected_vendors = sorted(actual_vendors - EXPECTED_VENDOR_DIRECTORIES)
    missing_vendors = sorted(EXPECTED_VENDOR_DIRECTORIES - actual_vendors)
    if unexpected_vendors:
        errors.append("unapproved vendored directories: " + ", ".join(unexpected_vendors))
    if missing_vendors:
        errors.append("required vendored directories are missing: " + ", ".join(missing_vendors))

    owners = cmake_source_owners()
    for source in SOURCE_ROOT.rglob("*.cpp"):
        count = owners[source.resolve()]
        if count != 1:
            errors.append(
                f"{relative(source)} must have exactly one CMake owner; found {count}"
            )

    backend_cmake = (SOURCE_ROOT / "backend" / "CMakeLists.txt").read_text(encoding="utf-8")
    render_cmake = (SOURCE_ROOT / "render" / "CMakeLists.txt").read_text(encoding="utf-8")
    if "add_library(sirius_backend_cpu STATIC" not in backend_cmake:
        errors.append("the CPU tracer is not owned by the sirius_backend_cpu target")
    if "cpu/geodesic_tracer.cpp" not in backend_cmake:
        errors.append("sirius_backend_cpu does not compile cpu/geodesic_tracer.cpp")
    if "geodesic_tracer.cpp" in render_cmake:
        errors.append("the render target directly compiles a backend implementation")
    if "sirius_backend_cpu" not in render_cmake:
        errors.append("the render target has no explicit CPU-backend dependency")

    for layer, allowed in LAYER_DEPENDENCIES.items():
        for source in (SOURCE_ROOT / layer).rglob("*"):
            if source.suffix not in {".h", ".cpp"}:
                continue
            text = source.read_text(encoding="utf-8")
            for dependency in PROJECT_INCLUDE.findall(text):
                if dependency not in allowed:
                    errors.append(
                        f"{relative(source)} reaches upward or sideways into sirius::{dependency}"
                    )

    for header in BOUNDARY_HEADERS:
        text = CPP_NON_CODE.sub(" ", header.read_text(encoding="utf-8"))
        for identifier in CAMEL_IDENTIFIER.findall(text):
            if identifier.startswith("k") and len(identifier) > 1 and identifier[1].isupper():
                continue
            errors.append(
                f"{relative(header)} exposes non-snake-case identifier '{identifier}'"
            )

    return errors


def main() -> int:
    errors = verify()
    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        return 1
    source_count = sum(1 for _ in SOURCE_ROOT.rglob("*.cpp"))
    print(
        f"repository structure owns {source_count} first-party translation units exactly once, "
        f"enforces {len(LAYER_DEPENDENCIES)} dependency layers, and permits "
        f"{len(EXPECTED_VENDOR_DIRECTORIES)} live vendor directories with immutable build inputs "
        "and strict test-volume identity"
    )
    return 0
