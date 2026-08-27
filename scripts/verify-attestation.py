#!/usr/bin/env python3
"""Validate external-domain attestations without converting them into local evidence."""

import argparse
import hashlib
import json
import math
import re
import shutil
import struct
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zlib
from datetime import datetime, timezone
from pathlib import Path


KNOWN_DOMAINS = {
    "physical-radeon-780m",
    "wsl2-dozen",
    "windows-native-build",
    "windows-native-vulkan",
    "macos-native-build",
    "macos-moltenvk",
    "physical-imax-5616x4096",
    "viewer-native-window-input",
}
SCENE_EVIDENCE_PREFIX = "[Session] Scene evidence: "
MINIMUM_SOURCE_AVAILABLE_TESTS = 700
QUALIFICATION_TESTED_ARTIFACTS = {
    "sirius", "sirius_app_tests", "sirius_backend_tests", "sirius_base_tests",
    "sirius_core_tests", "sirius_oracle_tests", "sirius_render_tests",
}
QUALIFICATION_PRODUCT_ARTIFACTS = {
    "alignment_receipt", "operating_model", "sirius", "starfield",
    "trace_fp32comp_spv", "trace_fp64_spv", "trace_spv",
    "viewer_rdsd003a_fragment", "viewer_rdsd003a_vertex",
    "viewer_rdsd004a_fragment", "viewer_rdsd004a_vertex",
    "viewer_rdsd005a_fragment", "viewer_rdsd005a_vertex",
}
QUALIFICATION_PRODUCT_EVIDENCE = {
    "operating_model": "qualification-product-operating_model",
    "starfield": "qualification-product-starfield",
    "trace_fp32comp_spv": "qualification-product-trace_fp32comp_spv",
    "trace_fp64_spv": "qualification-product-trace_fp64_spv",
    "trace_spv": "qualification-product-trace_spv",
    "viewer_rdsd003a_fragment": "qualification-product-viewer_rdsd003a_fragment",
    "viewer_rdsd003a_vertex": "qualification-product-viewer_rdsd003a_vertex",
    "viewer_rdsd004a_fragment": "qualification-product-viewer_rdsd004a_fragment",
    "viewer_rdsd004a_vertex": "qualification-product-viewer_rdsd004a_vertex",
    "viewer_rdsd005a_fragment": "qualification-product-viewer_rdsd005a_fragment",
    "viewer_rdsd005a_vertex": "qualification-product-viewer_rdsd005a_vertex",
}
QUALIFICATION_RUNTIME_RESOURCE_PATHS = {
    "operating_model": Path("model/operating_model.json"),
    "starfield": Path("assets/Starfield.png"),
    "trace_fp32comp_spv": Path("kernels/trace_fp32comp.spv"),
    "trace_fp64_spv": Path("kernels/trace_fp64.spv"),
    "trace_spv": Path("kernels/trace.spv"),
    "viewer_rdsd003a_fragment": Path("shaders/RDSD003A.frag"),
    "viewer_rdsd003a_vertex": Path("shaders/RDSD003A.vert"),
    "viewer_rdsd004a_fragment": Path("shaders/RDSD004A.frag"),
    "viewer_rdsd004a_vertex": Path("shaders/RDSD004A.vert"),
    "viewer_rdsd005a_fragment": Path("shaders/RDSD005A.frag"),
    "viewer_rdsd005a_vertex": Path("shaders/RDSD005A.vert"),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def text_contains(device, *needles):
    text = " ".join(
        str(device.get(field, ""))
        for field in ("name", "driver_name", "driver_info")
    ).casefold()
    return any(needle.casefold() in text for needle in needles)


def require_positive_integer(value, field):
    require(type(value) is int and value > 0,
            f"runtime device {field} must be a positive integer")


def paeth_predictor(left, above, upper_left):
    estimate = left + above - upper_left
    left_distance = abs(estimate - left)
    above_distance = abs(estimate - above)
    upper_left_distance = abs(estimate - upper_left)
    if left_distance <= above_distance and left_distance <= upper_left_distance:
        return left
    if above_distance <= upper_left_distance:
        return above
    return upper_left


def inspect_png(path):
    """Decode an 8-bit non-interlaced PNG and measure represented RGB structure."""
    payload = path.read_bytes()
    require(payload.startswith(b"\x89PNG\r\n\x1a\n"),
            f"artifact is not a PNG: {path.name}")
    cursor = 8
    header = None
    compressed = bytearray()
    saw_end = False
    while cursor < len(payload):
        require(cursor + 12 <= len(payload), f"truncated PNG chunk: {path.name}")
        length = struct.unpack(">I", payload[cursor:cursor + 4])[0]
        chunk_type = payload[cursor + 4:cursor + 8]
        chunk_end = cursor + 12 + length
        require(chunk_end <= len(payload), f"truncated PNG payload: {path.name}")
        chunk = payload[cursor + 8:cursor + 8 + length]
        recorded_crc = struct.unpack(">I", payload[cursor + 8 + length:chunk_end])[0]
        require(zlib.crc32(chunk_type + chunk) & 0xffffffff == recorded_crc,
                f"PNG chunk checksum mismatch: {path.name}")
        if chunk_type == b"IHDR":
            require(header is None and length == 13, f"invalid PNG IHDR: {path.name}")
            header = struct.unpack(">IIBBBBB", chunk)
        elif chunk_type == b"IDAT":
            compressed.extend(chunk)
        elif chunk_type == b"IEND":
            require(length == 0, f"invalid PNG IEND: {path.name}")
            saw_end = True
            break
        cursor = chunk_end

    require(header is not None and compressed and saw_end,
            f"PNG is missing IHDR, IDAT, or IEND: {path.name}")
    width, height, bit_depth, color_type, compression, filtering, interlace = header
    require(0 < width <= 8192 and 0 < height <= 8192,
            f"PNG dimensions are outside the represented render domain: {path.name}")
    require(bit_depth == 8 and color_type in {0, 2, 6}
            and compression == 0 and filtering == 0 and interlace == 0,
            f"PNG must be non-interlaced 8-bit grayscale, RGB, or RGBA: {path.name}")
    channels = {0: 1, 2: 3, 6: 4}[color_type]
    row_bytes = width * channels
    expected_bytes = height * (row_bytes + 1)
    try:
        decompressor = zlib.decompressobj()
        filtered = decompressor.decompress(bytes(compressed), expected_bytes + 1)
    except zlib.error as error:
        raise ValueError(f"PNG pixel stream is invalid: {path.name}: {error}") from error
    require(decompressor.eof and not decompressor.unused_data
            and not decompressor.unconsumed_tail and len(filtered) == expected_bytes,
            f"PNG pixel stream has the wrong size: {path.name}")

    previous = bytearray(row_bytes)
    minima = [255, 255, 255]
    maxima = [0, 0, 0]
    colours = set()
    dark_channels = 0
    bright_channels = 0
    bright_rows = 0
    row_signatures = set()
    dark_table = bytes(1 if value <= 64 else 0 for value in range(256))
    bright_table = bytes(1 if value >= 96 else 0 for value in range(256))
    for y in range(height):
        offset = y * (row_bytes + 1)
        filter_type = filtered[offset]
        require(filter_type <= 4, f"PNG uses an invalid row filter: {path.name}")
        row = bytearray(filtered[offset + 1:offset + 1 + row_bytes])
        if filter_type != 0:
            for index in range(row_bytes):
                left = row[index - channels] if index >= channels else 0
                above = previous[index]
                upper_left = previous[index - channels] if index >= channels else 0
                if filter_type == 1:
                    prediction = left
                elif filter_type == 2:
                    prediction = above
                elif filter_type == 3:
                    prediction = (left + above) // 2
                else:
                    prediction = paeth_predictor(left, above, upper_left)
                row[index] = (row[index] + prediction) & 0xff

        if channels == 1:
            values = row
            value_min = min(values)
            value_max = max(values)
            minima = [min(current, value_min) for current in minima]
            maxima = [max(current, value_max) for current in maxima]
            if len(colours) < 64:
                colours.update((value, value, value) for value in values)
            row_dark = values.translate(dark_table).count(1)
            row_bright = values.translate(bright_table).count(1)
            dark_channels += row_dark * 3
            bright_channels += row_bright * 3
        else:
            row_bright = 0
            for channel in range(3):
                values = row[channel::channels]
                minima[channel] = min(minima[channel], min(values))
                maxima[channel] = max(maxima[channel], max(values))
                dark_channels += values.translate(dark_table).count(1)
                channel_bright = values.translate(bright_table).count(1)
                bright_channels += channel_bright
                row_bright += channel_bright
            if len(colours) < 64:
                colours.update(zip(row[0::channels], row[1::channels], row[2::channels]))
        if row_bright > 0:
            bright_rows += 1
        if len(row_signatures) < 64:
            row_signatures.add(zlib.crc32(row))
        previous = row

    channel_samples = width * height * 3
    return {
        "dimensions": (width, height),
        "channel_ranges": tuple(high - low for low, high in zip(minima, maxima)),
        "distinct_colours": len(colours),
        "dark_channel_fraction": dark_channels / channel_samples,
        "bright_channel_fraction": bright_channels / channel_samples,
        "bright_row_fraction": bright_rows / height,
        "distinct_row_signatures": len(row_signatures),
    }


def require_matching_number(actual, expected, field, tolerance=1.0e-6):
    require(type(actual) in (int, float) and math.isfinite(actual),
            f"scene transcript {field} must be finite")
    require(type(expected) in (int, float) and math.isfinite(expected),
            f"attested scene {field} must be finite")
    require(math.isclose(float(actual), float(expected), rel_tol=tolerance,
                         abs_tol=tolerance),
            f"scene transcript {field} does not match the attested claim")


def verify_governed_scene_transcript(
        transcripts, scene, render, png_name, device, dimensions, label):
    """Bind one governed P3/P5 claim to a completed typed-session event."""
    candidates = []
    for path in transcripts:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        evidence_indices = [
            index for index, line in enumerate(lines)
            if line.startswith(SCENE_EVIDENCE_PREFIX)
        ]
        for position, index in enumerate(evidence_indices):
            payload = lines[index][len(SCENE_EVIDENCE_PREFIX):]
            try:
                evidence = json.loads(payload)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"hashed scene evidence is not valid JSON in {path.name}: {error}"
                ) from error
            if (evidence.get("width"), evidence.get("height")) == dimensions:
                end = evidence_indices[position + 1] if position + 1 < len(evidence_indices) \
                    else len(lines)
                session_star_line = lines[index - 1] if index > 0 else ""
                candidates.append((path, evidence, session_star_line,
                                   lines[index + 1:end]))

    require(len(candidates) == 1,
            f"{label} requires exactly one hashed typed-session scene event")
    path, evidence, session_star_line, segment = candidates[0]
    require(evidence.get("schema") == "sirius-render-scene-v1",
            f"{label} scene transcript has an unknown schema")
    require(evidence.get("backend") == "Vulkan",
            f"{label} scene transcript must identify the Vulkan backend")
    require(evidence.get("metric") == scene.get("metric"),
            "scene transcript metric does not match the attested claim")
    require_matching_number(evidence.get("spin"), scene.get("spin"), "spin")
    require(evidence.get("width") == render.get("width")
            and evidence.get("height") == render.get("height"),
            "scene transcript dimensions do not match the attested claim")
    require(evidence.get("samples_per_pixel") == scene.get("samples_per_pixel"),
            "scene transcript sample count does not match the attested claim")
    require_matching_number(evidence.get("field_of_view"), scene.get("field_of_view"),
                            "field_of_view")
    require(evidence.get("disk_enabled") is scene.get("disk_enabled"),
            "scene transcript disk state does not match the attested claim")
    require(evidence.get("ray_bundles") is scene.get("ray_bundles"),
            "scene transcript ray-bundle state does not match the attested claim")
    require(evidence.get("point_starfield") is scene.get("point_starfield"),
            "scene transcript point-star state does not match the attested claim")
    star_count = evidence.get("point_star_count")
    require(type(star_count) is int and star_count >= scene.get("star_catalogue_minimum", 0),
            "scene transcript does not prove the attested point-star catalogue floor")
    require_matching_number(evidence.get("point_brightness_scale"),
                            scene.get("point_brightness_scale"),
                            "point_brightness_scale")
    beta = evidence.get("camera_beta")
    claimed_beta = scene.get("camera_beta")
    require(isinstance(beta, list) and isinstance(claimed_beta, list)
            and len(beta) == len(claimed_beta) == 3,
            "scene transcript camera beta is malformed")
    for index, (actual, expected) in enumerate(zip(beta, claimed_beta)):
        require_matching_number(actual, expected, f"camera_beta[{index}]")
    require(evidence.get("lens") == scene.get("lens"),
            "scene transcript lens does not match the attested claim")
    for field in ("focal_length", "aperture", "focus_distance"):
        require_matching_number(evidence.get(field), scene.get(field), field)

    session_star_match = re.fullmatch(
        r"\[Session\] Point-source star field: (\d+) stars, \d+ KiB index, beams on",
        session_star_line,
    )
    require(session_star_match is not None
            and int(session_star_match.group(1)) == star_count,
            f"{label} transcript does not bind the session point-star catalogue "
            "to the scene event")
    vulkan_stars = []
    for line in segment:
        match = re.fullmatch(
            r"\[Vulkan\] point-source star field: (\d+) stars, \d+ KiB index", line
        )
        if match:
            vulkan_stars.append(int(match.group(1)))
    require(vulkan_stars == [star_count],
            f"{label} transcript does not bind the point-star catalogue to Vulkan dispatch")

    device_lines = [line for line in segment if line.startswith("[Vulkan] device[")]
    require(len(device_lines) == 1 and str(device.get("name", "")) in device_lines[0],
            f"{label} transcript device does not match the attested physical device")
    budget_lines = [line for line in segment if line.startswith("[Vulkan] budget: ")]
    require(len(budget_lines) == 1
            and re.match(r"\[Vulkan\] budget: 2048 MiB, tile: \d+x\d+ ", budget_lines[0]),
            f"{label} transcript does not prove the 2048 MiB Vulkan budget")
    require(segment.count("[Vulkan] precision: fp32 rung") == 1,
            f"{label} transcript does not identify the governed precision rung")

    completion_pattern = re.compile(
        r"\[Session\] Vulkan render complete: Kerr on (.+), (\d+) tile\(s\) of "
        r"(\d+)px in (\d+) governed dispatch\(es\), ([0-9]+(?:\.[0-9]+)?)s"
    )
    completions = []
    for index, line in enumerate(segment):
        match = completion_pattern.fullmatch(line)
        if match:
            completions.append((index, match))
    require(len(completions) == 1,
            f"{label} transcript lacks one exact Vulkan completion event")
    completion_index, completion = completions[0]
    require(completion.group(1) == str(device.get("name", "")),
            f"{label} completion device does not match the attested physical device")
    require(all(int(completion.group(index)) > 0 for index in (2, 3, 4)),
            f"{label} completion reports an empty tile or dispatch ledger")
    rendered_seconds = float(completion.group(5))
    claimed_seconds = float(render.get("wall_seconds"))
    require(rendered_seconds > 0.0 and claimed_seconds + 2.0 >= rendered_seconds
            and claimed_seconds <= rendered_seconds + 120.0,
            f"{label} measured wall time is inconsistent with the renderer completion event")

    wrote_indices = []
    for index, line in enumerate(segment):
        if line.startswith("[Session] Wrote: "):
            written_name = (
                line[len("[Session] Wrote: "):].replace("\\", "/").rsplit("/", 1)[-1]
            )
            if written_name == png_name:
                wrote_indices.append(index)
    state_indices = [
        index for index, line in enumerate(segment)
        if line == "[Session] Finished with state: Complete"
    ]
    require(len(wrote_indices) == 1 and len(state_indices) == 1
            and completion_index < wrote_indices[0] < state_indices[0],
            f"{label} transcript {path.name} does not bind completion to {png_name}")


def require_governed_scene_claim(render, dimensions, label):
    require(isinstance(render, dict), f"{label} render claim is required")
    require((render.get("width"), render.get("height")) == dimensions,
            f"{label} requires an exact {dimensions[0]}x{dimensions[1]} render")
    require(render.get("memory_budget_mib") == 2048,
            f"{label} requires the 2048 MiB governed budget")
    require(type(render.get("wall_seconds")) in (int, float)
            and math.isfinite(render["wall_seconds"]) and render["wall_seconds"] > 0,
            f"{label} requires a positive measured wall time")
    scene = render.get("scene", {})
    require(scene.get("metric") == "Kerr", f"{label} scene must use Kerr")
    require(type(scene.get("spin")) in (int, float)
            and math.isfinite(scene["spin"]) and 0.0 < scene["spin"] < 1.0,
            f"{label} scene requires finite nonzero sub-extremal spin")
    require(type(scene.get("samples_per_pixel")) is int
            and scene["samples_per_pixel"] > 0,
            f"{label} scene requires a positive exact sample count")
    require(scene.get("ray_bundles") is True,
            f"{label} scene must propagate ray bundles")
    require(scene.get("point_starfield") is True,
            f"{label} scene must render the filtered point-star catalogue")
    require(scene.get("disk_enabled") is False,
            f"{label} scene must isolate the filtered point-star output")
    require(type(scene.get("star_catalogue_minimum")) is int
            and scene["star_catalogue_minimum"] >= 100000,
            f"{label} scene requires at least 100000 point stars")
    require(type(scene.get("point_brightness_scale")) in (int, float)
            and math.isfinite(scene["point_brightness_scale"])
            and scene["point_brightness_scale"] > 0,
            f"{label} scene requires positive point-star display calibration")
    beta = scene.get("camera_beta")
    require(isinstance(beta, list) and len(beta) == 3
            and all(type(component) in (int, float) and math.isfinite(component)
                    for component in beta)
            and 0.0 < sum(component * component for component in beta) < 1.0,
            f"{label} scene requires a finite nonzero timelike camera beta")
    require(scene.get("lens") == "ThinLens",
            f"{label} scene must compose motion with the finite-aperture lens")
    for field in ("focal_length", "aperture", "focus_distance"):
        require(type(scene.get(field)) in (int, float)
                and math.isfinite(scene[field]) and scene[field] > 0,
                f"{label} scene requires a positive finite {field}")
    require(type(scene.get("field_of_view")) in (int, float)
            and math.isfinite(scene["field_of_view"])
            and 1.0 <= scene["field_of_view"] <= 170.0,
            f"{label} scene requires a finite represented field of view")


def require_sparse_point_image(inspection, label):
    require(min(inspection["channel_ranges"]) >= 32,
            f"{label} PNG has insufficient RGB dynamic range")
    require(inspection["distinct_colours"] >= 64,
            f"{label} PNG has insufficient colour structure")
    require(inspection["dark_channel_fraction"] >= 0.90,
            f"{label} PNG lacks the required sparse dark-field structure")
    require(1.0e-6 <= inspection["bright_channel_fraction"] <= 0.05,
            f"{label} PNG has an implausible resolved-source fraction")
    require(0.005 <= inspection["bright_row_fraction"] <= 0.95,
            f"{label} PNG lacks spatially distributed resolved sources")
    require(inspection["distinct_row_signatures"] >= 32,
            f"{label} PNG lacks two-dimensional row structure")


def require_runtime_probe_image(inspection, label):
    width, height = inspection["dimensions"]
    require(width >= 64 and height >= 64,
            f"{label} PNG is too small to witness a runtime frame")
    require(max(inspection["channel_ranges"]) >= 8,
            f"{label} PNG has no meaningful RGB dynamic range")
    require(inspection["distinct_colours"] >= 8,
            f"{label} PNG has no meaningful colour structure")
    require(inspection["distinct_row_signatures"] >= 2,
            f"{label} PNG is spatially collapsed")


def inspect_build_alignment_receipt(path, source_revision):
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError("native build alignment receipt is not valid JSON") from error
    require(
        receipt.get("schema_version") == 1
        and receipt.get("method") == "revision-bound-attestation-set",
        "native build alignment receipt has an unsupported schema or method",
    )
    require(receipt.get("source_revision") == source_revision,
            "native build receipt names a different source revision")
    require(receipt.get("source_tree_clean") is True,
            "native build receipt was not configured from a clean source tree")
    mode = receipt.get("alignment_mode")
    require(mode == "qualification",
            "external evidence must come from a strict qualification build")
    required = receipt.get("required_domains")
    require(isinstance(required, list) and len(required) == len(set(required))
            and set(required) == KNOWN_DOMAINS,
            "native build receipt does not name the exact release domains")
    model = receipt.get("operating_model", {})
    digest = model.get("sha256")
    require(model.get("schema_version") == 3
            and isinstance(digest, str)
            and re.fullmatch(r"[0-9a-f]{64}", digest) is not None
            and model.get("external_domains") == required,
            "native build receipt is not bound to the operating model")
    admitted = receipt.get("admitted_domains")
    pending = receipt.get("pending_domains")
    require(isinstance(admitted, list) and isinstance(pending, list)
            and len(admitted) == len(set(admitted))
            and len(pending) == len(set(pending))
            and not set(admitted).intersection(pending)
            and set(admitted).union(pending) == KNOWN_DOMAINS,
            "native build receipt domains do not partition the release obligation")
    require(isinstance(receipt.get("records"), list),
            "native build receipt has no admitted-record ledger")
    ideal = receipt.get("ultimate_ideal", {})
    satisfied = not pending
    require(ideal.get("required") is True
            and ideal.get("satisfied") is satisfied
            and ideal.get("state") == ("aligned" if satisfied else "attestation_required"),
            "native build receipt has no ultimate-ideal obligation")
    return receipt


def inspect_qualification_build_gate(
    path, source_revision, alignment_path, report_path, inventory_path, claims,
    artifacts,
):
    try:
        gate = json.loads(path.read_text(encoding="utf-8"))
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError("qualification build-gate evidence is not valid JSON") from error
    require(
        isinstance(gate, dict)
        and set(gate) == {
            "schema_version", "status", "alignment_mode", "source", "ctest",
            "inputs", "tested_artifacts", "product_artifacts",
        }
        and gate.get("schema_version") == 1
        and gate.get("status") == "passed"
        and gate.get("alignment_mode") == "qualification",
        "qualification build gate has an unsupported schema, state, or mode",
    )
    require(gate.get("source") == {"revision": source_revision, "clean": True},
            "qualification build gate is not bound to the clean attested revision")

    def artifact_record(record):
        return (
            isinstance(record, dict)
            and set(record) == {"root", "path", "bytes", "sha256"}
            and record.get("root") in {"source", "build"}
            and isinstance(record.get("path"), str) and record["path"]
            and type(record.get("bytes")) is int and record["bytes"] > 0
            and isinstance(record.get("sha256"), str)
            and re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None
        )

    tested = gate.get("tested_artifacts")
    products = gate.get("product_artifacts")
    require(isinstance(tested, dict) and set(tested) == QUALIFICATION_TESTED_ARTIFACTS
            and all(artifact_record(record) for record in tested.values()),
            "qualification build gate does not bind the exact test executables")
    require(isinstance(products, dict) and set(products) == QUALIFICATION_PRODUCT_ARTIFACTS
            and all(artifact_record(record) for record in products.values()),
            "qualification build gate does not bind the complete candidate product")

    alignment_payload = alignment_path.read_bytes()
    alignment_digest = hashlib.sha256(alignment_payload).hexdigest()
    alignment_record = products["alignment_receipt"]
    inputs = gate.get("inputs")
    require(isinstance(inputs, dict)
            and set(inputs) == {"operating_model_sha256", "alignment_receipt_sha256"}
            and inputs.get("alignment_receipt_sha256") == alignment_digest
            and alignment_record["sha256"] == alignment_digest
            and alignment_record["bytes"] == len(alignment_payload),
            "qualification build gate is not bound to the attested alignment receipt")
    alignment = json.loads(alignment_payload)
    require(inputs.get("operating_model_sha256")
            == alignment.get("operating_model", {}).get("sha256")
            == products["operating_model"]["sha256"],
            "qualification build gate differs from the attested operating model")
    for product_name, artifact_name in QUALIFICATION_PRODUCT_EVIDENCE.items():
        product_path = artifacts.get(artifact_name)
        product_payload = product_path.read_bytes() if product_path is not None else b""
        require(product_path is not None
                and products[product_name]["bytes"] == len(product_payload)
                and products[product_name]["sha256"]
                == hashlib.sha256(product_payload).hexdigest(),
                f"qualification product is absent or differs from its gate: {product_name}")

    report = inspect_junit(report_path)
    gate_junit_path = artifacts.get("qualification-gate-junit")
    gate_log_path = artifacts.get("qualification-gate-log")
    gate_junit_payload = (
        gate_junit_path.read_bytes() if gate_junit_path is not None else b""
    )
    gate_log_payload = gate_log_path.read_bytes() if gate_log_path is not None else b""
    gate_report = inspect_junit(gate_junit_path) if gate_junit_path is not None else {}
    ctest = gate.get("ctest")
    require(isinstance(ctest, dict)
            and set(ctest) == {
                "inventory_sha256", "junit", "log", "registered", "executed",
                "failures", "errors", "skipped",
            }
            and artifact_record(ctest["junit"])
            and artifact_record(ctest["log"])
            and gate_junit_path is not None
            and gate_log_path is not None
            and ctest["junit"]["bytes"] == len(gate_junit_payload)
            and ctest["junit"]["sha256"]
            == hashlib.sha256(gate_junit_payload).hexdigest()
            and ctest["log"]["bytes"] == len(gate_log_payload)
            and ctest["log"]["sha256"] == hashlib.sha256(gate_log_payload).hexdigest()
            and gate_report == report
            and junit_case_names(gate_junit_path) == junit_case_names(report_path)
            and ctest["registered"] == ctest["executed"] == gate_report["cases"]
            and report["cases"] >= MINIMUM_SOURCE_AVAILABLE_TESTS
            and ctest["failures"] == ctest["errors"] == ctest["skipped"] == 0,
            "qualification build gate does not prove the attested zero-skip estate")
    canonical_inventory = json.dumps(
        inventory, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    require(ctest["inventory_sha256"] == hashlib.sha256(canonical_inventory).hexdigest(),
            "qualification build gate differs from the attested CTest inventory")

    candidate = claims.get("qualification_executable")
    candidate_path = artifacts.get("qualification-sirius.bin")
    candidate_payload = candidate_path.read_bytes() if candidate_path is not None else b""
    require(isinstance(candidate, dict)
            and set(candidate) == {"artifact", "bytes", "sha256"}
            and candidate.get("artifact") == "qualification-sirius.bin"
            and candidate_path is not None
            and candidate == {
                "artifact": "qualification-sirius.bin",
                "bytes": products["sirius"]["bytes"],
                "sha256": products["sirius"]["sha256"],
            }
            and len(candidate_payload) == candidate["bytes"]
            and hashlib.sha256(candidate_payload).hexdigest() == candidate["sha256"]
            and tested["sirius"]["bytes"] == products["sirius"]["bytes"]
            and tested["sirius"]["sha256"] == products["sirius"]["sha256"],
            "attestation is not bound to the qualification executable")
    return gate


def verify_document(data, location):
    require(data.get("schema_version") == 1, "schema_version must be 1")
    require(data.get("status") == "pass", "status must be pass")
    domains = data.get("domains")
    require(isinstance(domains, list) and domains, "domains must be a non-empty list")
    require(all(isinstance(domain, str) and domain for domain in domains),
            "every domain must be a non-empty string")
    require(len(domains) == len(set(domains)), "domains must be unique")
    unknown = set(domains) - KNOWN_DOMAINS
    require(not unknown, f"unknown attestation domains: {sorted(unknown)}")

    source_revision = data.get("source_revision")
    require(isinstance(source_revision, str)
            and re.fullmatch(r"[0-9a-f]{40,64}", source_revision) is not None,
            "source_revision must be a full lowercase hexadecimal revision")
    completed = data.get("completed_utc")
    require(isinstance(completed, str), "completed_utc is required")
    completed_at = datetime.fromisoformat(completed.replace("Z", "+00:00"))
    require(completed_at.utcoffset() == timezone.utc.utcoffset(completed_at),
            "completed_utc must carry an explicit UTC offset")

    platform = data.get("platform", {})
    device = data.get("device", {})
    claims = data.get("claims", {})
    require(claims.get("test_estate_passed") is True, "test estate must have passed")

    physical_domains = {
        "physical-radeon-780m",
        "wsl2-dozen",
        "windows-native-vulkan",
        "macos-moltenvk",
        "physical-imax-5616x4096",
    }
    runtime_device_domains = physical_domains | {"viewer-native-window-input"}
    if runtime_device_domains.intersection(domains):
        require(device.get("kind") in {"integrated", "discrete"},
                "runtime domains require a physical Vulkan device")
        require(isinstance(device.get("name"), str) and device["name"].strip(),
                "runtime domains require a non-empty device name")
        require(isinstance(device.get("driver_name"), str) and device["driver_name"].strip(),
                "runtime domains require a non-empty driver name")
        require_positive_integer(device.get("vendor_id"), "vendor_id")
        require_positive_integer(device.get("device_id"), "device_id")
        require_positive_integer(device.get("api_version"), "api_version")
        require(type(device.get("driver_id")) is int and device["driver_id"] >= 0,
                "runtime device driver_id must be a non-negative integer")
        require_positive_integer(device.get("device_local_bytes"), "device_local_bytes")
        require_positive_integer(device.get("render_memory_bytes"), "render_memory_bytes")
        require(type(device.get("supports_fp64")) is bool,
                "runtime device supports_fp64 must be a boolean")
        require(claims.get("runtime_ready") is True,
                "runtime domains require a successful readiness probe")

    if "physical-radeon-780m" in domains or "viewer-native-window-input" in domains:
        require("radeon 780m" in str(device.get("name", "")).casefold(),
                "Radeon 780M domain requires an exact device-name witness")
    if "wsl2-dozen" in domains:
        require(platform.get("wsl2") is True, "Dozen domain requires WSL2")
        require(text_contains(device, "dozen", "dzn"),
                "Dozen domain requires driver identity containing Dozen/dzn")
    if "windows-native-build" in domains or "windows-native-vulkan" in domains:
        require(str(platform.get("os", "")).casefold() == "windows",
                "Windows domains require a native Windows host")
    if "windows-native-vulkan" in domains:
        require(not text_contains(device, "dozen", "dzn"),
                "native Windows Vulkan may not be attested through Dozen")
    if "macos-native-build" in domains or "macos-moltenvk" in domains:
        require(str(platform.get("os", "")).casefold() == "macos",
                "macOS domains require a native macOS host")
    if "macos-moltenvk" in domains:
        require(text_contains(device, "moltenvk"),
                "MoltenVK domain requires MoltenVK driver identity")
    if "viewer-native-window-input" in domains:
        viewer = claims.get("viewer", {})
        require(viewer.get("backend") == "Vulkan",
                "native viewer evidence must use the Vulkan backend")
        require(viewer.get("window_created") is True, "viewer window was not created")
        require(viewer.get("frame_observed") is True,
                "no progressive Vulkan viewer frame was observed")
        require(type(viewer.get("published_frame_count")) is int
                and viewer["published_frame_count"] > 0,
                "native viewer evidence requires a positive published-frame count")
        require(viewer.get("keyboard_event_observed") is True,
                "no native keyboard event was observed")
        require(viewer.get("pointer_event_observed") is True,
                "no native pointer event was observed")
    if "physical-imax-5616x4096" in domains:
        require_governed_scene_claim(
            claims.get("p3_render_1080p"), (1920, 1080), "1080p P3/P5 domain"
        )
        require_governed_scene_claim(
            claims.get("imax_render"), (5616, 4096), "IMAX P3/P5 domain"
        )

    artifacts = data.get("artifacts")
    require(isinstance(artifacts, dict) and artifacts, "artifacts must be a non-empty object")
    base = location.parent
    artifact_paths = []
    artifact_files = {}
    resolved_artifact_paths = set()
    for name, record in artifacts.items():
        require(isinstance(name, str) and name, "artifact names must be non-empty strings")
        require(isinstance(record, dict), f"artifact {name} must be an object")
        relative = record.get("path", name)
        require(isinstance(relative, str) and relative,
                f"artifact {name} path must be a non-empty string")
        path = Path(relative)
        require(not path.is_absolute() and ".." not in path.parts,
                f"artifact {name} path must stay inside the attestation bundle")
        path = base / path
        resolved_path = path.resolve()
        try:
            resolved_path.relative_to(base.resolve())
        except ValueError as error:
            raise ValueError(f"artifact {name} escapes the attestation bundle") from error
        require(resolved_path not in resolved_artifact_paths,
                f"artifact {name} aliases another artifact path")
        resolved_artifact_paths.add(resolved_path)
        require(path.is_file(), f"artifact {name} is missing: {path}")
        payload = path.read_bytes()
        recorded_bytes = record.get("bytes")
        require(type(recorded_bytes) is int and recorded_bytes > 0,
                f"artifact {name} byte count must be a positive integer")
        require(len(payload) == recorded_bytes, f"artifact {name} size mismatch")
        recorded_digest = record.get("sha256")
        require(isinstance(recorded_digest, str)
                and len(recorded_digest) == 64
                and all(character in "0123456789abcdef" for character in recorded_digest),
                f"artifact {name} SHA-256 must be 64 lowercase hexadecimal characters")
        digest = hashlib.sha256(payload).hexdigest()
        require(digest == recorded_digest, f"artifact {name} hash mismatch")
        artifact_paths.append(path)
        artifact_files[name] = path

    build_domains = {"windows-native-build", "macos-native-build"}
    if build_domains.intersection(domains):
        reports = [path for path in artifact_paths if path.suffix.casefold() == ".xml"]
        require(len(reports) == 1, "native build attestation requires exactly one JUnit XML report")
        report = inspect_junit(reports[0])
        require(report["cases"] >= MINIMUM_SOURCE_AVAILABLE_TESTS
                and report["skipped"] == 0,
                "native build JUnit evidence is not a non-skipping test estate")
        require(claims.get("test_report") == report,
                "native build JUnit summary does not match the attested claim")
        case_names = junit_case_names(reports[0])
        require(
            "AlignmentAuthority.CompiledReceiptMatchesTheStagedRuntimeAuthority"
            in case_names,
            "native build JUnit does not prove the compiled alignment authority",
        )
        json_artifacts = load_json_artifacts(artifact_paths)
        receipts = [
            path for path, document in json_artifacts.items()
            if isinstance(document, dict)
            and document.get("method") == "revision-bound-attestation-set"
        ]
        inventories = [
            path for path, document in json_artifacts.items()
            if isinstance(document, dict) and document.get("kind") == "ctestInfo"
        ]
        gates = [
            path for path, document in json_artifacts.items()
            if isinstance(document, dict)
            and document.get("status") == "passed"
            and "tested_artifacts" in document
            and "product_artifacts" in document
        ]
        require(len(json_artifacts) == 3,
                "native build attestation has an unclassified JSON artifact")
        require(len(receipts) == 1,
                "native build attestation requires exactly one alignment receipt")
        require(len(inventories) == 1,
                "native build attestation requires exactly one CTest inventory")
        require(len(gates) == 1,
                "native build attestation requires exactly one qualification build gate")
        inspect_build_alignment_receipt(receipts[0], source_revision)
        inspect_qualification_build_gate(
            gates[0], source_revision, receipts[0], reports[0], inventories[0], claims,
            artifact_files,
        )
        require(case_names == inspect_ctest_inventory(inventories[0]),
                "native build JUnit does not equal the registered CTest estate")
    hardware_floor_domains = {
        "physical-radeon-780m",
        "wsl2-dozen",
        "physical-imax-5616x4096",
        "viewer-native-window-input",
    }
    if runtime_device_domains.intersection(domains):
        reports = [path for path in artifact_paths if path.suffix.casefold() == ".xml"]
        require(len(reports) == 1,
                "runtime attestation requires exactly one JUnit XML report")
        report = inspect_junit(reports[0])
        require(report["cases"] > 0 and report["skipped"] == 0,
                "runtime JUnit evidence is not a non-skipping test estate")
        require(report["cases"] >= MINIMUM_SOURCE_AVAILABLE_TESTS,
                "runtime JUnit evidence is below the source-available estate floor")
        require(claims.get("test_report") == report,
                "runtime JUnit summary does not match the attested claim")
        runtime_case_names = junit_case_names(reports[0])
        require(
            "AlignmentAuthority.CompiledReceiptMatchesTheStagedRuntimeAuthority"
            in runtime_case_names,
            "runtime JUnit does not prove the compiled alignment authority",
        )
        if hardware_floor_domains.intersection(domains):
            hardware_case_names = {
                case.get("name") for case in ET.parse(reports[0]).getroot().iter("testcase")
            }
            require(
                "ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames"
                in hardware_case_names,
                "hardware JUnit does not prove interactive progressive Vulkan refinement",
            )

        json_artifacts = load_json_artifacts(artifact_paths)
        inventories = [
            (path, document) for path, document in json_artifacts.items()
            if isinstance(document, dict) and "backends" in document
        ]
        receipts = [
            path for path, document in json_artifacts.items()
            if isinstance(document, dict)
            and document.get("method") == "revision-bound-attestation-set"
        ]
        ctest_inventories = [
            path for path, document in json_artifacts.items()
            if isinstance(document, dict) and document.get("kind") == "ctestInfo"
        ]
        gates = [
            path for path, document in json_artifacts.items()
            if isinstance(document, dict)
            and document.get("status") == "passed"
            and "tested_artifacts" in document
            and "product_artifacts" in document
        ]
        require(len(json_artifacts) == 4,
                "runtime attestation has an unclassified JSON artifact")
        require(len(inventories) == 1,
                "runtime attestation requires one hashed device inventory")
        require(len(receipts) == 1,
                "runtime attestation requires one compiled-authority alignment receipt")
        require(len(ctest_inventories) == 1,
                "runtime attestation requires one registered CTest inventory")
        require(len(gates) == 1,
                "runtime attestation requires one qualification build gate")
        inventory = inventories[0][1]
        inspect_build_alignment_receipt(receipts[0], source_revision)
        inspect_qualification_build_gate(
            gates[0], source_revision, receipts[0], reports[0],
            ctest_inventories[0], claims, artifact_files
        )
        require(runtime_case_names == inspect_ctest_inventory(ctest_inventories[0]),
                "runtime JUnit does not equal the registered CTest estate")
        inventory_vulkan = inventory.get("backends", {}).get("vulkan", {})
        device_index = device.get("index")
        require(type(device_index) is int and device_index >= 0,
                "runtime device inventory requires a non-negative selected index")
        require(inventory_vulkan.get("selected_device_index") == device_index,
                "hashed inventory selected a different Vulkan device")
        inventory_devices = inventory_vulkan.get("devices")
        require(isinstance(inventory_devices, list)
                and sum(item == device for item in inventory_devices) == 1,
                "attested device does not match the hashed inventory")
        require(str(inventory.get("platform", "")).casefold()
                == str(platform.get("os", "")).casefold()
                and inventory.get("wsl2") is platform.get("wsl2"),
                "attested platform does not match the hashed inventory")

        transcripts = [
            path for path in artifact_paths if path.suffix.casefold() in {".log", ".txt"}
        ]
        transcript_lines = []
        for path in transcripts:
            transcript_lines.extend(
                path.read_text(encoding="utf-8", errors="replace").splitlines()
            )
        require(transcript_lines.count(f"== source revision: {source_revision}") == 1,
                "runtime transcript does not bind the source revision")
        require(
            transcript_lines.count(
                "== readiness: evidence_generation=true vulkan=true "
                f"selected_device_index={device_index}"
            ) == 1,
            "runtime transcript does not bind readiness to the selected device",
        )
        preset = claims.get("preset")
        if "windows-native-vulkan" in domains:
            expected_preset = "windows-msvc"
        elif "macos-moltenvk" in domains:
            expected_preset = "macos"
        else:
            expected_preset = "linux-ci"
        require(
            preset == expected_preset
            and transcript_lines.count(
                f"== [1/6] configure + build ({expected_preset})"
            ) == 1,
            f"runtime attestation did not execute the {expected_preset} profile",
        )
        test_summaries = []
        for line in transcript_lines:
            match = re.fullmatch(
                r"100% tests passed, 0 tests failed out of (\d+)", line.strip()
            )
            if match:
                test_summaries.append(int(match.group(1)))
        require(test_summaries.count(report["cases"]) == 1,
                "runtime transcript does not bind the complete green JUnit estate")
    png_inspections = {}
    if physical_domains.intersection(domains):
        pngs = [path for path in artifact_paths if path.suffix.casefold() == ".png"]
        transcripts = [
            path for path in artifact_paths if path.suffix.casefold() in {".log", ".txt"}
        ]
        require(pngs, "runtime attestation requires a hashed PNG output")
        require(transcripts, "runtime attestation requires a hashed execution transcript")
        for path in pngs:
            inspection = inspect_png(path)
            require_runtime_probe_image(inspection, "runtime probe")
            png_inspections[path] = inspection
    if "viewer-native-window-input" in domains:
        transcripts = [
            path for path in artifact_paths if path.suffix.casefold() in {".log", ".txt"}
        ]
        require(transcripts, "native viewer attestation requires a hashed input transcript")
        observed_frame_counts = []
        for path in transcripts:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            if (
                any(line.startswith("window-created ") for line in lines)
                and any(line.startswith("keyboard-callback ") for line in lines)
                and any(line.startswith("pointer-callback ") for line in lines)
            ):
                frames = [
                    line for line in lines
                    if re.fullmatch(
                        r"frame-published backend=Vulkan width=\d+ height=\d+", line
                    )
                ]
                if frames:
                    observed_frame_counts.append(len(frames))
        require(
            claims["viewer"]["published_frame_count"] in observed_frame_counts,
            "no single hashed transcript binds a Vulkan frame to native window and input callbacks",
        )
    if "physical-imax-5616x4096" in domains:
        pngs = [path for path in artifact_paths if path.suffix.casefold() == ".png"]
        require(pngs, "P3/P5 domain requires hashed 1080p and IMAX PNG artifacts")
        matching = {(1920, 1080): [], (5616, 4096): []}
        for path in pngs:
            inspection = png_inspections[path]
            if inspection["dimensions"] in matching:
                matching[inspection["dimensions"]].append((path, inspection))
        require(len(matching[(1920, 1080)]) == 1,
                "P3 domain requires exactly one hashed 1920x1080 PNG artifact")
        require(len(matching[(5616, 4096)]) == 1,
                "IMAX domain requires exactly one hashed 5616x4096 PNG artifact")
        p3_1080 = matching[(1920, 1080)][0]
        imax = matching[(5616, 4096)][0]
        require_sparse_point_image(p3_1080[1], "1080p P3/P5")
        require_sparse_point_image(imax[1], "IMAX P3/P5")
        transcripts = [
            path for path in artifact_paths if path.suffix.casefold() in {".log", ".txt"}
        ]
        verify_governed_scene_transcript(
            transcripts,
            claims["p3_render_1080p"]["scene"],
            claims["p3_render_1080p"],
            p3_1080[0].name,
            device,
            (1920, 1080),
            "1080p P3/P5",
        )
        verify_governed_scene_transcript(
            transcripts,
            claims["imax_render"]["scene"],
            claims["imax_render"],
            imax[0].name,
            device,
            (5616, 4096),
            "IMAX P3/P5",
        )


def verify_path(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    verify_document(data, path)


def load_json_artifacts(paths):
    documents = {}
    for path in paths:
        if path.suffix.casefold() != ".json":
            continue
        try:
            documents[path] = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            raise ValueError(f"hashed JSON artifact is invalid: {path.name}") from error
    return documents


def inspect_junit(path):
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as error:
        raise ValueError(f"build evidence is not valid JUnit XML: {error}") from error
    cases = list(root.iter("testcase"))
    require(cases, "build evidence contains no JUnit test cases")
    failures = sum(1 for case in cases if case.find("failure") is not None)
    errors = sum(1 for case in cases if case.find("error") is not None)
    skipped = sum(1 for case in cases if case.find("skipped") is not None)
    require(failures == 0 and errors == 0,
            f"build evidence contains {failures} failure(s) and {errors} error(s)")
    return {"cases": len(cases), "failures": failures, "errors": errors, "skipped": skipped}


def junit_case_names(path):
    names = [case.get("name") for case in ET.parse(path).getroot().iter("testcase")]
    require(all(isinstance(name, str) and name for name in names),
            "JUnit contains a testcase without a name")
    require(len(names) == len(set(names)), "JUnit testcase names are not unique")
    return set(names)


def inspect_ctest_inventory(path):
    try:
        inventory = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError("CTest inventory is not valid JSON") from error
    require(inventory.get("kind") == "ctestInfo"
            and inventory.get("version", {}).get("major") == 1,
            "CTest inventory has an unsupported schema")
    tests = inventory.get("tests")
    require(isinstance(tests, list) and tests,
            "CTest inventory contains no registered tests")
    enabled = []
    all_names = []
    for test in tests:
        require(isinstance(test, dict) and isinstance(test.get("name"), str)
                and test["name"], "CTest inventory contains a malformed test")
        name = test["name"]
        all_names.append(name)
        disabled = False
        properties = test.get("properties", [])
        require(isinstance(properties, list),
                f"CTest inventory properties are malformed for {name}")
        for prop in properties:
            require(isinstance(prop, dict),
                    f"CTest inventory property is malformed for {name}")
            if prop.get("name") == "DISABLED":
                disabled = prop.get("value") in (True, 1, "1", "ON", "TRUE")
        if not disabled:
            enabled.append(name)
    require(len(all_names) == len(set(all_names)),
            "CTest inventory test names are not unique")
    require(enabled, "CTest inventory contains no enabled tests")
    return set(enabled)


def png_chunk(chunk_type, payload):
    return (struct.pack(">I", len(payload)) + chunk_type + payload
            + struct.pack(">I", zlib.crc32(chunk_type + payload) & 0xffffffff))


def write_self_test_png(path, width, height, structured):
    if structured:
        background = bytes((38, 38, 38, 255)) * width
    else:
        background = bytes((249, 255, 255, 255)) * width
    compressor = zlib.compressobj(level=1)
    compressed = bytearray()
    for y in range(height):
        row = background
        if structured and y % 8 == 0:
            row = bytearray(background)
            for star in range(8):
                x = (y * 131 + star * 683) % width
                red = 96 + (y + star * 23) % 160
                green = 96 + (y * 3 + star * 17) % 144
                blue = 96 + (y * 7 + star * 11) % 128
                row[x * 4:x * 4 + 4] = bytes((red, green, blue, 255))
        compressed.extend(compressor.compress(b"\0" + row))
    compressed.extend(compressor.flush())
    header = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + png_chunk(b"IHDR", header)
        + png_chunk(b"IDAT", bytes(compressed))
        + png_chunk(b"IEND", b"")
    )


def self_test():
    with tempfile.TemporaryDirectory(prefix="sirius-attestation-") as directory:
        root = Path(directory)
        source_revision = "0" * 40
        required_domains = sorted(KNOWN_DOMAINS)
        self_test_products = {}
        for logical_name, artifact_name in QUALIFICATION_PRODUCT_EVIDENCE.items():
            product_path = root / artifact_name
            product_path.write_bytes(f"qualification product: {logical_name}\n".encode())
            self_test_products[logical_name] = product_path
        self_test_model_digest = hashlib.sha256(
            self_test_products["operating_model"].read_bytes()
        ).hexdigest()
        native_build_receipt = root / "alignment_receipt.json"
        native_build_receipt.write_text(
            json.dumps({
                "schema_version": 1,
                "method": "revision-bound-attestation-set",
                "alignment_mode": "qualification",
                "source_revision": source_revision,
                "source_tree_clean": True,
                "operating_model": {
                    "schema_version": 3,
                    "sha256": self_test_model_digest,
                    "external_domains": required_domains,
                },
                "required_domains": required_domains,
                "admitted_domains": [],
                "pending_domains": required_domains,
                "records": [],
                "ultimate_ideal": {
                    "required": True,
                    "satisfied": False,
                    "state": "attestation_required",
                },
            }),
            encoding="utf-8",
        )
        self_test_device = {
            "index": 0,
            "name": "AMD Radeon 780M",
            "kind": "integrated",
            "driver_name": "RADV",
            "driver_info": "self-test",
            "vendor_id": 4098,
            "device_id": 5567,
            "api_version": 4206592,
            "driver_id": 3,
            "device_local_bytes": 8589934592,
            "render_memory_bytes": 2147483648,
            "supports_fp64": True,
        }
        artifact = root / "frame.png"
        artifact_1080 = root / "frame-1080.png"
        write_self_test_png(artifact, 5616, 4096, structured=True)
        write_self_test_png(artifact_1080, 1920, 1080, structured=True)
        test_report = root / "hardware-tests.xml"
        e3_case = "ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames"
        alignment_case = (
            "AlignmentAuthority.CompiledReceiptMatchesTheStagedRuntimeAuthority"
        )
        test_case_names = [
            e3_case if index == 0 else alignment_case if index == 1 else f"case-{index}"
            for index in range(MINIMUM_SOURCE_AVAILABLE_TESTS)
        ]
        test_report.write_text(
            "<testsuites><testsuite>"
            + "".join(
                '<testcase classname="self-test" '
                f'name="{name}"/>'
                for name in test_case_names
            )
            + "</testsuite></testsuites>\n",
            encoding="utf-8",
        )
        ctest_inventory_path = root / "ctest-inventory.json"
        ctest_inventory_path.write_text(
            json.dumps({
                "kind": "ctestInfo",
                "version": {"major": 1, "minor": 0},
                "tests": [
                    {"name": name, "properties": []} for name in test_case_names
                ],
            }),
            encoding="utf-8",
        )
        qualification_executable = root / "qualification-sirius.bin"
        qualification_executable.write_bytes(b"strict qualification candidate\0")

        def write_qualification_gate(gate_path, report_path, inventory_file,
                                     candidate_path):
            def gate_record(name, *, root_name="build", payload=None, digest=None,
                            byte_count=None):
                content = name.encode("utf-8") if payload is None else payload
                return {
                    "root": root_name,
                    "path": name,
                    "bytes": len(content) if byte_count is None else byte_count,
                    "sha256": hashlib.sha256(content).hexdigest()
                    if digest is None else digest,
                }

            tested = {
                name: gate_record(f"tests/{name}")
                for name in QUALIFICATION_TESTED_ARTIFACTS
            }
            products = {
                name: gate_record(f"products/{name}")
                for name in QUALIFICATION_PRODUCT_ARTIFACTS
            }
            receipt_payload = native_build_receipt.read_bytes()
            products["alignment_receipt"] = gate_record(
                "generated/sirius/alignment_receipt.json", payload=receipt_payload
            )
            for logical_name, product_path in self_test_products.items():
                products[logical_name] = gate_record(
                    f"qualification/{logical_name}", payload=product_path.read_bytes()
                )
            candidate_payload = candidate_path.read_bytes()
            candidate_record = gate_record(
                "src/sirius/app/sirius", payload=candidate_payload
            )
            products["sirius"] = candidate_record
            tested["sirius"] = dict(candidate_record)
            report_payload = report_path.read_bytes()
            gate_junit_path = root / "qualification-gate-junit"
            gate_log_path = root / "qualification-gate-log"
            gate_junit_path.write_bytes(report_payload)
            gate_log_path.write_bytes(b"all mandatory tests passed\n")
            inventory = json.loads(inventory_file.read_text(encoding="utf-8"))
            canonical_inventory = json.dumps(
                inventory, sort_keys=True, separators=(",", ":")
            ).encode("utf-8")
            case_count = inspect_junit(report_path)["cases"]
            gate = {
                "schema_version": 1,
                "status": "passed",
                "alignment_mode": "qualification",
                "source": {"revision": source_revision, "clean": True},
                "ctest": {
                    "inventory_sha256": hashlib.sha256(canonical_inventory).hexdigest(),
                    "junit": gate_record(
                        "generated/sirius/mandatory_gate_junit.xml",
                        payload=gate_junit_path.read_bytes(),
                    ),
                    "log": gate_record(
                        "generated/sirius/mandatory_gate_ctest.log",
                        payload=gate_log_path.read_bytes(),
                    ),
                    "registered": case_count,
                    "executed": case_count,
                    "failures": 0,
                    "errors": 0,
                    "skipped": 0,
                },
                "inputs": {
                    "operating_model_sha256": self_test_model_digest,
                    "alignment_receipt_sha256": hashlib.sha256(
                        receipt_payload
                    ).hexdigest(),
                },
                "tested_artifacts": tested,
                "product_artifacts": products,
            }
            gate_path.write_text(json.dumps(gate), encoding="utf-8")

        mandatory_gate = root / "mandatory_gate.json"
        write_qualification_gate(
            mandatory_gate, test_report, ctest_inventory_path, qualification_executable
        )
        qualification_claim = {
            "artifact": qualification_executable.name,
            "bytes": qualification_executable.stat().st_size,
            "sha256": hashlib.sha256(qualification_executable.read_bytes()).hexdigest(),
        }
        inventory_path = root / "device.json"
        inventory_path.write_text(
            json.dumps({
                "platform": "linux",
                "wsl2": False,
                "backends": {
                    "vulkan": {
                        "selected_device_index": 0,
                        "devices": [self_test_device],
                    }
                },
            }),
            encoding="utf-8",
        )
        transcript = root / "transcript.log"
        transcript.write_text(
            f"== source revision: {source_revision}\n"
            "== [1/6] configure + build (linux-ci)\n"
            "== readiness: evidence_generation=true vulkan=true selected_device_index=0\n"
            f"100% tests passed, 0 tests failed out of {MINIMUM_SOURCE_AVAILABLE_TESTS}\n"
            '[Session] Point-source star field: 100000 stars, 902 KiB index, beams on\n'
            '[Session] Scene evidence: {"schema":"sirius-render-scene-v1",'
            '"backend":"Vulkan","metric":"Kerr","spin":0.9,"width":1920,'
            '"height":1080,"samples_per_pixel":4,"field_of_view":60.0,'
            '"disk_enabled":false,"ray_bundles":true,"point_starfield":true,'
            '"point_star_count":100000,"point_brightness_scale":100.0,'
            '"camera_beta":[0.1,0.02,-0.01],"lens":"ThinLens",'
            '"focal_length":50.0,"aperture":2.8,"focus_distance":30.0}\n'
            '[Vulkan] device[0]: AMD Radeon 780M (integrated)\n'
            '[Vulkan] budget: 2048 MiB, tile: 1920x1080 (working set 32400 KiB)\n'
            '[Vulkan] precision: fp32 rung\n'
            '[Vulkan] point-source star field: 100000 stars, 902 KiB index\n'
            '[Session] Vulkan render complete: Kerr on AMD Radeon 780M, 1 tile(s) of '
            '1920px in 3000 governed dispatch(es), 20.5s\n'
            '[Session] Wrote: frame-1080.png\n'
            '[Session] Finished with state: Complete\n'
            '[Session] Point-source star field: 100000 stars, 902 KiB index, beams on\n'
            '[Session] Scene evidence: {"schema":"sirius-render-scene-v1",'
            '"backend":"Vulkan","metric":"Kerr","spin":0.9,"width":5616,'
            '"height":4096,"samples_per_pixel":4,"field_of_view":60.0,'
            '"disk_enabled":false,"ray_bundles":true,"point_starfield":true,'
            '"point_star_count":100000,"point_brightness_scale":100.0,'
            '"camera_beta":[0.1,0.02,-0.01],"lens":"ThinLens",'
            '"focal_length":50.0,"aperture":2.8,"focus_distance":30.0}\n'
            '[Vulkan] device[0]: AMD Radeon 780M (integrated)\n'
            '[Vulkan] budget: 2048 MiB, tile: 4096x4096 (working set 262144 KiB)\n'
            '[Vulkan] precision: fp32 rung\n'
            '[Vulkan] point-source star field: 100000 stars, 902 KiB index\n'
            '[Session] Vulkan render complete: Kerr on AMD Radeon 780M, 2 tile(s) of '
            '4096px in 32736 governed dispatch(es), 100.5s\n'
            '[Session] Wrote: frame.png\n'
            '[Session] Finished with state: Complete\n',
            encoding="utf-8",
        )
        digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        digest_1080 = hashlib.sha256(artifact_1080.read_bytes()).hexdigest()
        transcript_digest = hashlib.sha256(transcript.read_bytes()).hexdigest()
        scene_claim = {
            "metric": "Kerr",
            "spin": 0.9,
            "samples_per_pixel": 4,
            "field_of_view": 60.0,
            "disk_enabled": False,
            "ray_bundles": True,
            "point_starfield": True,
            "star_catalogue_minimum": 100000,
            "point_brightness_scale": 100.0,
            "camera_beta": [0.1, 0.02, -0.01],
            "lens": "ThinLens",
            "focal_length": 50.0,
            "aperture": 2.8,
            "focus_distance": 30.0,
        }
        valid = {
            "schema_version": 1,
            "status": "pass",
            "domains": ["physical-radeon-780m", "physical-imax-5616x4096"],
            "source_revision": source_revision,
            "completed_utc": "2026-07-31T00:00:00+00:00",
            "platform": {"os": "linux", "wsl2": False},
            "device": self_test_device,
            "claims": {
                "preset": "linux-ci",
                "test_estate_passed": True,
                "test_report": {
                    "cases": MINIMUM_SOURCE_AVAILABLE_TESTS,
                    "failures": 0,
                    "errors": 0,
                    "skipped": 0,
                },
                "runtime_ready": True,
                "qualification_executable": qualification_claim,
                "p3_render_1080p": {
                    "width": 1920,
                    "height": 1080,
                    "memory_budget_mib": 2048,
                    "wall_seconds": 21.0,
                    "scene": scene_claim,
                },
                "imax_render": {
                    "width": 5616,
                    "height": 4096,
                    "memory_budget_mib": 2048,
                    "wall_seconds": 101.0,
                    "scene": scene_claim,
                },
            },
            "artifacts": {
                "frame.png": {
                    "path": "frame.png",
                    "bytes": artifact.stat().st_size,
                    "sha256": digest,
                },
                "frame-1080.png": {
                    "path": "frame-1080.png",
                    "bytes": artifact_1080.stat().st_size,
                    "sha256": digest_1080,
                },
                "transcript.log": {
                    "path": "transcript.log",
                    "bytes": transcript.stat().st_size,
                    "sha256": transcript_digest,
                },
                "hardware-tests.xml": {
                    "path": "hardware-tests.xml",
                    "bytes": test_report.stat().st_size,
                    "sha256": hashlib.sha256(test_report.read_bytes()).hexdigest(),
                },
                "device.json": {
                    "path": "device.json",
                    "bytes": inventory_path.stat().st_size,
                    "sha256": hashlib.sha256(inventory_path.read_bytes()).hexdigest(),
                },
                "alignment_receipt.json": {
                    "path": "alignment_receipt.json",
                    "bytes": native_build_receipt.stat().st_size,
                    "sha256": hashlib.sha256(
                        native_build_receipt.read_bytes()
                    ).hexdigest(),
                },
                "ctest-inventory.json": {
                    "path": "ctest-inventory.json",
                    "bytes": ctest_inventory_path.stat().st_size,
                    "sha256": hashlib.sha256(
                        ctest_inventory_path.read_bytes()
                    ).hexdigest(),
                },
                "mandatory_gate.json": {
                    "path": mandatory_gate.name,
                    "bytes": mandatory_gate.stat().st_size,
                    "sha256": hashlib.sha256(mandatory_gate.read_bytes()).hexdigest(),
                },
                "qualification-sirius.bin": {
                    "path": qualification_executable.name,
                    "bytes": qualification_executable.stat().st_size,
                    "sha256": hashlib.sha256(
                        qualification_executable.read_bytes()
                    ).hexdigest(),
                },
                "qualification-gate-junit": {
                    "path": "qualification-gate-junit",
                    "bytes": (root / "qualification-gate-junit").stat().st_size,
                    "sha256": hashlib.sha256(
                        (root / "qualification-gate-junit").read_bytes()
                    ).hexdigest(),
                },
                "qualification-gate-log": {
                    "path": "qualification-gate-log",
                    "bytes": (root / "qualification-gate-log").stat().st_size,
                    "sha256": hashlib.sha256(
                        (root / "qualification-gate-log").read_bytes()
                    ).hexdigest(),
                },
                **{
                    path.name: {
                        "path": path.name,
                        "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                    for path in self_test_products.values()
                },
            },
        }
        verify_document(valid, root / "attestation.json")

        qualification_controls = []
        missing_gate = json.loads(json.dumps(valid))
        missing_gate["artifacts"].pop("mandatory_gate.json")
        qualification_controls.append(("missing qualification gate", missing_gate))
        substituted_candidate = json.loads(json.dumps(valid))
        substituted_candidate["claims"]["qualification_executable"]["sha256"] = "f" * 64
        qualification_controls.append(("substituted qualification candidate",
                                       substituted_candidate))
        missing_gate_junit = json.loads(json.dumps(valid))
        missing_gate_junit["artifacts"].pop("qualification-gate-junit")
        qualification_controls.append(("missing gate-generated JUnit", missing_gate_junit))
        for label, candidate in qualification_controls:
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: {label}")

        original_receipt_text = native_build_receipt.read_text(encoding="utf-8")
        development_receipt = json.loads(original_receipt_text)
        development_receipt["alignment_mode"] = "development"
        native_build_receipt.write_text(json.dumps(development_receipt), encoding="utf-8")
        development_evidence = json.loads(json.dumps(valid))
        development_payload = native_build_receipt.read_bytes()
        development_evidence["artifacts"]["alignment_receipt.json"].update(
            bytes=len(development_payload),
            sha256=hashlib.sha256(development_payload).hexdigest(),
        )
        try:
            verify_document(development_evidence, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: development evidence producer")
        native_build_receipt.write_text(original_receipt_text, encoding="utf-8")

        missing_e3_report = root / "hardware-tests-without-e3.xml"
        missing_e3_report.write_text(
            test_report.read_text(encoding="utf-8").replace(e3_case, "case-0"),
            encoding="utf-8",
        )
        missing_e3_payload = missing_e3_report.read_bytes()
        missing_e3 = json.loads(json.dumps(valid))
        missing_e3["artifacts"]["hardware-tests.xml"].update(
            path=missing_e3_report.name,
            bytes=len(missing_e3_payload),
            sha256=hashlib.sha256(missing_e3_payload).hexdigest(),
        )
        try:
            verify_document(missing_e3, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: no progressive Vulkan viewer evidence")

        original_transcript = transcript.read_text(encoding="utf-8")
        transcript.write_text(
            original_transcript.replace('"point_star_count":100000',
                                        '"point_star_count":99999'),
            encoding="utf-8",
        )
        mismatched_payload = transcript.read_bytes()
        mismatched = json.loads(json.dumps(valid))
        mismatched["artifacts"]["transcript.log"].update(
            bytes=len(mismatched_payload),
            sha256=hashlib.sha256(mismatched_payload).hexdigest(),
        )
        try:
            verify_document(mismatched, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: scene claim/transcript mismatch")
        transcript.write_text(original_transcript, encoding="utf-8")

        viewer_transcript = root / "viewer-transcript.log"
        viewer_transcript.write_text(
            "window-created opengl-version=self-test\n"
            "frame-published backend=Vulkan width=64 height=64\n"
            "keyboard-callback key=87 action=press\n"
            "pointer-callback kind=cursor x=1 y=1 dragging=false\n",
            encoding="utf-8",
        )
        viewer_payload = viewer_transcript.read_bytes()
        viewer = json.loads(json.dumps(valid))
        viewer["domains"] = ["viewer-native-window-input"]
        viewer["claims"] = {
            "preset": "linux-ci",
            "test_estate_passed": True,
            "test_report": valid["claims"]["test_report"],
            "qualification_executable": qualification_claim,
            "runtime_ready": True,
            "viewer": {
                "backend": "Vulkan",
                "window_created": True,
                "frame_observed": True,
                "published_frame_count": 1,
                "keyboard_event_observed": True,
                "pointer_event_observed": True,
            },
        }
        viewer["artifacts"] = {
            "viewer-transcript.log": {
                "path": "viewer-transcript.log",
                "bytes": len(viewer_payload),
                "sha256": hashlib.sha256(viewer_payload).hexdigest(),
            },
            "hardware-tests.xml": valid["artifacts"]["hardware-tests.xml"],
            "device.json": valid["artifacts"]["device.json"],
            "alignment_receipt.json": valid["artifacts"]["alignment_receipt.json"],
            "mandatory_gate.json": valid["artifacts"]["mandatory_gate.json"],
            "ctest-inventory.json": valid["artifacts"]["ctest-inventory.json"],
            "qualification-sirius.bin": valid["artifacts"]["qualification-sirius.bin"],
            "qualification-gate-junit": valid["artifacts"]["qualification-gate-junit"],
            "qualification-gate-log": valid["artifacts"]["qualification-gate-log"],
            **{
                path.name: valid["artifacts"][path.name]
                for path in self_test_products.values()
            },
            "transcript.log": valid["artifacts"]["transcript.log"],
        }
        verify_document(viewer, root / "attestation.json")
        for field in (
            "window_created",
            "frame_observed",
            "keyboard_event_observed",
            "pointer_event_observed",
        ):
            candidate = json.loads(json.dumps(viewer))
            candidate["claims"]["viewer"][field] = False
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: viewer {field} false")

        for field, value in (("backend", "Cpu"), ("published_frame_count", 0)):
            candidate = json.loads(json.dumps(viewer))
            candidate["claims"]["viewer"][field] = value
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: viewer {field}={value!r}")

        viewer_transcript.write_text(
            "claimed callbacks without callback records\n", encoding="utf-8"
        )
        bogus_payload = viewer_transcript.read_bytes()
        bogus_viewer = json.loads(json.dumps(viewer))
        bogus_viewer["artifacts"]["viewer-transcript.log"].update(
            bytes=len(bogus_payload), sha256=hashlib.sha256(bogus_payload).hexdigest()
        )
        try:
            verify_document(bogus_viewer, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: viewer claims without callback transcript")

        # Native runtime domains use their own presets, but must meet the same
        # revision/readiness/JUnit/inventory/transcript chain as the Linux
        # hardware profiles. Positive controls prevent a fail-only verifier;
        # the paired removals prove metadata booleans cannot substitute for
        # hashed upstream evidence.
        native_runtime_controls = [
            (
                "windows-native-vulkan",
                "windows",
                "windows-msvc",
                {
                    **self_test_device,
                    "driver_name": "AMD Windows Vulkan",
                    "driver_info": "native self-test",
                },
            ),
            (
                "macos-moltenvk",
                "macos",
                "macos",
                {
                    **self_test_device,
                    "name": "Apple GPU",
                    "driver_name": "MoltenVK",
                    "driver_info": "native self-test",
                },
            ),
        ]
        for domain, platform_name, preset, native_device in native_runtime_controls:
            native_inventory = root / f"{domain}-device.json"
            native_inventory.write_text(
                json.dumps({
                    "platform": platform_name,
                    "wsl2": False,
                    "backends": {
                        "vulkan": {
                            "selected_device_index": 0,
                            "devices": [native_device],
                        }
                    },
                }),
                encoding="utf-8",
            )
            native_transcript = root / f"{domain}-transcript.log"
            native_transcript.write_text(
                f"== source revision: {source_revision}\n"
                f"== [1/6] configure + build ({preset})\n"
                "== readiness: evidence_generation=true vulkan=true selected_device_index=0\n"
                "100% tests passed, 0 tests failed out of "
                f"{MINIMUM_SOURCE_AVAILABLE_TESTS}\n",
                encoding="utf-8",
            )
            native_document = {
                "schema_version": 1,
                "status": "pass",
                "domains": [domain],
                "source_revision": source_revision,
                "completed_utc": "2026-07-31T00:00:00+00:00",
                "platform": {"os": platform_name, "wsl2": False},
                "device": native_device,
                "claims": {
                    "preset": preset,
                    "test_estate_passed": True,
                    "test_report": valid["claims"]["test_report"],
                    "qualification_executable": qualification_claim,
                    "runtime_ready": True,
                },
                "artifacts": {
                    "frame.png": valid["artifacts"]["frame.png"],
                    "hardware-tests.xml": valid["artifacts"]["hardware-tests.xml"],
                    "alignment_receipt.json": valid["artifacts"]["alignment_receipt.json"],
                    "mandatory_gate.json": valid["artifacts"]["mandatory_gate.json"],
                    "ctest-inventory.json": valid["artifacts"]["ctest-inventory.json"],
                    "qualification-sirius.bin": valid["artifacts"]["qualification-sirius.bin"],
                    "qualification-gate-junit": valid["artifacts"]["qualification-gate-junit"],
                    "qualification-gate-log": valid["artifacts"]["qualification-gate-log"],
                    **{
                        path.name: valid["artifacts"][path.name]
                        for path in self_test_products.values()
                    },
                    native_inventory.name: {
                        "path": native_inventory.name,
                        "bytes": native_inventory.stat().st_size,
                        "sha256": hashlib.sha256(native_inventory.read_bytes()).hexdigest(),
                    },
                    native_transcript.name: {
                        "path": native_transcript.name,
                        "bytes": native_transcript.stat().st_size,
                        "sha256": hashlib.sha256(native_transcript.read_bytes()).hexdigest(),
                    },
                },
            }
            verify_document(native_document, root / "attestation.json")
            for missing_artifact in (
                "hardware-tests.xml",
                "alignment_receipt.json",
                "mandatory_gate.json",
                "ctest-inventory.json",
                "qualification-sirius.bin",
                "qualification-gate-junit",
                "qualification-gate-log",
                *(path.name for path in self_test_products.values()),
                native_inventory.name,
                native_transcript.name,
            ):
                candidate = json.loads(json.dumps(native_document))
                candidate["artifacts"].pop(missing_artifact)
                try:
                    verify_document(candidate, root / "attestation.json")
                except ValueError:
                    continue
                raise ValueError(
                    f"negative control accepted: {domain} without {missing_artifact}"
                )

        corruptions = [
            ("software device", lambda doc: doc["device"].update(kind="software")),
            ("wrong dimensions", lambda doc: doc["claims"]["imax_render"].update(width=4096)),
            ("wrong 1080p dimensions",
             lambda doc: doc["claims"]["p3_render_1080p"].update(width=1919)),
            ("missing 1080p claim", lambda doc: doc["claims"].pop("p3_render_1080p")),
            ("missing 1080p artifact", lambda doc: doc["artifacts"].pop("frame-1080.png")),
            ("no ray bundles",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(ray_bundles=False)),
            ("texture star field",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(point_starfield=False)),
            ("disk-obscured point-star field",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(disk_enabled=True)),
            ("too few point stars",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(
                 star_catalogue_minimum=99999)),
            ("uncalibrated point stars",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(
                 point_brightness_scale=0.0)),
            ("stationary camera",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(
                 camera_beta=[0.0, 0.0, 0.0])),
            ("pinhole camera",
             lambda doc: doc["claims"]["imax_render"]["scene"].update(lens="Pinhole")),
            ("non-required hardware preset",
             lambda doc: doc["claims"].update(preset="linux-gcc")),
            ("mismatched hardware test summary",
             lambda doc: doc["claims"]["test_report"].update(
                 cases=MINIMUM_SOURCE_AVAILABLE_TESTS + 1)),
            ("missing hardware JUnit",
             lambda doc: doc["artifacts"].pop("hardware-tests.xml")),
            ("missing device inventory",
             lambda doc: doc["artifacts"].pop("device.json")),
            ("false test state", lambda doc: doc["claims"].update(test_estate_passed=False)),
            ("bad hash", lambda doc: doc["artifacts"]["frame.png"].update(sha256="0" * 64)),
            ("zero artifact size", lambda doc: doc["artifacts"]["frame.png"].update(bytes=0)),
            ("missing driver identity", lambda doc: doc["device"].update(driver_name="")),
            ("boolean vendor identity", lambda doc: doc["device"].update(vendor_id=True)),
            ("zero device memory", lambda doc: doc["device"].update(render_memory_bytes=0)),
            (
                "aliased artifact",
                lambda doc: doc["artifacts"].update(
                    alias={
                        "path": "frame.png",
                        "bytes": artifact.stat().st_size,
                        "sha256": digest,
                    }
                ),
            ),
            ("ambiguous source revision", lambda doc: doc.update(source_revision="main")),
            ("unbound source revision", lambda doc: doc.update(source_revision="1" * 40)),
            ("naive completion time", lambda doc: doc.update(completed_utc="2026-07-31T00:00:00")),
            ("non-UTC completion time",
             lambda doc: doc.update(completed_utc="2026-07-31T10:00:00+10:00")),
        ]
        for label, corrupt in corruptions:
            candidate = json.loads(json.dumps(valid))
            corrupt(candidate)
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: {label}")

        collapsed = root / "collapsed.png"
        write_self_test_png(collapsed, 5616, 4096, structured=False)
        collapsed_image = json.loads(json.dumps(valid))
        payload = collapsed.read_bytes()
        collapsed_image["artifacts"]["frame.png"].update(
            path="collapsed.png", bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest()
        )
        try:
            verify_document(collapsed_image, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: visually collapsed IMAX PNG")

        collapsed_1080 = root / "collapsed-1080.png"
        write_self_test_png(collapsed_1080, 1920, 1080, structured=False)
        collapsed_1080_image = json.loads(json.dumps(valid))
        payload_1080 = collapsed_1080.read_bytes()
        collapsed_1080_image["artifacts"]["frame-1080.png"].update(
            path=collapsed_1080.name,
            bytes=len(payload_1080),
            sha256=hashlib.sha256(payload_1080).hexdigest(),
        )
        try:
            verify_document(collapsed_1080_image, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: visually collapsed 1080p PNG")

        domain_controls = [
            (
                "Dozen outside WSL2",
                ["wsl2-dozen"],
                {"os": "linux", "wsl2": False},
                {"name": "AMD Radeon 780M", "kind": "integrated",
                 "driver_name": "Mesa Dozen", "driver_info": "dzn",
                 "vendor_id": 4098, "device_id": 5567, "api_version": 4206592,
                 "driver_id": 3, "device_local_bytes": 8589934592,
                 "render_memory_bytes": 2147483648, "supports_fp64": True},
            ),
            (
                "Dozen called native Windows",
                ["windows-native-vulkan"],
                {"os": "windows", "wsl2": False},
                {"name": "AMD Radeon 780M", "kind": "integrated",
                 "driver_name": "Mesa Dozen", "driver_info": "dzn",
                 "vendor_id": 4098, "device_id": 5567, "api_version": 4206592,
                 "driver_id": 3, "device_local_bytes": 8589934592,
                 "render_memory_bytes": 2147483648, "supports_fp64": True},
            ),
            (
                "non-MoltenVK macOS runtime",
                ["macos-moltenvk"],
                {"os": "macos", "wsl2": False},
                {"name": "Apple GPU", "kind": "integrated",
                 "driver_name": "unknown", "driver_info": "unknown",
                 "vendor_id": 4203, "device_id": 1, "api_version": 4206592,
                 "driver_id": 14, "device_local_bytes": 8589934592,
                 "render_memory_bytes": 2147483648, "supports_fp64": False},
            ),
        ]
        for label, domains, platform, device in domain_controls:
            candidate = json.loads(json.dumps(valid))
            candidate["domains"] = domains
            candidate["platform"] = platform
            candidate["device"] = device
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: {label}")

        native_build_report = root / "native-build-tests.xml"
        native_build_names = [alignment_case] + [
            f"native-build-case-{index}"
            for index in range(1, MINIMUM_SOURCE_AVAILABLE_TESTS)
        ]
        native_build_report.write_text(
            "<testsuite>"
            + "".join(f'<testcase name="{name}"/>' for name in native_build_names)
            + "</testsuite>\n",
            encoding="utf-8",
        )
        native_build_inventory = root / "native-build-ctest-inventory.json"
        native_build_inventory.write_text(
            json.dumps({
                "kind": "ctestInfo",
                "version": {"major": 1, "minor": 0},
                "tests": [
                    {"name": name, "properties": []} for name in native_build_names
                ],
            }),
            encoding="utf-8",
        )
        native_build_gate = root / "native-build-mandatory-gate.json"
        write_qualification_gate(
            native_build_gate,
            native_build_report,
            native_build_inventory,
            qualification_executable,
        )
        def self_test_artifact(path):
            payload = path.read_bytes()
            return {
                "path": path.name,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }

        native_build_document = {
            "schema_version": 1,
            "status": "pass",
            "domains": ["windows-native-build"],
            "source_revision": source_revision,
            "completed_utc": "2026-08-27T00:00:00+00:00",
            "platform": {"os": "windows", "wsl2": False},
            "device": {},
            "claims": {
                "test_estate_passed": True,
                "runtime_ready": False,
                "qualification_executable": qualification_claim,
                "test_report": {
                    "cases": MINIMUM_SOURCE_AVAILABLE_TESTS,
                    "failures": 0,
                    "errors": 0,
                    "skipped": 0,
                },
            },
            "artifacts": {
                native_build_report.name: self_test_artifact(native_build_report),
                native_build_receipt.name: self_test_artifact(native_build_receipt),
                native_build_inventory.name: self_test_artifact(native_build_inventory),
                native_build_gate.name: self_test_artifact(native_build_gate),
                qualification_executable.name: self_test_artifact(
                    qualification_executable
                ),
                "qualification-gate-junit": self_test_artifact(
                    root / "qualification-gate-junit"
                ),
                "qualification-gate-log": self_test_artifact(
                    root / "qualification-gate-log"
                ),
                **{
                    path.name: self_test_artifact(path)
                    for path in self_test_products.values()
                },
            },
        }
        verify_document(native_build_document, root / "attestation.json")
        for missing_artifact in (
            native_build_receipt.name,
            native_build_inventory.name,
            native_build_gate.name,
            qualification_executable.name,
            "qualification-gate-junit",
            "qualification-gate-log",
            *(path.name for path in self_test_products.values()),
        ):
            incomplete_build = json.loads(json.dumps(native_build_document))
            incomplete_build["artifacts"].pop(missing_artifact)
            try:
                verify_document(incomplete_build, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(
                f"negative control accepted: native build without {missing_artifact}"
            )
        stale_receipt_document = json.loads(json.dumps(native_build_document))
        stale_receipt = json.loads(native_build_receipt.read_text(encoding="utf-8"))
        stale_receipt["source_revision"] = "f" * 40
        native_build_receipt.write_text(json.dumps(stale_receipt), encoding="utf-8")
        stale_receipt_document["artifacts"][native_build_receipt.name] = \
            self_test_artifact(native_build_receipt)
        try:
            verify_document(stale_receipt_document, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: stale native build receipt")

        junit = root / "tests.xml"
        junit.write_text(
            '<testsuite tests="1" failures="0"><testcase name="pass"/></testsuite>\n',
            encoding="utf-8",
        )
        require(inspect_junit(junit)["cases"] == 1, "valid JUnit control was not counted")
        junit.write_text(
            '<testsuite tests="1" failures="1"><testcase name="fail">'
            '<failure message="intentional"/></testcase></testsuite>\n',
            encoding="utf-8",
        )
        try:
            inspect_junit(junit)
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: failing JUnit build evidence")

        clean_revision = "a" * 40
        require(
            validate_build_source_identity(clean_revision, "", clean_revision)
            == clean_revision,
            "clean native-build source identity was not accepted",
        )
        for label, revision, status, expected in (
            ("dirty native-build source", clean_revision, " M source.cpp", clean_revision),
            ("short native-build revision", "abcd", "", "abcd"),
            ("caller/source revision mismatch", clean_revision, "", "b" * 40),
        ):
            try:
                validate_build_source_identity(revision, status, expected)
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: {label}")


def validate_build_source_identity(revision, status, expected_revision=None):
    require(re.fullmatch(r"[0-9a-f]{40,64}", revision) is not None,
            "build source revision is not a full lowercase hexadecimal revision")
    require(status == "", "native build attestation requires a clean source tree")
    if expected_revision is not None:
        require(expected_revision == revision,
                "expected source revision differs from the tested Git worktree")
    return revision


def inspect_build_source(source_root, expected_revision=None):
    require(source_root is not None and source_root.is_dir(),
            "record-build requires a Git source root")
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=source_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=source_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"could not inspect build source worktree: {error}") from error
    return validate_build_source_identity(revision, status, expected_revision)


def record_build(args):
    require(args.domain in {"windows-native-build", "macos-native-build"},
            "record-build accepts only native build domains")
    revision = inspect_build_source(args.source_root, args.source_revision)
    require(args.artifact.is_file(), f"build evidence is missing: {args.artifact}")
    require(args.alignment_receipt.is_file(),
            f"build alignment receipt is missing: {args.alignment_receipt}")
    require(args.mandatory_gate.is_file(),
            f"qualification build gate is missing: {args.mandatory_gate}")
    require(args.qualification_executable.is_file(),
            f"qualification executable is missing: {args.qualification_executable}")
    require(args.mandatory_gate_junit.is_file(),
            f"qualification gate JUnit is missing: {args.mandatory_gate_junit}")
    require(args.mandatory_gate_log.is_file(),
            f"qualification gate log is missing: {args.mandatory_gate_log}")
    require(args.test_dir.is_dir(), f"native build test directory is missing: {args.test_dir}")
    report = inspect_junit(args.artifact)
    require(report["cases"] >= MINIMUM_SOURCE_AVAILABLE_TESTS
            and report["skipped"] == 0,
            "native build JUnit evidence is not a non-skipping test estate")
    inspect_build_alignment_receipt(args.alignment_receipt, revision)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    bundled = args.output.parent / args.artifact.name
    if args.artifact.resolve() != bundled.resolve():
        shutil.copy2(args.artifact, bundled)
    bundled_receipt = args.output.parent / "alignment_receipt.json"
    if args.alignment_receipt.resolve() != bundled_receipt.resolve():
        shutil.copy2(args.alignment_receipt, bundled_receipt)
    bundled_gate = args.output.parent / "mandatory_gate.json"
    if args.mandatory_gate.resolve() != bundled_gate.resolve():
        shutil.copy2(args.mandatory_gate, bundled_gate)
    bundled_executable = args.output.parent / "qualification-sirius.bin"
    if args.qualification_executable.resolve() != bundled_executable.resolve():
        shutil.copy2(args.qualification_executable, bundled_executable)
    bundled_gate_junit = args.output.parent / "qualification-gate-junit"
    if args.mandatory_gate_junit.resolve() != bundled_gate_junit.resolve():
        shutil.copy2(args.mandatory_gate_junit, bundled_gate_junit)
    bundled_gate_log = args.output.parent / "qualification-gate-log"
    if args.mandatory_gate_log.resolve() != bundled_gate_log.resolve():
        shutil.copy2(args.mandatory_gate_log, bundled_gate_log)
    bundled_products = {}
    qualification_resource_root = args.qualification_executable.parent / "resources"
    for logical_name, relative_path in QUALIFICATION_RUNTIME_RESOURCE_PATHS.items():
        source = qualification_resource_root / relative_path
        require(source.is_file(), f"qualification product is missing: {source}")
        destination = args.output.parent / QUALIFICATION_PRODUCT_EVIDENCE[logical_name]
        shutil.copy2(source, destination)
        bundled_products[logical_name] = destination
    inventory_command = [
        "ctest",
        "--test-dir",
        str(args.test_dir),
        "--show-only=json-v1",
    ]
    if args.platform == "windows":
        inventory_command[3:3] = ["-C", "Release"]
    try:
        inventory_result = subprocess.run(
            inventory_command,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"could not inventory native CTest registration: {error}") from error
    bundled_inventory = args.output.parent / "ctest-inventory.json"
    bundled_inventory.write_text(inventory_result.stdout, encoding="utf-8")
    require(junit_case_names(bundled) == inspect_ctest_inventory(bundled_inventory),
            "native build JUnit does not equal current CTest registration")
    payload = bundled.read_bytes()
    receipt_payload = bundled_receipt.read_bytes()
    inventory_payload = bundled_inventory.read_bytes()
    gate_payload = bundled_gate.read_bytes()
    executable_payload = bundled_executable.read_bytes()
    gate_junit_payload = bundled_gate_junit.read_bytes()
    gate_log_payload = bundled_gate_log.read_bytes()
    product_payloads = {
        logical_name: path.read_bytes()
        for logical_name, path in bundled_products.items()
    }
    document = {
        "schema_version": 1,
        "status": "pass",
        "domains": [args.domain],
        "source_revision": revision,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "platform": {"os": args.platform, "wsl2": False},
        "device": {},
        "claims": {
            "test_estate_passed": True,
            "runtime_ready": False,
            "test_report": report,
            "qualification_executable": {
                "artifact": bundled_executable.name,
                "bytes": len(executable_payload),
                "sha256": hashlib.sha256(executable_payload).hexdigest(),
            },
        },
        "artifacts": {
            bundled.name: {
                "path": bundled.name,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            },
            bundled_receipt.name: {
                "path": bundled_receipt.name,
                "bytes": len(receipt_payload),
                "sha256": hashlib.sha256(receipt_payload).hexdigest(),
            },
            bundled_inventory.name: {
                "path": bundled_inventory.name,
                "bytes": len(inventory_payload),
                "sha256": hashlib.sha256(inventory_payload).hexdigest(),
            },
            bundled_gate.name: {
                "path": bundled_gate.name,
                "bytes": len(gate_payload),
                "sha256": hashlib.sha256(gate_payload).hexdigest(),
            },
            bundled_executable.name: {
                "path": bundled_executable.name,
                "bytes": len(executable_payload),
                "sha256": hashlib.sha256(executable_payload).hexdigest(),
            },
            bundled_gate_junit.name: {
                "path": bundled_gate_junit.name,
                "bytes": len(gate_junit_payload),
                "sha256": hashlib.sha256(gate_junit_payload).hexdigest(),
            },
            bundled_gate_log.name: {
                "path": bundled_gate_log.name,
                "bytes": len(gate_log_payload),
                "sha256": hashlib.sha256(gate_log_payload).hexdigest(),
            },
            **{
                path.name: {
                    "path": path.name,
                    "bytes": len(product_payloads[logical_name]),
                    "sha256": hashlib.sha256(
                        product_payloads[logical_name]
                    ).hexdigest(),
                }
                for logical_name, path in bundled_products.items()
            },
        },
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    verify_document(document, args.output)
    require(inspect_build_source(args.source_root, args.source_revision) == revision,
            "source identity changed while issuing native build evidence")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("attestation", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--record-build", action="store_true")
    parser.add_argument("--domain")
    parser.add_argument("--platform")
    parser.add_argument("--source-revision")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--alignment-receipt", type=Path)
    parser.add_argument("--mandatory-gate", type=Path)
    parser.add_argument("--qualification-executable", type=Path)
    parser.add_argument("--mandatory-gate-junit", type=Path)
    parser.add_argument("--mandatory-gate-log", type=Path)
    parser.add_argument("--test-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("attestation verifier rejected every false-evidence control")
        elif args.record_build:
            require(args.platform in {"windows", "macos"},
                    "record-build platform must be windows or macos")
            require(args.artifact is not None and args.alignment_receipt is not None
                    and args.mandatory_gate is not None
                    and args.qualification_executable is not None
                    and args.mandatory_gate_junit is not None
                    and args.mandatory_gate_log is not None
                    and args.test_dir is not None and args.output is not None,
                    "record-build requires --artifact, --alignment-receipt, "
                    "--mandatory-gate, --mandatory-gate-junit, --mandatory-gate-log, "
                    "--qualification-executable, --test-dir, and --output")
            require(args.source_root is not None,
                    "record-build requires --source-root")
            record_build(args)
            print(f"build attestation recorded and verified: {args.output}")
        elif args.attestation is not None:
            verify_path(args.attestation)
            print(f"attestation verified: {args.attestation}")
        else:
            parser.error("provide an attestation path or --self-test")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.exit(1, f"attestation: {error}\n")


if __name__ == "__main__":
    main()
