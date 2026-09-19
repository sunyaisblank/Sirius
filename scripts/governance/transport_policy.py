"""Single-authority wiring for camera, disk and volume transport."""

from __future__ import annotations

import re
from pathlib import Path

from .common import (
    CPP_NON_CODE,
    SOURCE_ROOT,
    relative,
)


KERR_TRANSFER_AUTHORITY = SOURCE_ROOT / "core" / "relativistic_transfer.h"
KERR_TRANSFER_CPU_CONSUMER = SOURCE_ROOT / "backend" / "cpu" / "geodesic_tracer.cpp"
KERR_TRANSFER_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_disk.slang"
KERR_TRANSFER_DEVICE_CONSUMER = SOURCE_ROOT / "kernels" / "trace.slang"
KERR_TRANSFER_PARITY_PROBE = SOURCE_ROOT / "kernels" / "parity_probe.slang"
PAGE_THORNE_HOST_AUTHORITY = SOURCE_ROOT / "core" / "disk" / "novikov_thorne_disk.h"
PAGE_THORNE_CPU_CONSUMER = SOURCE_ROOT / "backend" / "cpu" / "geodesic_tracer.h"
PAGE_THORNE_DEVICE_AUTHORITY = KERR_TRANSFER_DEVICE_AUTHORITY
PAGE_THORNE_DEVICE_CONSUMER = KERR_TRANSFER_DEVICE_CONSUMER
PAGE_THORNE_PARITY_PROBE = KERR_TRANSFER_PARITY_PROBE
THIN_LENS_HOST_AUTHORITY = SOURCE_ROOT / "core" / "camera.h"
THIN_LENS_CPU_CONSUMER = SOURCE_ROOT / "backend" / "cpu" / "geodesic_tracer.cpp"
THIN_LENS_CPU_LAUNCH = SOURCE_ROOT / "core" / "camera_launch.h"
THIN_LENS_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_camera.slang"
THIN_LENS_DEVICE_CONSUMER = SOURCE_ROOT / "kernels" / "trace.slang"
THIN_LENS_PARITY_PROBE = SOURCE_ROOT / "kernels" / "parity_probe.slang"
THIN_LENS_SAMPLE_AUTHORITY = SOURCE_ROOT / "core" / "camera_sampling.h"
THIN_LENS_APP_BOUNDARY = SOURCE_ROOT / "app" / "config" / "config_loader.cpp"
THIN_LENS_SESSION_BOUNDARY = SOURCE_ROOT / "render" / "session" / "render_session.cpp"
THIN_LENS_PIXEL_CONSUMER = SOURCE_ROOT / "render" / "session" / "pixel_shading.cpp"
THIN_LENS_VULKAN_BOUNDARY = SOURCE_ROOT / "render" / "vulkan_renderer.cpp"
THIN_LENS_DISPATCH_BOUNDARY = SOURCE_ROOT / "render" / "dispatch_governor.cpp"
VOLUME_TRANSFER_HOST_AUTHORITY = KERR_TRANSFER_AUTHORITY
VOLUME_TRANSFER_CPU_CONSUMER = KERR_TRANSFER_CPU_CONSUMER
VOLUME_TRANSFER_DEVICE_AUTHORITY = KERR_TRANSFER_DEVICE_AUTHORITY
VOLUME_TRANSFER_DEVICE_CONSUMER = KERR_TRANSFER_DEVICE_CONSUMER
VOLUME_TRANSFER_PARITY_PROBE = KERR_TRANSFER_PARITY_PROBE


def page_thorne_edge_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    required = {
        PAGE_THORNE_HOST_AUTHORITY: (
            "config_.r_inner",
            "r_inner_ = (config_.r_inner > 0)",
            "FullPageThorneFlux",
            "kInnerEdgeBuffer",
        ),
        PAGE_THORNE_CPU_CONSUMER: (
            "page_thorne_disk_->Temperature",
            "page_thorne_disk_->InnerRadius",
        ),
        PAGE_THORNE_DEVICE_AUTHORITY: (
            "FullPageThorneFluxShape",
            "FullPageThorneTemperature",
            "rInner",
            "innerInM",
        ),
        PAGE_THORNE_DEVICE_CONSUMER: (
            "VolumeDiskTemperature",
            "FullPageThorneTemperature",
            "innerRadius",
        ),
        PAGE_THORNE_PARITY_PROBE: (
            "FullPageThorneFluxShape(r / p1, innerRadius / p1, p2)",
            "FullPageThorneTemperature(innerT, r, innerRadius, p1, p2)",
        ),
    }
    code_by_path: dict[Path, str] = {}
    for path, markers in required.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"Page-Thorne participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        code_by_path[path] = code
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits Page-Thorne edge marker {marker}")

    cpu = code_by_path.get(PAGE_THORNE_CPU_CONSUMER, "")
    if re.search(
        r"disk_config\.r_inner\s*=\s*config_\.disk_inner\s*/\s*cached_m_\s*;",
        cpu,
    ) is None:
        errors.append("the CPU Page-Thorne profile is not constructed from the declared edge")
    if re.search(
        r"kTemperatureReferenceRadiusRatio\s*\*\s*"
        r"page_thorne_disk_->InnerRadius\s*\(\s*\)",
        cpu,
    ) is None:
        errors.append("the CPU Page-Thorne scale is not normalised at the declared edge")
    if re.search(
        r"void\s+SetConfig\s*\([^)]*\)\s*\{.*?config_\s*=\s*config\s*;.*?"
        r"page_thorne_disk_\.reset\s*\(\s*\)\s*;.*?\}",
        cpu,
        re.DOTALL,
    ) is None:
        errors.append("the CPU Page-Thorne cache can survive a changed declared edge")

    device = code_by_path.get(PAGE_THORNE_DEVICE_CONSUMER, "")
    for radius in ("cylindricalR", "cr"):
        if re.search(
            rf"FullPageThorneTemperature\s*\(\s*innerTemperature\s*,\s*"
            rf"{radius}\s*,\s*innerRadius\s*,\s*M\s*,\s*aStar\s*\)",
            device,
        ) is None:
            errors.append(
                f"the Slang Page-Thorne {radius} path does not consume the declared edge"
            )
    if re.search(r"\brIsco\b|\bComputeKerrISCO\s*\(", device):
        errors.append("the Slang trace path recomputes ISCO instead of consuming its declared edge")
    return errors


def verify_page_thorne_edge_authority_policy() -> None:
    valid = {
        PAGE_THORNE_HOST_AUTHORITY: (
            "config_.r_inner r_inner_ = (config_.r_inner > 0) "
            "FullPageThorneFlux kInnerEdgeBuffer"
        ),
        PAGE_THORNE_CPU_CONSUMER: (
            "void SetConfig(const TracerConfig& config) { config_ = config; "
            "page_thorne_disk_.reset(); } "
            "disk_config.r_inner = config_.disk_inner / cached_m_; "
            "kTemperatureReferenceRadiusRatio * page_thorne_disk_->InnerRadius(); "
            "page_thorne_disk_->Temperature"
        ),
        PAGE_THORNE_DEVICE_AUTHORITY: (
            "FullPageThorneFluxShape FullPageThorneTemperature rInner innerInM"
        ),
        PAGE_THORNE_DEVICE_CONSUMER: (
            "VolumeDiskTemperature "
            "FullPageThorneTemperature(innerTemperature, cylindricalR, innerRadius, M, aStar); "
            "FullPageThorneTemperature(innerTemperature, cr, innerRadius, M, aStar);"
        ),
        PAGE_THORNE_PARITY_PROBE: (
            "FullPageThorneFluxShape(r / p1, innerRadius / p1, p2); "
            "FullPageThorneTemperature(innerT, r, innerRadius, p1, p2);"
        ),
    }
    if page_thorne_edge_authority_errors(valid):
        raise RuntimeError("Page-Thorne edge policy rejected the declared-edge wiring")

    isco_cpu = dict(valid)
    isco_cpu[PAGE_THORNE_CPU_CONSUMER] = isco_cpu[PAGE_THORNE_CPU_CONSUMER].replace(
        "config_.disk_inner / cached_m_", "disk.IscoRadius()"
    )
    if not page_thorne_edge_authority_errors(isco_cpu):
        raise RuntimeError("Page-Thorne edge policy accepted a CPU ISCO substitution")

    stale_cache = dict(valid)
    stale_cache[PAGE_THORNE_CPU_CONSUMER] = stale_cache[
        PAGE_THORNE_CPU_CONSUMER
    ].replace("page_thorne_disk_.reset();", "retain_cached_profile();")
    if not page_thorne_edge_authority_errors(stale_cache):
        raise RuntimeError("Page-Thorne edge policy accepted a stale CPU profile")

    isco_device = dict(valid)
    isco_device[PAGE_THORNE_DEVICE_CONSUMER] = isco_device[
        PAGE_THORNE_DEVICE_CONSUMER
    ].replace("innerRadius, M, aStar", "rIsco, M, aStar")
    if not page_thorne_edge_authority_errors(isco_device):
        raise RuntimeError("Page-Thorne edge policy accepted a Slang ISCO substitution")


def thin_lens_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    required = {
        THIN_LENS_HOST_AUTHORITY: (
            "ProjectThinLensSample",
            "ThinLensGeometryIssue",
            "ray.aperture_up = sample.pupil_up",
            "ray.aperture_right = sample.pupil_right",
            "config_.focus_distance, pupil_u, pupil_v",
        ),
        THIN_LENS_CPU_CONSUMER: (
            "LaunchCameraRay(observer_metric,",
            "outgoing_chart_->FromIngoing(launch->position)",
            "ray.position = launch->position",
            "ray.velocity = launch->tangent",
        ),
        THIN_LENS_CPU_LAUNCH: (
            "LaunchCameraRay", "camera.aperture_up", "camera.aperture_right",
            "central->spatial[1]", "central->spatial[2]",
            "frame_at(result.position)", "PastDirectedCameraRay",
        ),
        THIN_LENS_DEVICE_AUTHORITY: (
            "ThinLensProjectionSample",
            "ProjectThinLensSample",
            "pupilRight",
            "pupilUp",
        ),
        THIN_LENS_DEVICE_CONSUMER: (
            "ProjectThinLensSample",
            "ScaleVec4Cart(cameraFrame.up, Real(pupilUp))",
            "ScaleVec4Cart(cameraFrame.right, Real(pupilRight))",
            "state.x = pos0",
            "float pupilSampleU = params[66]",
            "float pupilSampleV = params[67]",
        ),
        THIN_LENS_PARITY_PROBE: (
            "OP_THIN_LENS_PROJECTION",
            "ProjectThinLensSample",
            "pupil.pupilRight",
            "pupil.pupilUp",
        ),
        THIN_LENS_SAMPLE_AUTHORITY: (
            "ForEachCameraSample",
            "RadicalInverse(ordinal, 5)",
            "RadicalInverse(ordinal, 7)",
        ),
        THIN_LENS_APP_BOUNDARY: ("ThinLensGeometryIssue",),
        THIN_LENS_SESSION_BOUNDARY: ("ThinLensGeometryIssue",),
        THIN_LENS_PIXEL_CONSUMER: (
            "ForEachCameraSample",
            "sample.image_u",
            "sample.image_v",
            "sample.pupil_u",
            "sample.pupil_v",
        ),
        THIN_LENS_VULKAN_BOUNDARY: (
            "ExecuteDispatchRegions",
            "config.samples_per_pixel",
            "params[44] = sample.image_u",
            "params[45] = sample.image_v",
            "params[66] = sample.pupil_u",
            "params[67] = sample.pupil_v",
        ),
        THIN_LENS_DISPATCH_BOUNDARY: (
            "ExecuteDispatchRegions",
            "ForEachCameraSample(samples_per_pixel",
            "submit(region, sample, sample_index)",
        ),
    }
    code_by_path: dict[Path, str] = {}
    for path, markers in required.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"thin-lens participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        code_by_path[path] = code
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits finite-pupil marker {marker}")

    host = code_by_path.get(THIN_LENS_HOST_AUTHORITY, "")
    if "pupil_radius > kMaximumThinLensPupilFraction * local_scale" not in host:
        errors.append("the host thin-lens authority omits its local tangent-plane bound")

    cpu = code_by_path.get(THIN_LENS_CPU_CONSUMER, "")
    # Follow the actual launch function, including its outgoing-chart change.
    # A store in an early failure branch, or matching dead code elsewhere in
    # the translation unit, cannot establish the successful launch dataflow.
    launch = ""
    signature = re.search(
        r"\bLightray\s+GeodesicTracer::InitializeLightray\s*\([^{}]*\)\s*\{",
        cpu,
        re.DOTALL,
    )
    if signature:
        depth = 1
        for index in range(signature.end(), len(cpu)):
            depth += (cpu[index] == "{") - (cpu[index] == "}")
            if depth == 0:
                launch = cpu[signature.end():index]
                break
    launch_flow = re.search(
        r"auto\s+launch\s*=\s*(?:step_executor_\s*\?\s*"
        r"step_executor_->Launch\(observer_metric,\s*cached_a_\s*\*\s*cached_m_,\s*camera_ray\)\s*:\s*)?"
        r"LaunchCameraRay\(observer_metric,\s*"
        r"cached_a_\s*\*\s*cached_m_,\s*camera_ray\)\s*;"
        r".*?const\s+auto\s+mapping\s*=\s*"
        r"outgoing_chart_->FromIngoing\(launch->position\)\s*;"
        r".*?launch->position\s*=\s*mapping->position\s*;"
        r"\s*launch->tangent\s*=\s*mapping->Apply\(launch->tangent\)\s*;"
        r"\s*launch->observer.time\s*=\s*mapping->Apply\(launch->observer.time\)\s*;"
        r"\s*for\s*\(auto&\s+axis\s*:\s*launch->observer.spatial\)"
        r"\s*axis\s*=\s*mapping->Apply\(axis\)\s*;"
        r"(?P<after_mapping>.*?)ray.position\s*=\s*launch->position\s*;"
        r"\s*ray.velocity\s*=\s*launch->tangent\s*;",
        launch, re.DOTALL,
    )
    if launch_flow is None:
        errors.append("the CPU tracer does not commit the displaced launch and chart transform")
    elif re.search(r"launch->(?:position|tangent|observer)\s*=(?!=)",
                   launch_flow.group("after_mapping")):
        errors.append("the CPU tracer discards the displaced or mapped pupil state")

    helper = code_by_path.get(THIN_LENS_CPU_LAUNCH, "")
    pupil_flow = re.search(
        r"\bLaunchCameraRay\s*\([^{}]*\)\s*\{"
        r".*?const\s+auto\s+central\s*=\s*frame_at\(result.position\)\s*;"
        r".*?result.position\s*\+=\s*central->spatial\[2\]\s*\*\s*camera.aperture_right"
        r"\s*-\s*central->spatial\[1\]\s*\*\s*camera.aperture_up\s*;"
        r"(?P<before_rebuild>.*?)const\s+auto\s+displaced\s*=\s*frame_at\(result.position\)\s*;"
        r"(?P<before_frame>.*?)result.observer\s*=\s*\*displaced\s*;"
        r"(?P<before_tangent>.*?)const\s+auto\s+tangent\s*=\s*"
        r"relativity::PastDirectedCameraRay\(result.observer,\s*direction\)\s*;"
        r"(?P<before_commit>.*?)result.tangent\s*=\s*\*tangent\s*;",
        helper, re.DOTALL,
    )
    if pupil_flow is None:
        errors.append("the CPU camera authority does not rebuild the displaced observer frame")
    else:
        gaps = " ".join(pupil_flow.group(name) for name in (
            "before_rebuild", "before_frame", "before_tangent", "before_commit"))
        if re.search(r"result\.(?:position|observer)(?:\s*\([^)]*\))?\s*=(?!=)", gaps):
            errors.append("the CPU camera authority discards the displaced pupil state")

    device = code_by_path.get(THIN_LENS_DEVICE_CONSUMER, "")
    if re.search(
        r"ProjectThinLensSample\s*\([^;]*?pupilSampleU\s*,\s*pupilSampleV\s*\)",
        device,
        re.DOTALL,
    ) is None:
        errors.append("the Slang thin lens reuses film coordinates for its pupil")
    if re.search(
        r"pos0\s*=\s*AddVec4Cart\s*\(.*?cameraFrame\.up.*?pupilUp.*?"
        r"cameraFrame\.right.*?pupilRight.*?state\.x\s*=\s*pos0",
        device,
        re.DOTALL,
    ) is None:
        errors.append("the Slang tracer does not move the live launch event across the pupil")
    return errors


def verify_thin_lens_authority_policy() -> None:
    valid = {
        THIN_LENS_HOST_AUTHORITY: (
            "ProjectThinLensSample ThinLensGeometryIssue "
            "ray.aperture_up = sample.pupil_up; "
            "ray.aperture_right = sample.pupil_right; "
            "config_.focus_distance, pupil_u, pupil_v; "
            "pupil_radius > kMaximumThinLensPupilFraction * local_scale"
        ),
        THIN_LENS_CPU_CONSUMER: (
            "Lightray GeodesicTracer::InitializeLightray(const CameraRay& camera_ray) { "
            "auto launch = LaunchCameraRay(observer_metric, cached_a_ * cached_m_, camera_ray); "
            "if (!launch) return ray; if (outgoing_chart_) { "
            "const auto mapping = outgoing_chart_->FromIngoing(launch->position); "
            "if (!mapping) return ray; launch->position = mapping->position; "
            "launch->tangent = mapping->Apply(launch->tangent); "
            "launch->observer.time = mapping->Apply(launch->observer.time); "
            "for (auto& axis : launch->observer.spatial) axis = mapping->Apply(axis); } "
            "ray.position = launch->position; ray.velocity = launch->tangent; return ray; }"
        ),
        THIN_LENS_CPU_LAUNCH: (
            "optional<CameraLaunch> LaunchCameraRay(IMetric& metric, double spin, const CameraRay& camera) { "
            "const auto central = frame_at(result.position); "
            "result.position += central->spatial[2] * camera.aperture_right - "
            "central->spatial[1] * camera.aperture_up; "
            "const auto displaced = frame_at(result.position); result.observer = *displaced; "
            "const auto tangent = relativity::PastDirectedCameraRay(result.observer, direction); "
            "result.tangent = *tangent; return result; }"
        ),
        THIN_LENS_DEVICE_AUTHORITY: (
            "ThinLensProjectionSample ProjectThinLensSample pupilRight pupilUp"
        ),
        THIN_LENS_DEVICE_CONSUMER: (
            "float pupilSampleU = params[66]; float pupilSampleV = params[67]; "
            "ProjectThinLensSample(imageX, imageY, tanHalfFov, focalLength, aperture, "
            "focusDistance, pupilSampleU, pupilSampleV); pos0 = AddVec4Cart("
            "ScaleVec4Cart(cameraFrame.up, Real(pupilUp)), "
            "ScaleVec4Cart(cameraFrame.right, Real(pupilRight))); state.x = pos0;"
        ),
        THIN_LENS_PARITY_PROBE: (
            "OP_THIN_LENS_PROJECTION ProjectThinLensSample "
            "pupil.pupilRight pupil.pupilUp"
        ),
        THIN_LENS_SAMPLE_AUTHORITY: (
            "ForEachCameraSample RadicalInverse(ordinal, 5) RadicalInverse(ordinal, 7)"
        ),
        THIN_LENS_APP_BOUNDARY: "ThinLensGeometryIssue",
        THIN_LENS_SESSION_BOUNDARY: "ThinLensGeometryIssue",
        THIN_LENS_PIXEL_CONSUMER: (
            "ForEachCameraSample sample.image_u sample.image_v "
            "sample.pupil_u sample.pupil_v"
        ),
        THIN_LENS_VULKAN_BOUNDARY: (
            "ExecuteDispatchRegions config.samples_per_pixel params[44] = sample.image_u; "
            "params[45] = sample.image_v; params[66] = sample.pupil_u; "
            "params[67] = sample.pupil_v;"
        ),
        THIN_LENS_DISPATCH_BOUNDARY: (
            "ExecuteDispatchRegions ForEachCameraSample(samples_per_pixel, callback); "
            "submit(region, sample, sample_index);"
        ),
    }
    if thin_lens_authority_errors(valid):
        raise RuntimeError("thin-lens policy rejected the finite-pupil wiring")

    for path, original, replacement in (
        (THIN_LENS_DISPATCH_BOUNDARY, "ForEachCameraSample", "ForEachImageSample"),
        (THIN_LENS_DISPATCH_BOUNDARY, "submit(region, sample, sample_index)",
         "submit(region, CameraSample{}, sample_index)"),
        (THIN_LENS_VULKAN_BOUNDARY, "ExecuteDispatchRegions", "SubmitUnsampledRegion"),
        (THIN_LENS_VULKAN_BOUNDARY, "params[66] = sample.pupil_u",
         "params[66] = sample.image_u"),
    ):
        changed = dict(valid)
        changed[path] = changed[path].replace(original, replacement)
        if not thin_lens_authority_errors(changed):
            raise RuntimeError(f"thin-lens policy accepted broken dispatch sampling: {original}")

    cpu_mutations = (
        (THIN_LENS_CPU_LAUNCH, "result.position +=", "discarded_position +="),
        (THIN_LENS_CPU_LAUNCH, "const auto displaced = frame_at(result.position)",
         "const auto displaced = central"),
        (THIN_LENS_CPU_LAUNCH, "result.observer = *displaced", "result.observer = *central"),
        (THIN_LENS_CPU_LAUNCH, "result.tangent = *tangent", "discarded_tangent = *tangent"),
        (THIN_LENS_CPU_LAUNCH, "const auto displaced =", "result.position = pinhole; const auto displaced ="),
        (THIN_LENS_CPU_CONSUMER, "FromIngoing(launch->position)", "FromIngoing(pinhole_position)"),
        (THIN_LENS_CPU_CONSUMER, "launch->position = mapping->position", "discarded_position = mapping->position"),
        (THIN_LENS_CPU_CONSUMER, "launch->tangent = mapping->Apply(launch->tangent)", "launch->tangent = pinhole_tangent"),
        (THIN_LENS_CPU_CONSUMER, "launch->observer.time = mapping->Apply(launch->observer.time)", "discarded_axis = mapping->Apply(launch->observer.time)"),
        (THIN_LENS_CPU_CONSUMER, "ray.position = launch->position", "ray.position = pinhole_position"),
        (THIN_LENS_CPU_CONSUMER, "ray.velocity = launch->tangent", "discarded_velocity = launch->tangent"),
        (THIN_LENS_CPU_CONSUMER, "GeodesicTracer::InitializeLightray", "GeodesicTracer::UnusedLaunchExample"),
    )
    for path, before, after in cpu_mutations:
        changed = dict(valid)
        assert before in changed[path]
        changed[path] = changed[path].replace(before, after)
        if not thin_lens_authority_errors(changed):
            raise RuntimeError(f"thin-lens policy accepted broken launch wiring: {before}")

    pinhole_device = dict(valid)
    pinhole_device[THIN_LENS_DEVICE_CONSUMER] = pinhole_device[
        THIN_LENS_DEVICE_CONSUMER
    ].replace("pos0 = AddVec4Cart(", "discardedPosition = AddVec4Cart(")
    if not thin_lens_authority_errors(pinhole_device):
        raise RuntimeError("thin-lens policy accepted a Slang direction-only pupil")

    unbounded = dict(valid)
    unbounded[THIN_LENS_HOST_AUTHORITY] = unbounded[THIN_LENS_HOST_AUTHORITY].replace(
        "pupil_radius > kMaximumThinLensPupilFraction * local_scale", "false"
    )
    if not thin_lens_authority_errors(unbounded):
        raise RuntimeError("thin-lens policy accepted an unbounded local pupil")

    diagonal_host = dict(valid)
    diagonal_host[THIN_LENS_PIXEL_CONSUMER] = diagonal_host[
        THIN_LENS_PIXEL_CONSUMER
    ].replace("sample.pupil_u sample.pupil_v", "sample.image_u sample.image_v")
    if not thin_lens_authority_errors(diagonal_host):
        raise RuntimeError("thin-lens policy accepted collapsed CPU film/pupil dimensions")

    diagonal_device = dict(valid)
    diagonal_device[THIN_LENS_DEVICE_CONSUMER] = diagonal_device[
        THIN_LENS_DEVICE_CONSUMER
    ].replace("pupilSampleU, pupilSampleV", "sampleU, sampleV")
    if not thin_lens_authority_errors(diagonal_device):
        raise RuntimeError("thin-lens policy accepted collapsed Slang film/pupil dimensions")


def kerr_zamo_transfer_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    required = {
        KERR_TRANSFER_AUTHORITY: (
            "TryKerrStationaryFrameFrequencyTransfer",
            "KerrZamoFrequencyTransfer",
            "KerrDiskTransfer",
        ),
        KERR_TRANSFER_CPU_CONSUMER: (
            "relativity::TryKerrStationaryFrameFrequencyTransfer",
            "relativity::KerrZamoFrequencyTransfer",
            "zamo_transfer->frame_frequency",
        ),
        KERR_TRANSFER_DEVICE_AUTHORITY: (
            "TryKerrStationaryFrameFrequencyTransfer",
            "KerrZamoFrequencyTransfer",
            "KerrPhotonKillingQuantitiesCart",
        ),
        KERR_TRANSFER_DEVICE_CONSUMER: (
            "TryKerrStationaryFrameFrequencyTransfer",
            "KerrZamoFrequencyTransfer",
            "zamoTransfer.y",
        ),
        KERR_TRANSFER_PARITY_PROBE: (
            "OP_KERR_ZAMO_TRANSFER",
            "OP_KERR_DISK_TRANSFER",
            "KerrZamoFrequencyTransfer",
            "ComputeKerrDiskTransferCart",
        ),
    }
    for path, markers in required.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"Kerr ZAMO transfer participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits Kerr ZAMO marker {marker}")

    cpu = CPP_NON_CODE.sub(" ", documents.get(KERR_TRANSFER_CPU_CONSUMER, ""))
    if "eulerian_frequency" in cpu:
        errors.append("the CPU volume path restored a Kerr-Schild Eulerian frequency substitute")
    device = CPP_NON_CODE.sub(" ", documents.get(KERR_TRANSFER_DEVICE_CONSUMER, ""))
    if "eulerianFrequency" in device:
        errors.append("the Slang volume path restored a Kerr-Schild Eulerian frequency substitute")
    return errors


def verify_kerr_zamo_transfer_authority_policy() -> None:
    valid = {
        KERR_TRANSFER_AUTHORITY: (
            "TryKerrStationaryFrameFrequencyTransfer KerrZamoFrequencyTransfer "
            "KerrDiskTransfer"
        ),
        KERR_TRANSFER_CPU_CONSUMER: (
            "relativity::TryKerrStationaryFrameFrequencyTransfer "
            "relativity::KerrZamoFrequencyTransfer zamo_transfer->frame_frequency"
        ),
        KERR_TRANSFER_DEVICE_AUTHORITY: (
            "TryKerrStationaryFrameFrequencyTransfer KerrZamoFrequencyTransfer "
            "KerrPhotonKillingQuantitiesCart"
        ),
        KERR_TRANSFER_DEVICE_CONSUMER: (
            "TryKerrStationaryFrameFrequencyTransfer KerrZamoFrequencyTransfer zamoTransfer.y"
        ),
        KERR_TRANSFER_PARITY_PROBE: (
            "OP_KERR_ZAMO_TRANSFER OP_KERR_DISK_TRANSFER "
            "KerrZamoFrequencyTransfer ComputeKerrDiskTransferCart"
        ),
    }
    if kerr_zamo_transfer_authority_errors(valid):
        raise RuntimeError("Kerr ZAMO policy rejected the governed host/device authority")
    detached = dict(valid)
    detached[KERR_TRANSFER_CPU_CONSUMER] = "independent_volume_frequency();"
    if not kerr_zamo_transfer_authority_errors(detached):
        raise RuntimeError("Kerr ZAMO policy accepted a detached CPU volume consumer")
    slicing_normal = dict(valid)
    slicing_normal[KERR_TRANSFER_DEVICE_CONSUMER] += " eulerianFrequency"
    if not kerr_zamo_transfer_authority_errors(slicing_normal):
        raise RuntimeError("Kerr ZAMO policy accepted a Kerr-Schild slicing-normal substitute")


def volumetric_transfer_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    required = {
        VOLUME_TRANSFER_HOST_AUTHORITY: (
            "GreyLayerAbsorbedFraction",
            "std::expm1",
            "AccumulateObserverToSourceLayer",
        ),
        VOLUME_TRANSFER_CPU_CONSUMER: (
            "AccumulateVolumetricEmission",
            "transfer.optical_depth > 0.0",
        ),
        VOLUME_TRANSFER_DEVICE_AUTHORITY: ("GreyLayerAbsorbedFraction",),
        VOLUME_TRANSFER_DEVICE_CONSUMER: (
            "AccumulateVolumeSegment",
            "GreyLayerAbsorbedFraction",
            "volumeOpticalDepth > 0.0f",
        ),
        VOLUME_TRANSFER_PARITY_PROBE: (
            "OP_GREY_LAYER_ABSORPTION",
            "GreyLayerAbsorbedFraction",
        ),
    }
    for path, markers in required.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"volumetric transfer participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits volumetric transfer marker {marker}")

    cpu = CPP_NON_CODE.sub(" ", documents.get(VOLUME_TRANSFER_CPU_CONSUMER, ""))
    first_accumulation = cpu.find("AccumulateVolumetricEmission(")
    if first_accumulation < 0:
        errors.append("the CPU path has no accepted-segment volume accumulator call")
    elif "IsInVolumetricDisk" in cpu[:first_accumulation]:
        errors.append(
            "the CPU volume path prefilters accepted segments by endpoint membership"
        )
    if re.search(r"optical_depth\s*>\s*0\.01", cpu):
        errors.append("the CPU volume path discards represented optically thin transfer")

    device = CPP_NON_CODE.sub(" ", documents.get(VOLUME_TRANSFER_DEVICE_CONSUMER, ""))
    if re.search(r"volumeOpticalDepth\s*>\s*0\.01", device):
        errors.append("the Slang volume path discards represented optically thin transfer")
    return errors


def verify_volumetric_transfer_authority_policy() -> None:
    valid = {
        VOLUME_TRANSFER_HOST_AUTHORITY: (
            "GreyLayerAbsorbedFraction std::expm1 AccumulateObserverToSourceLayer"
        ),
        VOLUME_TRANSFER_CPU_CONSUMER: (
            "AccumulateVolumetricEmission( transfer.optical_depth > 0.0"
        ),
        VOLUME_TRANSFER_DEVICE_AUTHORITY: "GreyLayerAbsorbedFraction",
        VOLUME_TRANSFER_DEVICE_CONSUMER: (
            "AccumulateVolumeSegment GreyLayerAbsorbedFraction "
            "volumeOpticalDepth > 0.0f"
        ),
        VOLUME_TRANSFER_PARITY_PROBE: (
            "OP_GREY_LAYER_ABSORPTION GreyLayerAbsorbedFraction"
        ),
    }
    if volumetric_transfer_authority_errors(valid):
        raise RuntimeError("volumetric-transfer policy rejected the complete transfer seam")

    endpoint_prefilter = dict(valid)
    endpoint_prefilter[VOLUME_TRANSFER_CPU_CONSUMER] = (
        "if (IsInVolumetricDisk(previous_r, previous_z)) "
        "AccumulateVolumetricEmission( transfer.optical_depth > 0.0"
    )
    if not volumetric_transfer_authority_errors(endpoint_prefilter):
        raise RuntimeError("volumetric-transfer policy accepted an endpoint prefilter")

    thin_cutoff = dict(valid)
    thin_cutoff[VOLUME_TRANSFER_CPU_CONSUMER] = (
        "AccumulateVolumetricEmission( transfer.optical_depth > 0.01"
    )
    thin_cutoff[VOLUME_TRANSFER_DEVICE_CONSUMER] = (
        "AccumulateVolumeSegment GreyLayerAbsorbedFraction "
        "volumeOpticalDepth > 0.01f"
    )
    if len(volumetric_transfer_authority_errors(thin_cutoff)) < 2:
        raise RuntimeError("volumetric-transfer policy accepted an optically thin cutoff")
