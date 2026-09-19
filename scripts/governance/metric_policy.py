"""Live metric, curvature and matter-equation source contracts."""

from __future__ import annotations

import re
from pathlib import Path

from .common import (
    CPP_NON_CODE,
    OPERATING_MODEL,
    ROOT,
    SOURCE_ROOT,
    relative,
)
from .transport_policy import (
    KERR_TRANSFER_CPU_CONSUMER,
    KERR_TRANSFER_DEVICE_CONSUMER,
    KERR_TRANSFER_PARITY_PROBE,
)


MORRIS_HOST_AUTHORITY = SOURCE_ROOT / "core" / "metrics" / "morris_thorne_family.h"
MORRIS_HOST_EVENT_AUTHORITY = SOURCE_ROOT / "core" / "trace_boundary.h"
MORRIS_CPU_CONSUMER = KERR_TRANSFER_CPU_CONSUMER
MORRIS_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_metrics.slang"
MORRIS_DEVICE_EVENT_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_trace_event.slang"
MORRIS_DEVICE_TOPOLOGY_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_ellis_topology.slang"
MORRIS_DEVICE_CONSUMER = KERR_TRANSFER_DEVICE_CONSUMER
MORRIS_PARITY_PROBE = KERR_TRANSFER_PARITY_PROBE
ALCUBIERRE_REGISTRY_AUTHORITY = SOURCE_ROOT / "core" / "metrics" / "registry.h"
ALCUBIERRE_HOST_AUTHORITY = SOURCE_ROOT / "core" / "metrics" / "warp_drive_family.h"
ALCUBIERRE_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_metrics.slang"
ALCUBIERRE_DEVICE_INTEGRATOR = SOURCE_ROOT / "kernels" / "gr_integrator.slang"
ALCUBIERRE_PARITY_PROBE = KERR_TRANSFER_PARITY_PROBE
KERR_SCHILD_HOST_AUTHORITY = SOURCE_ROOT / "core" / "metrics" / "kerr_schild_family.h"
KERR_SCHILD_CURVATURE_ORACLE = ROOT / "tests" / "core" / "metric_curvature_oracle.h"
KERR_SCHILD_FIELD_EQUATION_GATE = (
    ROOT / "tests" / "core" / "kerr_schild_field_equation_test.cpp"
)
KERR_SCHILD_TEST_CMAKE = ROOT / "tests" / "core" / "CMakeLists.txt"
KERR_SCHILD_OPERATING_MODEL = OPERATING_MODEL


def morris_thorne_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    required = {
        MORRIS_HOST_AUTHORITY: (
            "MorrisThorneCartesian",
            "params.Phi0 == 0.0",
            "represented_parameter",
            "g(0, 0) = Dual<double>(-1.0)",
            "IsotropicThroatRadius",
            "IsotropicEllisThroatRadius",
            "EllisInvertedIsotropicRadius",
            "MapEllisSecondSheetSkyDirection",
            "conformal_base",
        ),
        MORRIS_HOST_EVENT_AUTHORITY: (
            "FindPolynomialRootsOnUnitInterval",
            "SphericalSegmentPolynomial",
            "FindSphericalBoundaryEvent",
            "SphericalBoundarySense",
            "FindSphericalCaptureEvent",
        ),
        MORRIS_CPU_CONSUMER: (
            "IsotropicEllisThroatRadius",
            "WormholeTopology::TwoSheet",
            "FindSphericalCaptureEvent",
            "FindSphericalBoundaryEvent",
            "SphericalBoundarySense::DecreasingRadius",
            "terminal_throat_boundary",
            "terminal_opposite_infinity",
            "TraceResult::Outcome::Throat",
            "TraceResult::AsymptoticSheet::Opposite",
            "MapEllisSecondSheetSkyDirection",
        ),
        MORRIS_DEVICE_AUTHORITY: (
            "IsMorrisThorneCartesianEventRepresented",
            "Phi0 != 0.0f",
            "GetMorrisThorneCartesianMetric",
            "GetMorrisThorneCartesianChristoffel",
            "g[0][0] = -1.0f",
            "conformalBase",
        ),
        MORRIS_DEVICE_EVENT_AUTHORITY: (
            "IsFiniteAcceptedSegmentValue",
            "FindSexticRoots",
            "SphericalSegmentPolynomial",
            "FindSphericalBoundaryEvent",
            "kSphericalBoundaryDecreasingRadius",
            "FindSphericalCaptureEvent",
        ),
        MORRIS_DEVICE_TOPOLOGY_AUTHORITY: (
            "EllisOppositeEscapeRadius",
            "MapEllisSecondSheetSkyDirection",
            "TraceEllisTwoSheet",
            "InitEllisTwoSheet",
            "AdvanceEllisTwoSheet",
            "kEllisTraceOppositeInfinity",
            "FindSphericalBoundaryEvent",
        ),
        MORRIS_DEVICE_CONSUMER: (
            "FindSphericalCaptureEvent",
            "terminalCaptureSurface",
            "InitEllisTwoSheet",
            "AdvanceEllisTwoSheet",
            "oppositeSheet",
            "MapEllisSecondSheetSkyDirection",
        ),
        MORRIS_PARITY_PROBE: (
            "OP_SPHERICAL_CAPTURE_EVENT",
            "FindSphericalCaptureEvent",
            "OP_ELLIS_TWO_SHEET_TRACE",
            "TraceEllisTwoSheet",
            "MapEllisSecondSheetSkyDirection",
            "GetMorrisThorneCartesianMetric",
            "GetMorrisThorneCartesianChristoffel",
        ),
    }
    code_by_path: dict[Path, str] = {}
    for path, markers in required.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"Morris-Thorne authority participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        code_by_path[path] = code
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits Morris-Thorne marker {marker}")

    stale = re.compile(
        r"\b(?:GetMorrisThorneMetric|GetMorrisThorneChristoffel|"
        r"WormholeShapeFunction|WormholeShapeDerivative)\b|"
        r"b0\s*\*\s*1\.001|\bb\s*/\s*\(\s*r\s*-\s*b\s*\)"
    )
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        match = stale.search(code)
        if match:
            errors.append(
                f"{relative(path)} restores a clamped/areal Cartesian Morris-Thorne path: "
                f"{match.group(0)}"
            )
        if "one_sheet_topology" in code:
            errors.append(f"{relative(path)} restores the lossy wormhole topology boolean")

    host_metric = code_by_path.get(MORRIS_HOST_AUTHORITY, "")
    if "InsideCaptureSurface" in host_metric or "SphericalCaptureRadius" in host_metric:
        errors.append("the regular Ellis throat is exposed as an intrinsic capture surface")

    topology = code_by_path.get(MORRIS_DEVICE_TOPOLOGY_AUTHORITY, "")
    if topology.count("TraceEllisTwoSheet(") != 1:
        errors.append("the device Ellis topology module has no single trace authority")
    for name, return_type in (
        ("InitEllisTwoSheet", "EllisTopologyContinuation"),
        ("AdvanceEllisTwoSheet", "void"),
    ):
        definition = rf"\bpublic\s+{return_type}\s+{name}\s*\("
        if (
            len(re.findall(definition, topology)) != 1
            or len(re.findall(rf"\b{name}\s*\(", topology)) != 2
        ):
            errors.append(f"the device Ellis topology module has no single shared {name} authority")
        if len(re.findall(rf"\b{name}\s*\(", code_by_path.get(MORRIS_DEVICE_CONSUMER, ""))) != 1:
            errors.append(
                f"{relative(MORRIS_DEVICE_CONSUMER)} does not consume shared {name} exactly once"
            )
    wrapper = re.search(
        r"\bpublic\s+EllisTopologyTraceResult\s+TraceEllisTwoSheet\s*\([^)]*\)\s*"
        r"\{(?P<body>.*?)return\s+continuation\.trace\s*;\s*\}",
        topology,
        re.DOTALL,
    )
    if wrapper is None or any(
        len(re.findall(rf"\b{name}\s*\(", wrapper.group("body"))) != 1
        for name in ("InitEllisTwoSheet", "AdvanceEllisTwoSheet")
    ):
        errors.append("the device Ellis parity wrapper bypasses shared initialization/advancement")
    if code_by_path.get(MORRIS_PARITY_PROBE, "").count("TraceEllisTwoSheet(") != 1:
        errors.append("the device Ellis parity probe does not consume the shared trace wrapper")

    for path, consumers in (
        (MORRIS_CPU_CONSUMER, ("StepBundle(", "AdvancePolarisationFrame(",
                               "AccumulateVolumetricEmission(", "FindDiskIntersection(")),
        (MORRIS_DEVICE_CONSUMER, ("AccumulateVolumeSegment(",)),
    ):
        code = code_by_path.get(path, "")
        trace_start = code.find("TraceResult GeodesicTracer::Trace") if path == MORRIS_CPU_CONSUMER else 0
        capture = code.find("FindSphericalCaptureEvent(", max(trace_start, 0))
        if capture < 0:
            continue
        for consumer in consumers:
            consumer_call = code.find(consumer, capture)
            if consumer_call < 0:
                errors.append(
                    f"{relative(path)} does not expose {consumer} after throat localisation"
                )
    return errors


def verify_morris_thorne_authority_policy() -> None:
    valid = {
        MORRIS_HOST_AUTHORITY: (
            "MorrisThorneCartesian params.Phi0 == 0.0 represented_parameter "
            "g(0, 0) = Dual<double>(-1.0) "
            "IsotropicThroatRadius IsotropicEllisThroatRadius "
            "EllisInvertedIsotropicRadius MapEllisSecondSheetSkyDirection conformal_base"
        ),
        MORRIS_HOST_EVENT_AUTHORITY: (
            "FindPolynomialRootsOnUnitInterval SphericalSegmentPolynomial "
            "FindSphericalBoundaryEvent SphericalBoundarySense FindSphericalCaptureEvent"
        ),
        MORRIS_CPU_CONSUMER: (
            "TraceResult GeodesicTracer::Trace IsotropicEllisThroatRadius "
            "WormholeTopology::TwoSheet FindSphericalCaptureEvent "
            "FindSphericalBoundaryEvent SphericalBoundarySense::DecreasingRadius "
            "terminal_throat_boundary terminal_opposite_infinity TraceResult::Outcome::Throat "
            "TraceResult::AsymptoticSheet::Opposite MapEllisSecondSheetSkyDirection StepBundle( "
            "AdvancePolarisationFrame( AccumulateVolumetricEmission( FindDiskIntersection("
        ),
        MORRIS_DEVICE_AUTHORITY: (
            "IsMorrisThorneCartesianEventRepresented Phi0 != 0.0f "
            "GetMorrisThorneCartesianMetric GetMorrisThorneCartesianChristoffel "
            "g[0][0] = -1.0f conformalBase"
        ),
        MORRIS_DEVICE_EVENT_AUTHORITY: (
            "IsFiniteAcceptedSegmentValue FindSexticRoots SphericalSegmentPolynomial "
            "FindSphericalBoundaryEvent kSphericalBoundaryDecreasingRadius "
            "FindSphericalCaptureEvent"
        ),
        MORRIS_DEVICE_TOPOLOGY_AUTHORITY: (
            "EllisOppositeEscapeRadius MapEllisSecondSheetSkyDirection "
            "kEllisTraceOppositeInfinity FindSphericalBoundaryEvent "
            "public EllisTopologyContinuation InitEllisTwoSheet() {} "
            "public void AdvanceEllisTwoSheet() {} "
            "public EllisTopologyTraceResult TraceEllisTwoSheet() { "
            "continuation = InitEllisTwoSheet(); AdvanceEllisTwoSheet(); "
            "return continuation.trace; }"
        ),
        MORRIS_DEVICE_CONSUMER: (
            "FindSphericalCaptureEvent terminalCaptureSurface "
            "InitEllisTwoSheet( AdvanceEllisTwoSheet( "
            "oppositeSheet MapEllisSecondSheetSkyDirection AccumulateVolumeSegment("
        ),
        MORRIS_PARITY_PROBE: (
            "OP_SPHERICAL_CAPTURE_EVENT FindSphericalCaptureEvent "
            "OP_ELLIS_TWO_SHEET_TRACE TraceEllisTwoSheet( MapEllisSecondSheetSkyDirection "
            "GetMorrisThorneCartesianMetric GetMorrisThorneCartesianChristoffel"
        ),
    }
    if morris_thorne_authority_errors(valid):
        raise RuntimeError("Morris-Thorne policy rejected the exact event-synchronised seam")

    clamped = dict(valid)
    clamped[MORRIS_DEVICE_AUTHORITY] += " float r = max(pos.r, b0 * 1.001f);"
    if not morris_thorne_authority_errors(clamped):
        raise RuntimeError("Morris-Thorne policy accepted a fabricated throat clamp")

    endpoint_only = dict(valid)
    endpoint_only[MORRIS_CPU_CONSUMER] = endpoint_only[MORRIS_CPU_CONSUMER].replace(
        "FindSphericalCaptureEvent", "InsideCaptureSurface"
    )
    if not morris_thorne_authority_errors(endpoint_only):
        raise RuntimeError("Morris-Thorne policy accepted endpoint-only throat capture")

    intrinsic_capture = dict(valid)
    intrinsic_capture[MORRIS_HOST_AUTHORITY] += " InsideCaptureSurface"
    if not morris_thorne_authority_errors(intrinsic_capture):
        raise RuntimeError("Morris-Thorne policy accepted an intrinsic throat capture surface")

    unnormalised_lapse = dict(valid)
    unnormalised_lapse[MORRIS_HOST_AUTHORITY] = unnormalised_lapse[
        MORRIS_HOST_AUTHORITY
    ].replace("params.Phi0 == 0.0", "std::isfinite(params.Phi0)")
    unnormalised_lapse[MORRIS_DEVICE_AUTHORITY] = unnormalised_lapse[
        MORRIS_DEVICE_AUTHORITY
    ].replace("Phi0 != 0.0f", "!isfinite(Phi0)")
    if len(morris_thorne_authority_errors(unnormalised_lapse)) < 2:
        raise RuntimeError("Morris-Thorne policy accepted an unnormalised live lapse")

    independent_device_trace = dict(valid)
    independent_device_trace[MORRIS_PARITY_PROBE] = independent_device_trace[
        MORRIS_PARITY_PROBE
    ].replace("TraceEllisTwoSheet(", "IndependentEllisTrace(")
    if not morris_thorne_authority_errors(independent_device_trace):
        raise RuntimeError("Morris-Thorne policy accepted an independent parity trajectory")

    for name in ("InitEllisTwoSheet", "AdvanceEllisTwoSheet"):
        bypassed = dict(valid)
        bypassed[MORRIS_DEVICE_CONSUMER] = bypassed[MORRIS_DEVICE_CONSUMER].replace(
            name + "(", "Independent" + name + "("
        )
        if not morris_thorne_authority_errors(bypassed):
            raise RuntimeError(f"Morris-Thorne policy accepted production bypass of {name}")
        bypassed_wrapper = dict(valid)
        authority = bypassed_wrapper[MORRIS_DEVICE_TOPOLOGY_AUTHORITY]
        split = authority.index("public EllisTopologyTraceResult TraceEllisTwoSheet")
        bypassed_wrapper[MORRIS_DEVICE_TOPOLOGY_AUTHORITY] = (
            authority[:split] + authority[split:].replace(name + "(", "Independent" + name + "(")
        )
        if not morris_thorne_authority_errors(bypassed_wrapper):
            raise RuntimeError(f"Morris-Thorne policy accepted parity-wrapper bypass of {name}")


def alcubierre_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    required = {
        ALCUBIERRE_REGISTRY_AUTHORITY: (
            "AlcubierreScaleIssue",
            "AlcubierreParameterIssue",
            "kMinAlcubierreSigmaRadius",
            "kMaxAlcubierreSigmaRadius",
        ),
        ALCUBIERRE_HOST_AUTHORITY: (
            "AlcubierreParameterIssue(params.vs, params.R, params.sigma)",
            "AlcubierreParameterIssue(proposed.vs, proposed.R, proposed.sigma)",
            "WarpDriveFamily::InverseMetric",
            "g_inv(0, 0) = Dual<double>(-1.0)",
            "g_inv(1, 1) = Dual<double>(1.0 - vsf * vsf)",
        ),
        ALCUBIERRE_DEVICE_AUTHORITY: (
            "IsWarpDriveEventRepresented",
            "sigmaR >= 0.1f && sigmaR <= 100.0f",
            "GetWarpDriveMetric",
            "GetWarpDriveChristoffel",
            "g_inv[0][0] = -1.0f",
            "g_inv[1][1] = 1.0f - vsf2",
        ),
        ALCUBIERRE_DEVICE_INTEGRATOR: (
            "TryGeodesicAccelerationWarpCart",
            "IsWarpDriveEventRepresented(state.x",
            "IntegrateGeodesicRK4WarpCartDelta",
            "return IsWarpDriveEventRepresented(candidate",
        ),
        ALCUBIERRE_PARITY_PROBE: (
            "GetWarpDriveMetric",
            "GetWarpDriveChristoffel",
            "IntegrateGeodesicRK4WarpCart",
        ),
    }
    code_by_path: dict[Path, str] = {}
    for path, markers in required.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"Alcubierre authority participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        code_by_path[path] = code
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits Alcubierre marker {marker}")

    device = code_by_path.get(ALCUBIERRE_DEVICE_AUTHORITY, "")
    if "det_tx" in device:
        errors.append("the exact Alcubierre inverse was replaced by a repaired determinant")
    if device.count("IsWarpDriveEventRepresented(") != 3:
        errors.append("the device Alcubierre metric/connection do not share one domain predicate")
    integrator = code_by_path.get(ALCUBIERRE_DEVICE_INTEGRATOR, "")
    if integrator.count("TryGeodesicAccelerationWarpCart(") != 5:
        errors.append("not every Alcubierre RK4 stage consumes the fallible acceleration")
    return errors


def verify_alcubierre_authority_policy() -> None:
    valid = {
        ALCUBIERRE_REGISTRY_AUTHORITY: (
            "AlcubierreScaleIssue AlcubierreParameterIssue "
            "kMinAlcubierreSigmaRadius kMaxAlcubierreSigmaRadius"
        ),
        ALCUBIERRE_HOST_AUTHORITY: (
            "AlcubierreParameterIssue(params.vs, params.R, params.sigma) "
            "AlcubierreParameterIssue(proposed.vs, proposed.R, proposed.sigma) "
            "WarpDriveFamily::InverseMetric g_inv(0, 0) = Dual<double>(-1.0) "
            "g_inv(1, 1) = Dual<double>(1.0 - vsf * vsf)"
        ),
        ALCUBIERRE_DEVICE_AUTHORITY: (
            "IsWarpDriveEventRepresented( IsWarpDriveEventRepresented( "
            "IsWarpDriveEventRepresented( sigmaR >= 0.1f && sigmaR <= 100.0f "
            "GetWarpDriveMetric GetWarpDriveChristoffel g_inv[0][0] = -1.0f "
            "g_inv[1][1] = 1.0f - vsf2"
        ),
        ALCUBIERRE_DEVICE_INTEGRATOR: (
            "TryGeodesicAccelerationWarpCart( TryGeodesicAccelerationWarpCart( "
            "TryGeodesicAccelerationWarpCart( TryGeodesicAccelerationWarpCart( "
            "TryGeodesicAccelerationWarpCart( IsWarpDriveEventRepresented(state.x "
            "IntegrateGeodesicRK4WarpCartDelta return IsWarpDriveEventRepresented(candidate"
        ),
        ALCUBIERRE_PARITY_PROBE: (
            "GetWarpDriveMetric GetWarpDriveChristoffel IntegrateGeodesicRK4WarpCart"
        ),
    }
    if alcubierre_authority_errors(valid):
        raise RuntimeError("Alcubierre policy rejected the exact fail-closed seam")

    host_bypass = dict(valid)
    host_bypass[ALCUBIERRE_HOST_AUTHORITY] = host_bypass[ALCUBIERRE_HOST_AUTHORITY].replace(
        "AlcubierreParameterIssue(params.vs, params.R, params.sigma)",
        "std::isfinite(params.vs)",
    )
    if not alcubierre_authority_errors(host_bypass):
        raise RuntimeError("Alcubierre policy accepted a direct host-domain bypass")

    unresolved_device = dict(valid)
    unresolved_device[ALCUBIERRE_DEVICE_AUTHORITY] = unresolved_device[
        ALCUBIERRE_DEVICE_AUTHORITY
    ].replace("sigmaR >= 0.1f && sigmaR <= 100.0f", "sigmaR > 0.0f")
    if not alcubierre_authority_errors(unresolved_device):
        raise RuntimeError("Alcubierre policy accepted an unresolved device wall")

    repaired_inverse = dict(valid)
    repaired_inverse[ALCUBIERRE_DEVICE_AUTHORITY] += " Real det_tx = max(abs(det), 1e-10f);"
    if not alcubierre_authority_errors(repaired_inverse):
        raise RuntimeError("Alcubierre policy accepted a repaired non-exact inverse")

    infallible_stage = dict(valid)
    infallible_stage[ALCUBIERRE_DEVICE_INTEGRATOR] = infallible_stage[
        ALCUBIERRE_DEVICE_INTEGRATOR
    ].replace("TryGeodesicAccelerationWarpCart(", "GeodesicAccelerationWarpCart(", 1)
    if not alcubierre_authority_errors(infallible_stage):
        raise RuntimeError("Alcubierre policy accepted an infallible RK stage")


def kerr_schild_field_equation_errors(documents: dict[Path, str]) -> list[str]:
    """Keep the live metric and its independent field-equation gate inseparable."""
    errors: list[str] = []
    required_code = {
        KERR_SCHILD_HOST_AUTHORITY: (
            "IsRepresentedKerrSchildParameters",
            "ComputeNullVector",
            "result = (Scalar(2.0) * Scalar(p.M) * radius - Scalar(p.Q) * Scalar(p.Q)) / sigma",
            "result = result + Scalar(p.Lambda) * r2 / Scalar(3.0)",
        ),
        KERR_SCHILD_CURVATURE_ORACLE: (
            "RicciFromConnectionFiniteDifference",
            "TensorOps::Christoffel(g, dg)",
            "TensorOps::Inverse(g)",
            "d_gamma",
        ),
        KERR_SCHILD_FIELD_EQUATION_GATE: (
            "KerrSchildFamily metric",
            "KerrSchildParams::Minkowski()",
            "KerrSchildParams::Schwarzschild(",
            "KerrSchildParams::Kerr(",
            "KerrSchildParams::ReissnerNordstrom(",
            "KerrSchildParams::KerrNewman(",
            "KerrSchildParams::DeSitter(",
            "KerrSchildParams{1.0, 0.0, 0.0, 0.01}",
            "require_extremal",
            "require_horizonless",
            "KerrNewmanPotential(",
            "PotentialDerivative(",
            "EightPiMaxwellStress(",
            "RicciFromConnectionFiniteDifference(metric",
            "EXPECT_NEAR(ricci.scalar, 4.0 * sample.parameters.Lambda",
            "EXPECT_NEAR(left, matter",
            "EXPECT_GT(matter_scale, 10.0 * tolerance)",
            "EXPECT_GT(cosmological_scale, 10.0 * tolerance)",
            "KerrNewmanPotentialIsSourceFreeOutsideTheRing",
            "TensorOps::Determinant",
            "volume_density",
        ),
    }
    code_by_path: dict[Path, str] = {}
    for path, markers in required_code.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"Kerr-Schild field-equation participant is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        code_by_path[path] = code
        for marker in markers:
            if marker not in code:
                errors.append(f"{relative(path)} omits Kerr-Schild field-equation marker {marker}")

    oracle = code_by_path.get(KERR_SCHILD_CURVATURE_ORACLE, "")
    gate = code_by_path.get(KERR_SCHILD_FIELD_EQUATION_GATE, "")
    if "metric.InverseMetric(" in oracle:
        errors.append("the curvature oracle reuses the specialised inverse under test")
    if "ComputeH(" in gate:
        errors.append("the independent Maxwell matter oracle reuses production H")

    cmake = documents.get(KERR_SCHILD_TEST_CMAKE, "")
    if "kerr_schild_field_equation_test.cpp" not in cmake:
        errors.append("the Kerr-Schild field-equation gate is not compiled")
    model = documents.get(KERR_SCHILD_OPERATING_MODEL, "")
    for witness in (
        "KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell",
        "KerrSchildFieldEquations.KerrNewmanPotentialIsSourceFreeOutsideTheRing",
    ):
        if witness not in model:
            errors.append(f"the operating model omits Kerr-Schild witness {witness}")
    return errors


def verify_kerr_schild_field_equation_policy() -> None:
    valid = {
        KERR_SCHILD_HOST_AUTHORITY: (
            "IsRepresentedKerrSchildParameters ComputeNullVector "
            "result = (Scalar(2.0) * Scalar(p.M) * radius - Scalar(p.Q) * Scalar(p.Q)) / sigma result = result + Scalar(p.Lambda) * r2 / Scalar(3.0)"
        ),
        KERR_SCHILD_CURVATURE_ORACLE: (
            "RicciFromConnectionFiniteDifference TensorOps::Christoffel(g, dg) "
            "TensorOps::Inverse(g) d_gamma"
        ),
        KERR_SCHILD_FIELD_EQUATION_GATE: (
            "KerrSchildFamily metric KerrSchildParams::Minkowski() "
            "KerrSchildParams::Schwarzschild( KerrSchildParams::Kerr( "
            "KerrSchildParams::ReissnerNordstrom( KerrSchildParams::KerrNewman( "
            "KerrSchildParams::DeSitter( KerrSchildParams{1.0, 0.0, 0.0, 0.01} "
            "require_extremal require_horizonless KerrNewmanPotential( PotentialDerivative( "
            "EightPiMaxwellStress( RicciFromConnectionFiniteDifference(metric "
            "EXPECT_NEAR(ricci.scalar, 4.0 * sample.parameters.Lambda "
            "EXPECT_NEAR(left, matter EXPECT_GT(matter_scale, 10.0 * tolerance) "
            "EXPECT_GT(cosmological_scale, 10.0 * tolerance) "
            "KerrNewmanPotentialIsSourceFreeOutsideTheRing "
            "TensorOps::Determinant volume_density"
        ),
        KERR_SCHILD_TEST_CMAKE: "kerr_schild_field_equation_test.cpp",
        KERR_SCHILD_OPERATING_MODEL: (
            "KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell "
            "KerrSchildFieldEquations.KerrNewmanPotentialIsSourceFreeOutsideTheRing"
        ),
    }
    if kerr_schild_field_equation_errors(valid):
        raise RuntimeError("Kerr-Schild policy rejected the independent field-equation seam")

    specialised_inverse = dict(valid)
    specialised_inverse[KERR_SCHILD_CURVATURE_ORACLE] = specialised_inverse[
        KERR_SCHILD_CURVATURE_ORACLE
    ].replace("TensorOps::Inverse(g)", "metric.InverseMetric(position, inverse)")
    if not kerr_schild_field_equation_errors(specialised_inverse):
        raise RuntimeError("Kerr-Schild policy accepted a coupled specialised inverse")

    reused_h = dict(valid)
    reused_h[KERR_SCHILD_FIELD_EQUATION_GATE] += " metric.ComputeH("
    if not kerr_schild_field_equation_errors(reused_h):
        raise RuntimeError("Kerr-Schild policy accepted production H as a matter oracle")

    oracle_substitution = dict(valid)
    oracle_substitution[KERR_SCHILD_FIELD_EQUATION_GATE] = oracle_substitution[
        KERR_SCHILD_FIELD_EQUATION_GATE
    ].replace("KerrSchildFamily metric", "KerrBoyerLindquist metric")
    if not kerr_schild_field_equation_errors(oracle_substitution):
        raise RuntimeError("Kerr-Schild policy accepted an oracle-only metric substitution")

    missing_cosmological_source = dict(valid)
    missing_cosmological_source[KERR_SCHILD_FIELD_EQUATION_GATE] = (
        missing_cosmological_source[KERR_SCHILD_FIELD_EQUATION_GATE].replace(
            "EXPECT_NEAR(ricci.scalar, 4.0 * sample.parameters.Lambda", ""
        )
    )
    if not kerr_schild_field_equation_errors(missing_cosmological_source):
        raise RuntimeError("Kerr-Schild policy accepted an omitted cosmological equation")

    missing_maxwell_equation = dict(valid)
    missing_maxwell_equation[KERR_SCHILD_FIELD_EQUATION_GATE] = missing_maxwell_equation[
        KERR_SCHILD_FIELD_EQUATION_GATE
    ].replace("KerrNewmanPotentialIsSourceFreeOutsideTheRing", "")
    if not kerr_schild_field_equation_errors(missing_maxwell_equation):
        raise RuntimeError("Kerr-Schild policy accepted an omitted source-free Maxwell gate")

    missing_model_evidence = dict(valid)
    missing_model_evidence[KERR_SCHILD_OPERATING_MODEL] = missing_model_evidence[
        KERR_SCHILD_OPERATING_MODEL
    ].replace("KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell", "")
    if not kerr_schild_field_equation_errors(missing_model_evidence):
        raise RuntimeError("Kerr-Schild policy accepted missing operating-model evidence")
