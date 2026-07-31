#!/usr/bin/env python3
"""Regenerate tests/labels/CTestLabels.cmake from the test sources.

The labels file must reference tests by their exact registered names
(Suite.TestName), which historically drifted whenever suites were renamed:
labels silently stopped attaching and the Mandatory gate lost coverage
without any signal. This generator makes the source tree the single
authority: it extracts every TEST/TEST_F, applies the suite-to-label policy
below, and rewrites the file. CI runs it with --check to fail on drift.

Usage:
  python3 scripts/generate-ctest-labels.py           # rewrite the file
  python3 scripts/generate-ctest-labels.py --check   # exit 1 if stale
"""

import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NEW_TEST_DIR = ROOT / "tests"
NEW_OUTPUT = NEW_TEST_DIR / "labels" / "CTestLabels.cmake"

MANDATORY = "Mandatory;Correctness"
STABILITY = "Mandatory;Stability"
CORRECTNESS = "Correctness"
PERFORMANCE = "Performance"
OPERATIONAL = "Mandatory;Operational"

# Suite -> labels. An absent suite is a hard policy error: assigning a default
# would let a new operating surface silently escape the build gate.
SUITE_LABELS = {
    # Mathematical foundations (build gate)
    "DualNumberTests": MANDATORY,
    "TensorTests": MANDATORY,
    "TensorInverseTests": MANDATORY,
    "ChristoffelTransformTests": MANDATORY,
    "MetricRegistryTests": MANDATORY,
    "SchwarzschildTests": MANDATORY,
    "ChristoffelTests": MANDATORY,
    "KerrTests": MANDATORY,
    "KerrSchildTests": MANDATORY,
    "ReissnerNordstromTests": MANDATORY,
    "EllisDrainholeTests": MANDATORY,
    "MorrisThorneCartesianTests": MANDATORY,
    "AlcubierreMetricTests": MANDATORY,
    "NumericalMetricTests": MANDATORY,
    "MetricLoaderChainTests": MANDATORY,
    "MetricDerivativeTests": MANDATORY,
    "MetricValidationTests": MANDATORY,
    "CoordinateTransformTests": MANDATORY,
    "GeodesicTests": MANDATORY,
    "RK45IntegratorTests": MANDATORY,
    "ConservationLawTests": MANDATORY,
    "LivePathConservationTests": MANDATORY,
    "CaptureSurfaceTests": MANDATORY,
    "DeterminismTests": MANDATORY,
    "MandatoryKillingTests": MANDATORY,
    "SymplecticIntegratorTest": MANDATORY,
    "KerrMetricDTest": MANDATORY,
    "Vec4dTests": MANDATORY,
    "Mat4dTests": MANDATORY,
    "Vec4dTest": MANDATORY,
    "Mat4dTest": MANDATORY,
    "GeodesicStateDTests": MANDATORY,
    "HamiltonianStateDTests": MANDATORY,
    "ConversionTests": MANDATORY,
    "AnalyticValidationTest": MANDATORY,
    "AnalyticValidationTests": MANDATORY,
    "DNGRParityTest": MANDATORY,
    "StokesVectorTests": MANDATORY,
    "MuellerMatrixTests": MANDATORY,
    "ParallelTransportTests": MANDATORY,
    "PinholeCameraTest": MANDATORY,
    "ThinLensCameraTest": MANDATORY,
    "FisheyeCameraTest": MANDATORY,
    "ICameraTest": MANDATORY,
    "CameraFactoryTest": MANDATORY,
    "SpectralEmissionTest": MANDATORY,
    "BlackbodyParityTest": MANDATORY,
    # Numerical stability (build gate)
    "NumericalStabilityTests": STABILITY,
    "NumericalStabilityTest": STABILITY,
    "NaNInfDetectionTests": STABILITY,
    "GPUConservationTests": STABILITY,
    "CPUGPUParityTests": STABILITY,
    # Correctness, not gating (feature-level behaviour)
    "GeodesicPathTests": MANDATORY,
    "BeamPropagationTest": MANDATORY,
    "EinsteinRingTest": MANDATORY,
    "AccretionDiskTest": MANDATORY,
    "TurbulenceTest": CORRECTNESS,
    "CoronaConfigTests": CORRECTNESS,
    "CoronaEmissivityTests": CORRECTNESS,
    "CoronaGeometryTests": CORRECTNESS,
    "JetGeometryTests": CORRECTNESS,
    "JetEmissionTests": CORRECTNESS,
    "JetDopplerTests": CORRECTNESS,
    "JetRayMarchTests": CORRECTNESS,
    "SpectralRadianceTest": CORRECTNESS,
    "SpectralValidationTests": CORRECTNESS,
    "SpectralUtilsTests": CORRECTNESS,
    "PolarisedEmissionTests": CORRECTNESS,
    "GeodesicDeviationTests": CORRECTNESS,
    "FIDOFrameTests": CORRECTNESS,
    "GFactorTests": CORRECTNESS,
    "SMBHParamsTest": CORRECTNESS,
    "RenderPipelineTests": CORRECTNESS,
    "GeodesicTracerTest": OPERATIONAL,
    "GeodesicTracerVolumetric": OPERATIONAL,
    "MorrisThorneTracerTest": CORRECTNESS,
    "PNGWriterTest": OPERATIONAL,
    "EXRRoundTripTests": OPERATIONAL,
    "FilmSimulationTest": CORRECTNESS,
    "StarEntryTests": CORRECTNESS,
    "StarfieldConfigTests": CORRECTNESS,
    "StarfieldGeneratorTests": CORRECTNESS,
    "TonemapTests": CORRECTNESS,
    "BloomFilterTests": CORRECTNESS,
    "ColourGradingTests": CORRECTNESS,
    # Performance (warning only)
    "FPSBenchmarks": PERFORMANCE,
    "FPSThresholdTests": PERFORMANCE,
    "MemoryUsageTests": PERFORMANCE,
    "GeodesicBenchmarks": PERFORMANCE,
    "ChristoffelBenchmark": PERFORMANCE,
    # New-tree suites (no legacy counterpart)
    "Contracts": MANDATORY,
    "ContractsDeathTest": MANDATORY,
    "Error": MANDATORY,
    "KernelParity": MANDATORY,
    "KernelTrace": OPERATIONAL,
    "VulkanBackend": OPERATIONAL,
    "RenderSessionProbe": OPERATIONAL,
    "TileScheduler": OPERATIONAL,
    "DisplayBuffer": OPERATIONAL,
    "OracleConnection": MANDATORY,
    "WalkerPenrose": MANDATORY,
    "WalkerPenroseLivePath": MANDATORY,
    "CpuGeodesicReferenceTests": MANDATORY,
    "ConfigValidation": MANDATORY,
    "ConfigEnvironment": OPERATIONAL,
    "ConfigLoading": OPERATIONAL,
    "ConfigSchema": MANDATORY,
    "RenderCommandParse": OPERATIONAL,
    "ViewCommandOperational": OPERATIONAL,
    "CommandRouter": OPERATIONAL,
    "PlatformPaths": OPERATIONAL,
    "MemoryGovernor": OPERATIONAL,
    "DispatchGovernor": OPERATIONAL,
    "VulkanRenderSession": OPERATIONAL,
    "RayBundleTest": MANDATORY,
    "CameraAberrationTest": MANDATORY,
    "DopplerToggleTest": OPERATIONAL,
    "StarfieldPointTest": OPERATIONAL,
    "KernelBeam": CORRECTNESS,
    "KernelPortability": OPERATIONAL,
    "PixelSampling": OPERATIONAL,
    "ShadowBoundary": MANDATORY,
}


def collect_new_tree_tests(test_dir=NEW_TEST_DIR):
    tests = {}
    pattern = re.compile(r"TEST(?:_F)?\(\s*(\w+)\s*,\s*(\w+)\s*\)")
    unsupported = re.compile(
        r"^\s*(TEST_P|TYPED_TEST|TYPED_TEST_P)\s*\(",
        re.MULTILINE,
    )
    for path in sorted(test_dir.rglob("*_test.cpp")):
        content = path.read_text(encoding="utf-8")
        macro = unsupported.search(content)
        if macro:
            try:
                display_path = path.relative_to(ROOT)
            except ValueError:
                display_path = path
            raise ValueError(
                f"{display_path} uses {macro.group(1)}, whose discovered "
                "names are not covered by the Mandatory-label generator"
            )
        for suite, name in pattern.findall(content):
            tests.setdefault(suite, []).append(name)
    return tests


def render_new_tree(tests):
    # The new tree registers tests at build time (gtest_discover_tests), so
    # this file is loaded by ctest AFTER each module's discovery include, via
    # TEST_INCLUDE_FILES per module directory (the legacy mechanism). Entries
    # for tests a directory does not own are inert: ctest-time properties for
    # never-registered names attach to nothing. No if(TEST) guards; that
    # predicate does not evaluate during ctest processing.
    lines = [
        "# GENERATED FILE - do not edit by hand.",
        "# Regenerate with: python3 scripts/generate-ctest-labels.py",
        "# Labels policy lives in that script; CI fails on a stale copy.",
        "",
    ]
    unknown = []
    for suite in sorted(tests):
        labels = SUITE_LABELS.get(suite)
        if labels is None:
            unknown.append(suite)
            continue
        names = [
            f"    {suite}.{name}"
            for name in tests[suite]
            if not name.startswith("DISABLED_")
        ]
        if not names:
            continue
        lines.append("set_tests_properties(")
        lines.extend(names)
        lines.append(f'    PROPERTIES LABELS "{labels}"')
        lines.append(")")
        lines.append("")
    return "\n".join(lines) + "\n", unknown


def main():
    if "--self-test" in sys.argv:
        with tempfile.TemporaryDirectory(prefix="sirius-label-governance-") as directory:
            fixture = Path(directory) / "escape_test.cpp"
            fixture.write_text(
                "TEST_P(EscapeSuite, Ungoverned) {}\\n",
                encoding="utf-8",
            )
            try:
                collect_new_tree_tests(Path(directory))
            except ValueError as error:
                if "TEST_P" in str(error):
                    print("source governance rejected an ungoverned parameterised test")
                    return 0
            print("error: source governance accepted an ungoverned parameterised test",
                  file=sys.stderr)
            return 1

    try:
        new_tests = collect_new_tree_tests()
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    new_content, new_unknown = render_new_tree(new_tests)
    for suite in new_unknown:
        print(f"error: suite '{suite}' has no explicit label policy", file=sys.stderr)
    if new_unknown:
        return 1

    if "--check" in sys.argv:
        new_current = NEW_OUTPUT.read_text(encoding="utf-8") if NEW_OUTPUT.exists() else ""
        if new_current != new_content:
            print("CTest label file is stale; run scripts/generate-ctest-labels.py",
                  file=sys.stderr)
            return 1
        print("CTest label file is up to date")
        return 0

    NEW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    NEW_OUTPUT.write_text(new_content, encoding="utf-8")
    print(f"wrote {NEW_OUTPUT} ({len(new_tests)} suites)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
