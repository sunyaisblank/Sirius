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

import sys
import tempfile
from pathlib import Path

from test_source_governance import collect_source_tests

ROOT = Path(__file__).resolve().parent.parent
NEW_TEST_DIR = ROOT / "tests"
NEW_OUTPUT = NEW_TEST_DIR / "labels" / "CTestLabels.cmake"

MANDATORY = "Mandatory;Correctness"
STABILITY = "Mandatory;Stability"
CORRECTNESS = "Mandatory;Correctness"
OPERATIONAL = "Mandatory;Operational"

# Suite -> labels. An absent suite is a hard policy error: assigning a default
# would let a new operating surface silently escape the build gate.
SUITE_LABELS = {
    # Mathematical foundations (build gate)
    "DualNumberTests": MANDATORY,
    "TensorTests": MANDATORY,
    "TensorInverseTests": MANDATORY,
    "MetricRegistryTests": MANDATORY,
    "SchwarzschildTests": MANDATORY,
    "KottlerHorizonTests": MANDATORY,
    "ChristoffelTests": MANDATORY,
    "KerrTests": MANDATORY,
    "KerrSchildTests": MANDATORY,
    "ReissnerNordstromTests": MANDATORY,
    "EllisDrainholeTests": MANDATORY,
    "MorrisThorneCartesianTests": MANDATORY,
    "AlcubierreMetricTests": MANDATORY,
    "MetricLoaderChainTests": MANDATORY,
    "MetricDerivativeTests": MANDATORY,
    "MetricValidationTests": MANDATORY,
    "CoordinateTransformTests": MANDATORY,
    "GeodesicTests": MANDATORY,
    "RK45IntegratorTests": MANDATORY,
    "ConservationLawTests": MANDATORY,
    "LivePathConservationTests": MANDATORY,
    "ObserverFrameTests": MANDATORY,
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
    # Numerical stability (build gate)
    "NumericalStabilityTests": STABILITY,
    "NumericalStabilityTest": STABILITY,
    "NaNInfDetectionTests": STABILITY,
    "GPUConservationTests": STABILITY,
    "CPUGPUParityTests": STABILITY,
    # Feature-level correctness. Every enabled source test is build-gating;
    # category labels remain useful for focused local and CI runs.
    "GeodesicPathTests": MANDATORY,
    "BeamPropagationTest": MANDATORY,
    "AccretionDiskTest": MANDATORY,
    "DiskCoordinateTest": MANDATORY,
    "TurbulenceTest": CORRECTNESS,
    "CoronaConfigTests": CORRECTNESS,
    "CoronaEmissivityTests": CORRECTNESS,
    "CoronaGeometryTests": CORRECTNESS,
    "JetGeometryTests": CORRECTNESS,
    "JetEmissionTests": CORRECTNESS,
    "JetDopplerTests": CORRECTNESS,
    "SpectralRadianceTest": CORRECTNESS,
    "SpectralValidationTests": CORRECTNESS,
    "SpectralUtilsTests": CORRECTNESS,
    "PolarisedEmissionTests": CORRECTNESS,
    "FIDOFrameTests": CORRECTNESS,
    "GFactorTests": CORRECTNESS,
    "RenderPipelineTests": CORRECTNESS,
    "GeodesicTracerTest": OPERATIONAL,
    "GeodesicTracerVolumetric": OPERATIONAL,
    "GeodesicTracerRedshift": OPERATIONAL,
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
    # New-tree suites (no legacy counterpart)
    "Contracts": MANDATORY,
    "ContractsDeathTest": MANDATORY,
    "Error": MANDATORY,
    "Sha256Tests": MANDATORY,
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
    "AlignmentAuthority": OPERATIONAL,
    "BuildGateAuthority": OPERATIONAL,
    "MemoryGovernor": OPERATIONAL,
    "DispatchGovernor": OPERATIONAL,
    "VulkanRenderSession": OPERATIONAL,
    "RayBundleTest": MANDATORY,
    "CameraWorldlineTest": MANDATORY,
    "DopplerToggleTest": OPERATIONAL,
    "StarfieldPointTest": OPERATIONAL,
    "KernelBeam": MANDATORY,
    "KernelPortability": OPERATIONAL,
    "PixelSampling": OPERATIONAL,
    "ShadowBoundary": MANDATORY,
}


def collect_new_tree_tests(test_dir=NEW_TEST_DIR):
    tests = {}
    for record in collect_source_tests(test_dir, ROOT).values():
        tests.setdefault(record.suite, []).append(record.case)
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
                else:
                    print(f"error: unexpected rejection: {error}", file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted an ungoverned parameterised test",
                      file=sys.stderr)
                return 1

            fixture.write_text(
                "TEST(EscapeSuite, NoOp) { EXPECT_TRUE(true); }\\n",
                encoding="utf-8",
            )
            try:
                collect_new_tree_tests(Path(directory))
            except ValueError as error:
                if "no-op assertion" in str(error):
                    print("source governance rejected an obvious no-op assertion")
                else:
                    print(f"error: unexpected rejection: {error}", file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted an obvious no-op assertion",
                      file=sys.stderr)
                return 1

            fixture.write_text(
                "TEST(EscapeSuite, AssertionFree) { observe(); }\\n",
                encoding="utf-8",
            )
            try:
                collect_new_tree_tests(Path(directory))
            except ValueError as error:
                if "no direct" in str(error):
                    print("source governance rejected an assertion-free test")
                else:
                    print(f"error: unexpected rejection: {error}", file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted an assertion-free test", file=sys.stderr)
                return 1

            fixture.write_text(
                "TEST(EscapeSuite, DISABLED_Hidden) { EXPECT_TRUE(observe()); }\n",
                encoding="utf-8",
            )
            try:
                collect_new_tree_tests(Path(directory))
            except ValueError as error:
                if "disabled test" in str(error):
                    print("source governance rejected a disabled test")
                else:
                    print(f"error: unexpected rejection: {error}", file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted a disabled test", file=sys.stderr)
                return 1

            fixture.write_text(
                "#if 0\nTEST(EscapeSuite, Inactive) { EXPECT_TRUE(observe()); }\n#endif\n",
                encoding="utf-8",
            )
            try:
                collect_new_tree_tests(Path(directory))
            except ValueError as error:
                if "inactive #if 0" in str(error):
                    print("source governance rejected an inactive test block")
                    return 0
                print(f"error: unexpected rejection: {error}", file=sys.stderr)
                return 1
            print("error: source governance accepted an inactive test block", file=sys.stderr)
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
