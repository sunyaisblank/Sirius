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

from test_source_governance import APP_RENDERING_TESTS, collect_source_tests

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
    "KerrOrbitAuthority": MANDATORY,
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
    "KerrZamoTransfer": MANDATORY,
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
    "VolumetricDiskClosure": MANDATORY,
    "TurbulenceTest": CORRECTNESS,
    "CoronaConfigTests": CORRECTNESS,
    "CoronaEmissivityTests": CORRECTNESS,
    "CoronaGeometryTests": CORRECTNESS,
    "JetGeometryTests": CORRECTNESS,
    "JetEmissionTests": CORRECTNESS,
    "JetDopplerTests": CORRECTNESS,
    "SpectralRadianceTest": CORRECTNESS,
    "Cie1931ObserverAuthority": CORRECTNESS,
    "SrgbTransferAuthority": CORRECTNESS,
    "XyzSrgbAuthority": CORRECTNESS,
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
    "ViewerDisplayContract": MANDATORY,
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
    "CpuJacobiOracle": MANDATORY,
    "CpuTraceBoundary": MANDATORY,
    "KernelPortability": OPERATIONAL,
    "PixelSampling": OPERATIONAL,
    "ShadowBoundary": MANDATORY,
}


def collect_new_tree_tests(test_dir=NEW_TEST_DIR):
    tests = {}
    rendering = set()
    for record in collect_source_tests(test_dir, ROOT).values():
        tests.setdefault(record.suite, []).append(record.case)
        try:
            relative = record.path.relative_to(ROOT)
        except ValueError:
            relative = record.path
        if len(relative.parts) >= 2 and relative.parts[0:2] == ("tests", "render"):
            rendering.add(record.name)
    return tests, rendering


def render_new_tree(tests, rendering_tests):
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
            if not name.startswith("DISABLED_") and f"{suite}.{name}" not in rendering_tests
        ]
        if not names:
            continue
        lines.append("set_tests_properties(")
        lines.extend(names)
        lines.append(f'    PROPERTIES LABELS "{labels}"')
        lines.append(")")
        lines.append("")

    rendering_by_labels = {}
    for name in sorted(rendering_tests):
        suite, _, _ = name.partition(".")
        labels = SUITE_LABELS.get(suite)
        if labels is None:
            continue
        rendering_by_labels.setdefault(f"{labels};Rendering", []).append(name)
    for labels, names in sorted(rendering_by_labels.items()):
        lines.append("set_tests_properties(")
        lines.extend(f"    {name}" for name in names)
        lines.append(f'    PROPERTIES LABELS "{labels}"')
        lines.append(")")
        lines.append("")
    return "\n".join(lines) + "\n", unknown, sorted(APP_RENDERING_TESTS - rendering_tests)


def main():
    if "--self-test" in sys.argv:
        with tempfile.TemporaryDirectory(prefix="sirius-label-governance-") as directory:
            label_fixture, unknown_fixture, _ = render_new_tree(
                {"RenderPipelineTests": ["PureBoundary", "RenderedBoundary"]},
                {"RenderPipelineTests.RenderedBoundary"},
            )
            if (unknown_fixture or
                    "RenderPipelineTests.PureBoundary\n"
                    "    PROPERTIES LABELS \"Mandatory;Correctness\"" not in label_fixture or
                    "RenderPipelineTests.RenderedBoundary\n"
                    "    PROPERTIES LABELS \"Mandatory;Correctness;Rendering\""
                    not in label_fixture):
                print("error: rendering-directory label inheritance is not structural",
                      file=sys.stderr)
                return 1

            boundary_root = Path(directory)
            app_fixture = boundary_root / "tests" / "app" / "render_escape_test.cpp"
            app_fixture.parent.mkdir(parents=True)
            app_fixture.write_text(
                "TEST(ViewCommandOperational, HeadlessRefinementProducesASynchronisedFrame) {\n"
                "  InteractiveViewer viewer; ASSERT_TRUE(viewer.Start());\n}\n",
                encoding="utf-8",
            )
            try:
                collect_source_tests(boundary_root / "tests", boundary_root)
            except ValueError as error:
                if "outside tests/render" not in str(error):
                    print(f"error: unexpected render-boundary rejection: {error}",
                          file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted a rendering test in tests/app",
                      file=sys.stderr)
                return 1

            render_fixture = boundary_root / "tests" / "render" / app_fixture.name
            render_fixture.parent.mkdir(parents=True)
            app_fixture.replace(render_fixture)
            records = collect_source_tests(boundary_root / "tests", boundary_root)
            if set(records) != {
                    "ViewCommandOperational.HeadlessRefinementProducesASynchronisedFrame"}:
                print("error: source governance rejected the rendering test in tests/render",
                      file=sys.stderr)
                return 1
            render_fixture.unlink()

            app_fixture.write_text(
                "TEST(RenderCommandParse, NewlyAddedExecution) {\n"
                "  RenderCommand command; GlobalOptions globals; SiriusConfig config;\n"
                "  EXPECT_EQ(command.Execute({}, globals, config), 0);\n}\n",
                encoding="utf-8",
            )
            try:
                collect_source_tests(boundary_root / "tests", boundary_root)
            except ValueError as error:
                if "unclassified RenderCommand execution" not in str(error):
                    print(f"error: unexpected command-boundary rejection: {error}",
                          file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted unclassified app-side execution",
                      file=sys.stderr)
                return 1
            app_fixture.unlink()

            app_fixture.write_text(
                "int HiddenRender() {\n"
                "  RenderCommand command; GlobalOptions globals; SiriusConfig config;\n"
                "  return command.Execute({}, globals, config);\n}\n"
                "TEST(RenderCommandParse, HelperEscape) { EXPECT_EQ(HiddenRender(), 0); }\n",
                encoding="utf-8",
            )
            try:
                collect_source_tests(boundary_root / "tests", boundary_root)
            except ValueError as error:
                if "hides app-side Start/Execute" not in str(error):
                    print(f"error: unexpected helper-boundary rejection: {error}",
                          file=sys.stderr)
                    return 1
            else:
                print("error: source governance accepted a hidden app-side execution helper",
                      file=sys.stderr)
                return 1
            app_fixture.unlink()

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
        new_tests, rendering_tests = collect_new_tree_tests()
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    new_content, new_unknown, missing_rendering = render_new_tree(new_tests, rendering_tests)
    for suite in new_unknown:
        print(f"error: suite '{suite}' has no explicit label policy", file=sys.stderr)
    if new_unknown:
        return 1
    for name in missing_rendering:
        print(f"error: governed rendering test is missing: {name}", file=sys.stderr)
    if missing_rendering:
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
