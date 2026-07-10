#!/usr/bin/env python3
"""Regenerate src/Sirius.Test/CTestLabels.cmake from the test sources.

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
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEST_DIR = ROOT / "src" / "Sirius.Test"
OUTPUT = TEST_DIR / "CTestLabels.cmake"

MANDATORY = "Mandatory;Correctness"
STABILITY = "Mandatory;Stability"
CORRECTNESS = "Correctness"
PERFORMANCE = "Performance"

# Suite -> labels. A suite absent from this map gets CORRECTNESS and a
# warning, so new suites are never silently ungated; DISABLED_ tests are
# skipped by gtest itself and need no exclusion here.
SUITE_LABELS = {
    # Mathematical foundations (build gate)
    "DualNumberTests": MANDATORY,
    "TensorTests": MANDATORY,
    "TensorInverseTests": MANDATORY,
    "MetricRegistryTests": MANDATORY,
    "SchwarzschildTests": MANDATORY,
    "ChristoffelTests": MANDATORY,
    "KerrTests": MANDATORY,
    "KerrSchildTests": MANDATORY,
    "ReissnerNordstromTests": MANDATORY,
    "EllisDrainholeTests": MANDATORY,
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
    "GeodesicPathTests": CORRECTNESS,
    "BeamPropagationTest": CORRECTNESS,
    "EinsteinRingTest": CORRECTNESS,
    "AccretionDiskTest": CORRECTNESS,
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
    "GeodesicTracerTest": CORRECTNESS,
    "PNGWriterTest": CORRECTNESS,
    "EXRRoundTripTests": CORRECTNESS,
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
}


def collect_tests():
    tests = {}
    pattern = re.compile(r"TEST(?:_F)?\(\s*(\w+)\s*,\s*(\w+)\s*\)")
    for path in sorted(TEST_DIR.rglob("*.cpp")):
        for suite, name in pattern.findall(path.read_text(encoding="utf-8")):
            tests.setdefault(suite, []).append(name)
    return tests


def render(tests):
    lines = [
        "# =============================================================================",
        "# CTestLabels.cmake - Test labels for mandatory gating",
        "# =============================================================================",
        "# GENERATED FILE - do not edit by hand.",
        "# Regenerate with: python3 scripts/generate-ctest-labels.py",
        "# CI verifies freshness with --check; suite renames that would silently",
        "# detach labels now fail the build instead.",
        "#",
        "# Labels:",
        "#   Mandatory   - Build fails if these tests fail",
        "#   Correctness - Mathematical/behavioural invariant tests",
        "#   Stability   - Numerical stability tests",
        "#   Performance - FPS/resource metrics (warning only)",
        "# =============================================================================",
        "",
    ]
    unknown = []
    for suite in sorted(tests):
        labels = SUITE_LABELS.get(suite)
        if labels is None:
            unknown.append(suite)
            labels = CORRECTNESS
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
    tests = collect_tests()
    content, unknown = render(tests)
    for suite in unknown:
        print(f"warning: suite '{suite}' not in SUITE_LABELS; defaulted to Correctness",
              file=sys.stderr)

    if "--check" in sys.argv:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != content:
            print("CTestLabels.cmake is stale; run scripts/generate-ctest-labels.py",
                  file=sys.stderr)
            return 1
        print("CTestLabels.cmake is up to date")
        return 0

    OUTPUT.write_text(content, encoding="utf-8")
    print(f"wrote {OUTPUT} ({len(tests)} suites)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
