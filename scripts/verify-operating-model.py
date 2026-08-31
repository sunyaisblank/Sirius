#!/usr/bin/env python3
"""Fail when a required operating-model dimension lacks build-gating evidence."""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

from test_source_governance import SourceTest, collect_source_tests

ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "tests" / "operating_model.json"
LABELS = ROOT / "tests" / "labels" / "CTestLabels.cmake"
REQUIRED_DIMENSIONS = {
    "attestation_admission_and_release_alignment",
    "compile_time_contracts",
    "test_evidence_integrity",
    "test_registration_completeness",
    "configuration_and_input_boundary",
    "install_and_volume_initialisation",
    "operator_script_exit_and_artifact_integrity",
    "runtime_resource_integrity",
    "cpu_render_path",
    "session_lifecycle_and_cancellation",
    "vulkan_capability_and_render_path",
    "sampling_semantics",
    "device_identity_and_allocation_domain",
    "physics_oracles_and_near_extremal_boundary",
    "natural_scaled_trace_domain",
    "polarised_transport_to_film",
    "page_thorne_and_volumetric_transfer",
    "ray_bundles_and_filtered_point_catalogue",
    "camera_worldline_lens_and_film",
    "interactive_viewer_projection",
    "output_encoding",
    "memory_and_dispatch_governance",
    "kernel_portability",
    "metric_catalogue_and_declines",
}
REQUIRED_CAPABILITIES = {
    "revision_bound_release_alignment",
    "polarised_thin_disk_cpu",
    "polarised_volumetric_transfer",
    "polarised_temporal_motion_blur",
    "polarised_vulkan_render",
    "cpu_scalar_motion_blur",
    "vulkan_scalar_motion_blur",
    "covariant_relativistic_jet_transfer",
    "ray_bundles_outside_stationary_black_holes",
    "inverse_compton_corona_transfer",
    "narrowband_line_transfer",
    "newtonian_shakura_sunyaev_temperature_profile",
    "spherical_kottler_cosmological_sector",
    "doppler_suppression_diagnostic",
    "phenomenological_volumetric_disk",
    "disk_emission_outside_schwarzschild_kerr",
    "morris_thorne_one_sheet",
    "morris_thorne_two_sheet",
    "alcubierre_exterior_resolved_profile",
    "vulkan_volumetric_sample_bound",
    "viewer_input_state_logic",
    "viewer_native_window_input",
    "native_p2900_contracts",
    "native_p2996_reflection",
    "physical_radeon_780m_runtime",
    "wsl2_dozen_runtime",
    "windows_native_build_and_runtime",
    "macos_native_build_and_moltenvk",
    "physical_imax_5616x4096",
}
REQUIRED_ACCEPTANCE_CRITERIA = {
    "P1": "build_gated",
    "P2": "build_gated",
    "P3": "attestation_required",
    "P4": "build_gated",
    "P5": "attestation_required",
    "P6": "build_gated",
    "E1": "build_gated",
    "E2": "build_gated",
    "E3": "attestation_required",
    "E4": "build_gated",
}
REQUIRED_ACCEPTANCE_EVIDENCE = {
    "P1": {
        "ctest:OperationalP1.NearExtremalBurnIn",
        "gtest:AnalyticValidationTest.NearExtremalKerrConservesEnergyAngularMomentumAndCarter",
        "gtest:AnalyticValidationTest.SymplecticChartExitAndNullProjectionFailureRemainDistinct",
        "gtest:KernelParity.NearExtremalKerrLiveRenderIntegratorConservesEnergyAngularMomentumAndCarter",
        "gtest:KerrTests.CartesianMetricIsScaleCovariantBelowTheFormerSpinFloor",
        "gtest:KerrTests.SingularKerrDiskDeclinesInsteadOfReceivingAnEpsilonRadius",
        "gtest:CoordinateTransformTests.KerrRadiusPreservesExactZeroAndScaleCovariance",
        "gtest:RK45IntegratorTests.UnrepresentedStageShrinksBeforeMetricEvaluation",
        "gtest:KernelParity.RepresentedSubThresholdKerrMetricIsScaleCovariant",
        "gtest:KernelParity.UnrepresentedKerrStageShrinksBeforeMetricEvaluation",
        "gtest:KernelParity.NullProjectionPreservesConeBranchAndFailsClosed",
        "gtest:KernelParity.PrecisionProbeArtifactsCarryOnlyTheirDeclaredFloat64Capability",
        "gtest:KernelParity.PrecisionRungsConserveNearExtremalKerrWithoutImageComparison",
        "gtest:LivePathConservationTests.NearExtremalKerrEnergyAngularMomentumAndCarter",
        "gtest:ShadowBoundary.KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p",
        "gtest:VulkanRenderSession.KerrNearExtremalBardeenBoundaryAt1080p",
    },
    "P2": {
        "gtest:BeamPropagationTest.SchwarzschildCircularPhotonCongruenceMatchesClosedFormToOnePartPerMillion",
        "gtest:BeamPropagationTest.SchwarzschildRadialCongruenceMatchesClosedFormToOnePartPerMillion",
        "gtest:CpuJacobiOracle.CurvatureScalarMatchesAnalyticKerrOffEquator",
        "gtest:CpuJacobiOracle.RadialPointSourceCongruenceMatchesClosedForm",
        "gtest:CpuJacobiOracle.TidalContractionMatchesAnalyticSchwarzschildAtMatchedEvents",
        "gtest:CpuTraceBoundary.JacobiBundleTerminatesAtTheSameCausalEvent",
        "gtest:KernelBeam.BeamFlagWiresDeviationWithoutMovingDefault",
        "gtest:KernelParity.BeamEllipseRetainsBothAxesAndOutputOrientation",
        "gtest:KernelParity.CelestialTangentBasisIsSharedByBeamAndPointFilter",
        "gtest:KernelParity.DeviceRadialPointSourceCongruenceMatchesClosedForm",
        "gtest:KernelParity.DeviceTidalContractionMatchesAnalyticKerrAtMatchedEvents",
        "gtest:KernelParity.DeviceTidalContractionMatchesAnalyticSchwarzschildAtMatchedEvents",
        "gtest:RayBundleTest.MagnificationComesOnlyFromJacobiMap",
        "gtest:StarfieldGeneratorTests.EllipticalFilterUsesTheBeamSachsBasis",
        "gtest:StarfieldPointTest.BeamFootprintSuppressesStarFlicker",
    },
    "P3": {
        "ctest:OperationalAttestation.FalseExternalEvidenceIsRejected",
        "ctest:OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice",
        "gtest:RenderSessionProbe.SceneEvidenceBindsCanonicalTypedConfiguration",
        "gtest:StarfieldPointTest.CatalogueMeetsSizeFloorAndIsFinite",
        "gtest:StarfieldPointTest.ImaxCatalogueIndexFitsTheTwoGigabyteOperatingEnvelope",
        "gtest:StarfieldPointTest.SpatialIndexMatchesExhaustiveBeamOracle",
        "gtest:StarfieldGeneratorTests.EllipticalFilterUsesTheBeamSachsBasis",
        "gtest:KernelParity.CelestialTangentBasisIsSharedByBeamAndPointFilter",
        "gtest:VulkanRenderSession.CombinedParitySceneRetainsResolvedImageStructure",
        "gtest:VulkanRenderSession.IndexedPointCatalogueReachesLiveKernel",
    },
    "P4": {
        "gtest:AnalyticValidationTest.PageThorneFluxApproachesNewtonianCubicFalloff",
        "gtest:AnalyticValidationTest.PageThorneFluxMatchesIndependentQuadrature",
        "gtest:AnalyticValidationTest.TruncatedPageThorneDiskUsesDeclaredZeroTorqueEdge",
        "gtest:DopplerToggleTest.SuppressionCollapsesDiskAsymmetry",
        "gtest:GeodesicTracerRedshift.NearExtremalInnerDiskEmissionRemainsFinite",
        "gtest:GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile",
        "gtest:KerrOrbitAuthority.EveryCpuConsumerSharesSignedIscoAndCircularEmitterLaw",
        "gtest:KerrOrbitAuthority.UnrepresentedInputsDeclineInsteadOfChangingTheOrbit",
        "gtest:KernelParity.RepresentedSmallKerrSpinDoesNotAliasToSchwarzschildIsco",
        "gtest:AccretionDiskTest.ShakuraSunyaevProfileHasZeroTorqueEdgeAndDeclaredScale",
        "gtest:GeodesicTracerTest.LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile",
        "gtest:KernelParity.ShakuraSunyaevZeroTorqueTemperatureMatchesCoreModel",
        "gtest:GeodesicTracerVolumetric.RedshiftAndDopplerReachTheLiveVolumeSource",
        "gtest:GeodesicTracerVolumetric.TransferAccumulatesAcrossEveryTraversedSegment",
        "gtest:VolumetricDiskClosure.TruncatedGaussianColumnEqualsDeclaredOpticalDepth",
        "gtest:KernelParity.TruncatedGaussianOpacityMatchesFiniteColumnCoreClosure",
        "gtest:KernelParity.DiskEmissionAppliesExactlyOneGFourthFactor",
        "gtest:KernelParity.FullPageThorneDiskTemperatureMatchesIndependentCoreModel",
        "gtest:SpectralValidationTests.KerrDiskTransferMatchesKillingFieldContraction",
        "gtest:SpectralValidationTests.ZamoBranchRemainsTimelikeInsideTheErgosphere",
        "gtest:SpectralValidationTests.PlanckSpectralRadianceMatchesIndependentCodataEquation",
        "gtest:SpectralValidationTests.WienPeakMatchesIndependentNumericalPlanckMaximum",
        "gtest:SpectralValidationTests.StefanBoltzmannExitanceMatchesIndependentHemisphericPlanckIntegral",
        "gtest:SpectralValidationTests.BlackbodyLawsRejectUnrepresentedDomains",
        "gtest:SpectralRadianceTest.BlackbodyBinsDelegateToPlanckAuthorityAndRejectInvalidTemperature",
        "gtest:KernelParity.BlackbodyMatchesIntegratedCoreSpectrum",
        "gtest:SpectralEmissionTest.BolometricDiskAuthorityAppliesExactlyOneGFourthFactor",
        "gtest:VulkanRenderSession.ThinAndVolumetricDopplerSuppressionAffectLiveEmission",
        "gtest:VulkanRenderSession.ProceduralVolumetricTurbulenceReachesLiveKernel",
    },
    "P5": {
        "ctest:OperationalAttestation.FalseExternalEvidenceIsRejected",
        "ctest:OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice",
        "gtest:CameraWorldlineTest.RestScreenRayAndWorldlineComposeOverLensModels",
        "gtest:MemoryGovernor.TwoGigabyteBudgetSeatsAWorkableTile",
        "gtest:ObserverFrameTests.KerrCameraRayIsPastNullAndLaunchFrequencyIsOne",
        "gtest:ObserverFrameTests.MovingCameraLaunchMatchesIndependentLorentzTransform",
        "gtest:RenderSessionProbe.SceneEvidenceBindsCanonicalTypedConfiguration",
        "gtest:ThinLensCameraTest.CentreApertureSharesPinholeProjectionAndFieldOfView",
        "gtest:VulkanRenderSession.CombinedParitySceneRetainsResolvedImageStructure",
        "gtest:VulkanRenderSession.NonSquareMultisamplingCameraAndLensReachLiveKernel",
    },
    "P6": {
        "gtest:EXRRoundTripTests.HDRGradientSurvivesWriteAndRead",
        "gtest:EXRRoundTripTests.PpmRgbaBoundaryAppliesExactlyOneSrgbEncode",
        "gtest:PNGWriterTest.DecodeRoundTripMatchesSRGBEncoding",
        "gtest:RenderSessionProbe.CpuKerrRenderProducesValidPpmThroughTheOwnedWriter",
        "gtest:RenderSessionProbe.FilmAffectsDisplayOutputButNeverLinearExr",
        "gtest:SrgbTransferAuthority.CurveMatchesIndependentIecOracleAcrossBothBranches",
        "gtest:SrgbTransferAuthority.EightBitQuantisationClipsAndDeclinesNonfiniteInputs",
        "gtest:SrgbTransferAuthority.SpectralFacadeDelegatesWithoutChangingClippingSemantics",
        "gtest:XyzSrgbAuthority.HostSpectralFacadesDelegateToTheExactTransform",
        "gtest:XyzSrgbAuthority.MatchesIndependentPrimaryWhiteDerivation",
        "gtest:XyzSrgbAuthority.PreservesExtendedGamutAndPropagatesNonfiniteInputs",
        "gtest:KernelParity.XyzD65ToLinearSrgbMatchesHostAuthority",
    },
    "E1": {
        "gtest:ConfigValidation.MetricMassAndObserverCoordinateRadiusAreIdentityAware",
        "gtest:ConfigValidation.DiskRequestDeclinesForEveryMetricWithoutAnEmissionModel",
        "gtest:CpuTraceBoundary.CentralEventIsInvariantUnderBundleFeatureToggle",
        "gtest:CpuTraceBoundary.EveryAdvertisedCpuMetricConstructsAndTracesOneRay",
        "gtest:MetricRegistryTests.BackendSupportMatchesImplementations",
        "gtest:MetricRegistryTests.EveryCanonicalNameParsesToItsOwnId",
        "gtest:MetricRegistryTests.PositiveLambdaObserverAndHorizonShareOneCausalDomain",
        "gtest:RenderSessionProbe.CpuMorrisThorneRenderCompletes",
        "gtest:RenderSessionProbe.EveryRegisteredCpuMetricCompletesAFrame",
        "gtest:VulkanRenderSession.CpuVulkanAgreeOnMorrisThorneGeometryWithinStatisticalBounds",
    },
    "E2": {
        "gtest:ConfigValidation.PolarisationRequiresRepresentedThinBlackHoleDisk",
        "gtest:GeodesicTracerTest.LiveDiskCrossingCarriesTransportedPhysicalStokesOrientation",
        "gtest:RenderSessionProbe.CpuPolarisationModeConsumesTransportedDiskStokes",
        "gtest:RenderSessionProbe.PolarisedAndTwoSheetRequestsDeclineAtTheTypedBoundary",
        "gtest:WalkerPenrose.BoyerLindquistInitialDataAndAxisExitDeclineWithoutSubstitution",
        "gtest:WalkerPenrose.SchwarzschildEquatorialPerpendicularTransportsRigidly",
        "gtest:WalkerPenroseLivePath.AgreesWithOracleAcrossCharts",
        "gtest:WalkerPenroseLivePath.ConservesConstantAndOrthonormality",
    },
    "E3": {
        "ctest:OperationalAttestation.FalseExternalEvidenceIsRejected",
        "ctest:OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice",
        "gtest:MemoryGovernor.TwoGigabyteBudgetSeatsAWorkableTile",
        "gtest:RenderSessionProbe.CpuKerrRenderProducesValidPngAndExr",
        "gtest:ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames",
        "gtest:VulkanRenderSession.Kerr160x120CompletesAcrossMultipleGovernedTiles",
    },
    "E4": {
        "ctest:OperationalBuildGate.ReleaseInstallRequiresExactPassedEstate",
        "ctest:OperationalBuildPolicy.ReleaseCannotWeakenGates",
        "ctest:OperationalEvidence.RequiredProfileRejectsGoogleTestSkip",
        "ctest:OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims",
        "gtest:KernelParity.KerrSchildMetricMatchesLegacyToOnePartInMillion",
        "gtest:BuildGateAuthority.ReleaseReceiptBindsEveryInstalledProductAtInitialisation",
        "gtest:OracleConnection.AnalyticConnectionAgreesWithMetricDerivatives",
    },
}
CAPABILITY_STATES = {
    "supported",
    "bounded",
    "fail_closed",
    "substituted",
    "attestation_required",
}
SOURCE_AVAILABLE_TEST_FLOOR = 700
REQUIRED_ATTESTATION_PROFILES = {
    "viewer_native_window_input": {"viewer-native-window-input"},
    "physical_radeon_780m_runtime": {"physical-radeon-780m"},
    "wsl2_dozen_runtime": {"wsl2-dozen"},
    "windows_native_build_and_runtime": {
        "windows-native-build",
        "windows-native-vulkan",
    },
    "macos_native_build_and_moltenvk": {
        "macos-native-build",
        "macos-moltenvk",
    },
    "physical_imax_5616x4096": {"physical-imax-5616x4096"},
}
REQUIRED_EXTERNAL_DOMAINS = set().union(*REQUIRED_ATTESTATION_PROFILES.values())
INTERNAL_ACCEPTANCE_PROFILES = {"compile", "cpu", "vulkan"}
REQUIRED_SECTION_POLICY_DIGESTS = {
    "acceptance_criteria": "c3752010fbfb99bcf13ee9cefd61ce8bc412fa4f6e75172b7f94b0621a5e0f0e",
    "required_dimensions": "0227e0061d3521201d9edc94980ee0c51d947bafe9126df2d739c9e1fa69eb47",
    "capability_contracts": "9bab37294cff0e24143efafcf7f83dfc00e2b22c08697086717e8816b475cd87",
}
CONDITIONAL_SOURCE_PATHS = {Path("tests/render/vulkan_render_test.cpp")}


def source_tests():
    return collect_source_tests(ROOT / "tests", ROOT)


def mandatory_tests():
    content = LABELS.read_text(encoding="utf-8")
    found = set()
    block_pattern = re.compile(
        r"set_tests_properties\(\n(?P<names>.*?)"
        r'PROPERTIES LABELS "(?P<labels>[^"]+)"\n\)',
        re.DOTALL,
    )
    for match in block_pattern.finditer(content):
        if "Mandatory" not in match.group("labels").split(";"):
            continue
        for name in re.findall(r"^\s+(\w+\.\w+)\s*$", match.group("names"), re.MULTILINE):
            found.add(name)
    return found


def manual_ctest_labels():
    content = (ROOT / "tests" / "CMakeLists.txt").read_text(encoding="utf-8")
    labels = {}
    pattern = re.compile(r"set_tests_properties\((?P<body>.*?)\)", re.DOTALL)
    for match in pattern.finditer(content):
        body = match.group("body")
        before_properties, separator, properties = body.partition("PROPERTIES")
        if not separator:
            continue
        label_match = re.search(r'LABELS\s+"([^"]+)"', properties)
        if not label_match:
            continue
        for name in re.findall(r"^\s*([A-Za-z0-9_.-]+)\s*$", before_properties, re.MULTILINE):
            labels[name] = set(label_match.group(1).split(";"))
    return labels


def evidence_shape_error(evidence):
    if not isinstance(evidence, list) or not evidence:
        return "no evidence"
    seen = set()
    for item in evidence:
        if not isinstance(item, dict):
            return "malformed evidence record"
        kind = item.get("kind")
        name = item.get("name", "")
        if kind not in {"gtest", "ctest"} or not isinstance(name, str) or not name:
            return "malformed evidence identity"
        identity = (kind, name)
        if identity in seen:
            return f"duplicate evidence {kind}:{name}"
        seen.add(identity)
    return None


def evidence_policy_error(evidence, required_identities):
    if evidence_shape_error(evidence) is not None:
        return None
    actual = {f"{item['kind']}:{item['name']}" for item in evidence}
    missing = sorted(required_identities - actual)
    unexpected = sorted(actual - required_identities)
    if missing or unexpected:
        return f"evidence policy differs: missing={missing}, unexpected={unexpected}"
    return None


def gtest_source_error(record):
    if not record.has_direct_postcondition:
        return "no direct EXPECT_/ASSERT_ postcondition"
    return None


def source_floor_error(test_count):
    if test_count < SOURCE_AVAILABLE_TEST_FLOOR:
        return (
            f"only {test_count} source GoogleTests remain; policy requires at least "
            f"{SOURCE_AVAILABLE_TEST_FLOOR}"
        )
    return None


def mandatory_coverage_error(source_names, mandatory_names):
    uncovered = sorted(set(source_names) - set(mandatory_names))
    if uncovered:
        preview = ", ".join(uncovered[:5])
        suffix = "" if len(uncovered) <= 5 else f" (+{len(uncovered) - 5} more)"
        return f"enabled source GoogleTests are not Mandatory: {preview}{suffix}"
    return None


def ctest_inventory_errors(document, tests, declared_ctests, test_floor):
    if not isinstance(document, dict) or document.get("kind") != "ctestInfo":
        return ["CTest inventory has an unsupported shape"]
    entries = document.get("tests")
    if not isinstance(entries, list):
        return ["CTest inventory has no test list"]

    errors = []
    names = []
    labels_by_name = {}
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("name"), str):
            errors.append("CTest inventory contains a malformed test")
            continue
        name = entry["name"]
        names.append(name)
        labels = set()
        for prop in entry.get("properties", []):
            if isinstance(prop, dict) and prop.get("name") == "LABELS":
                value = prop.get("value")
                if isinstance(value, list):
                    labels.update(item for item in value if isinstance(item, str))
        labels_by_name[name] = labels

    actual = set(names)
    if len(names) != len(actual):
        errors.append("CTest inventory contains duplicate test identities")
    if len(actual) < test_floor:
        errors.append(
            f"CTest inventory has only {len(actual)} tests; policy requires {test_floor}"
        )
    unmandatory = sorted(name for name in actual if "Mandatory" not in labels_by_name[name])
    if unmandatory:
        errors.append(f"live CTest registrations are not Mandatory: {unmandatory[:5]}")

    source_names = set(tests)
    unexpected = sorted(actual - source_names - set(declared_ctests))
    if unexpected:
        errors.append(f"live CTest inventory has undeclared registrations: {unexpected[:5]}")
    non_operational = sorted(
        name
        for name in actual - source_names
        if "Operational" not in labels_by_name[name]
    )
    if non_operational:
        errors.append(
            f"non-source CTest registrations are not Operational: {non_operational[:5]}"
        )

    missing = source_names - actual
    conditional = {
        name
        for name, record in tests.items()
        if record.path.relative_to(ROOT) in CONDITIONAL_SOURCE_PATHS
    }
    unexpected_missing = sorted(missing - conditional)
    if unexpected_missing:
        errors.append(
            f"source GoogleTests are absent from live CTest: {unexpected_missing[:5]}"
        )
    conditional_missing = missing & conditional
    if conditional_missing and conditional_missing != conditional:
        errors.append("live CTest contains only part of a conditional source suite")
    return errors


def attestation_profile_errors(capabilities):
    errors = []
    external_authorities = {}
    for capability in capabilities:
        capability_id = capability.get("id") if isinstance(capability, dict) else None
        expected = REQUIRED_ATTESTATION_PROFILES.get(capability_id)
        state = capability.get("state") if isinstance(capability, dict) else None
        if (state == "attestation_required") != (expected is not None):
            errors.append(
                f"capability {capability_id}: attestation state differs from policy"
            )
            continue
        if expected is None:
            continue
        profiles = capability.get("profiles")
        if (
            not isinstance(profiles, list)
            or any(not isinstance(profile, str) for profile in profiles)
            or set(profiles) != expected
        ):
            errors.append(
                f"capability {capability_id}: external profiles must be "
                f"{sorted(expected)}"
            )
            continue
        for profile in profiles:
            previous = external_authorities.setdefault(profile, capability_id)
            if previous != capability_id:
                errors.append(
                    f"external profile {profile}: multiple capability authorities"
                )
    if set(external_authorities) != REQUIRED_EXTERNAL_DOMAINS:
        errors.append(
            "attestation capability profiles do not equal the exact release domains"
        )
    return errors


def criterion_profile_binding_error(declared_profiles, dependency_profiles):
    allowed_profiles = INTERNAL_ACCEPTANCE_PROFILES | REQUIRED_EXTERNAL_DOMAINS
    unknown_profiles = declared_profiles - allowed_profiles
    if unknown_profiles:
        return f"profiles contain unknown domains {sorted(unknown_profiles)}"
    declared_external = declared_profiles - INTERNAL_ACCEPTANCE_PROFILES
    if declared_external != dependency_profiles:
        return (
            "external profiles and dependency domains differ: "
            f"declared={sorted(declared_external)}, "
            f"dependencies={sorted(dependency_profiles)}"
        )
    return None


def semantic_policy_errors(model):
    errors = []
    for section, required_digest in REQUIRED_SECTION_POLICY_DIGESTS.items():
        payload = json.dumps(
            model.get(section), sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        actual_digest = hashlib.sha256(payload).hexdigest()
        if actual_digest != required_digest:
            errors.append(
                f"{section}: semantic policy digest differs from verifier policy"
            )
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--ctest-inventory", type=Path)
    args = parser.parse_args()
    if args.self_test:
        policy_model = json.loads(MODEL.read_text(encoding="utf-8"))
        drifted_policy_model = json.loads(json.dumps(policy_model))
        drifted_policy_model["required_dimensions"][0]["evidence"] = []
        valid = [{"kind": "gtest", "name": "Suite.Witness"}]
        duplicate = valid + valid
        malformed = [{"kind": "unknown", "name": "Suite.Witness"}]
        direct = SourceTest("Suite", "Witness", Path("witness.cpp"), "EXPECT_EQ(actual, 7);")
        assertion_free = SourceTest("Suite", "Witness", Path("witness.cpp"), "observe();")
        inventory_tests = {
            "Suite.Witness": SourceTest(
                "Suite", "Witness", ROOT / "tests" / "witness_test.cpp", "EXPECT_TRUE(ok);"
            )
        }
        valid_inventory = {
            "kind": "ctestInfo",
            "tests": [
                {
                    "name": "Suite.Witness",
                    "properties": [{"name": "LABELS", "value": ["Mandatory"]}],
                },
                {
                    "name": "Operational.Probe",
                    "properties": [
                        {"name": "LABELS", "value": ["Mandatory", "Operational"]}
                    ],
                },
            ],
        }
        incomplete_inventory = json.loads(json.dumps(valid_inventory))
        incomplete_inventory["tests"] = incomplete_inventory["tests"][1:]
        valid_capabilities = [
            {
                "id": capability_id,
                "state": "attestation_required",
                "profiles": sorted(profiles),
            }
            for capability_id, profiles in REQUIRED_ATTESTATION_PROFILES.items()
        ]
        drifted_capabilities = json.loads(json.dumps(valid_capabilities))
        drifted_capabilities[0]["profiles"] = ["renamed-but-unadmitted-domain"]
        bound_domain = next(iter(REQUIRED_EXTERNAL_DOMAINS))
        if (evidence_shape_error(valid) is not None
                or evidence_shape_error(duplicate) is None
                or evidence_shape_error(malformed) is None
                or evidence_policy_error(valid, {"gtest:Suite.Witness"})
                or evidence_policy_error(valid, {"gtest:Suite.Other"}) is None
                or gtest_source_error(direct) is not None
                or gtest_source_error(assertion_free) is None
                or source_floor_error(SOURCE_AVAILABLE_TEST_FLOOR) is not None
                or source_floor_error(SOURCE_AVAILABLE_TEST_FLOOR - 1) is None
                or mandatory_coverage_error({"Suite.Witness"}, {"Suite.Witness"})
                or mandatory_coverage_error({"Suite.Witness"}, set()) is None
                or ctest_inventory_errors(
                    valid_inventory,
                    inventory_tests,
                    {"Operational.Probe": {"Mandatory", "Operational"}},
                    0,
                )
                or not ctest_inventory_errors(
                    incomplete_inventory,
                    inventory_tests,
                    {"Operational.Probe": {"Mandatory", "Operational"}},
                    0,
                )
                or attestation_profile_errors(valid_capabilities)
                or not attestation_profile_errors(drifted_capabilities)
                or criterion_profile_binding_error({"cpu", bound_domain}, {bound_domain})
                or criterion_profile_binding_error({"cpu"}, set())
                or criterion_profile_binding_error({"cpu", bound_domain}, set()) is None
                or criterion_profile_binding_error({"cpu", "invented"}, set()) is None
                or semantic_policy_errors(policy_model)
                or not semantic_policy_errors(drifted_policy_model)):
            print("operating model: evidence-shape negative control failed", file=sys.stderr)
            return 1
        print(
            "operating model rejected malformed evidence, live-inventory drift, "
            "semantic-policy drift, and external-profile drift"
        )
        return 0

    model = json.loads(MODEL.read_text(encoding="utf-8"))
    if model.get("schema_version") != 3:
        print("operating model: unsupported schema_version", file=sys.stderr)
        return 1
    if model.get("source_available_test_floor") != SOURCE_AVAILABLE_TEST_FLOOR:
        print("operating model: source-available test floor differs from policy",
              file=sys.stderr)
        return 1

    try:
        tests = source_tests()
    except ValueError as error:
        print(f"operating model: {error}", file=sys.stderr)
        return 1
    mandatory = mandatory_tests()
    ctest_labels = manual_ctest_labels()
    failures = []
    failures.extend(semantic_policy_errors(model))
    if args.ctest_inventory is not None:
        try:
            inventory = json.loads(args.ctest_inventory.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            failures.append(f"CTest inventory is unreadable: {error}")
        else:
            failures.extend(
                ctest_inventory_errors(
                    inventory, tests, ctest_labels, SOURCE_AVAILABLE_TEST_FLOOR
                )
            )
    floor_error = source_floor_error(len(tests))
    if floor_error is not None:
        failures.append(floor_error)
    coverage_error = mandatory_coverage_error(tests, mandatory)
    if coverage_error is not None:
        failures.append(coverage_error)

    def check_evidence(owner, evidence):
        shape_error = evidence_shape_error(evidence)
        if shape_error is not None:
            failures.append(f"{owner}: {shape_error}")
            return
        for item in evidence:
            kind = item.get("kind")
            name = item.get("name", "")
            if kind == "gtest":
                if name not in tests:
                    failures.append(f"{owner}: missing test {name}")
                elif name not in mandatory:
                    failures.append(f"{owner}: test is not Mandatory: {name}")
                elif gtest_source_error(tests[name]) is not None:
                    failures.append(
                        f"{owner}: {gtest_source_error(tests[name])}: {name}"
                    )
            elif kind == "ctest":
                if name not in ctest_labels:
                    failures.append(f"{owner}: missing CTest {name}")
                elif not {"Mandatory", "Operational"}.issubset(ctest_labels[name]):
                    failures.append(
                        f"{owner}: CTest is not Mandatory and Operational: {name}"
                    )
            else:
                failures.append(f"{owner}: unknown evidence kind {kind!r}")

    criteria = model.get("acceptance_criteria", [])
    criterion_ids = [criterion.get("id") for criterion in criteria]
    if len(criterion_ids) != len(set(criterion_ids)) or any(
        not item for item in criterion_ids
    ):
        failures.append("acceptance-criterion ids must be present and unique")
    missing_criteria = set(REQUIRED_ACCEPTANCE_CRITERIA) - set(criterion_ids)
    unexpected_criteria = set(criterion_ids) - set(REQUIRED_ACCEPTANCE_CRITERIA)
    if missing_criteria:
        failures.append(f"missing required acceptance criteria: {sorted(missing_criteria)}")
    if unexpected_criteria:
        failures.append(
            f"unexpected acceptance criteria require verifier policy: "
            f"{sorted(unexpected_criteria)}"
        )

    for criterion in criteria:
        criterion_id = criterion.get("id")
        owner = f"acceptance criterion {criterion_id}"
        if criterion.get("required") is not True:
            failures.append(f"{owner}: required must be true")
        expected_verification = REQUIRED_ACCEPTANCE_CRITERIA.get(criterion_id)
        if criterion.get("verification") != expected_verification:
            failures.append(
                f"{owner}: verification must be {expected_verification!r}"
            )
        profiles = criterion.get("profiles")
        if (
            not isinstance(profiles, list)
            or not profiles
            or any(
                not isinstance(profile, str) or not profile.strip()
                for profile in profiles
            )
            or len(profiles) != len(set(profiles))
        ):
            failures.append(f"{owner}: profiles must be unique non-empty strings")
        if not isinstance(criterion.get("requirement"), str) or not criterion[
            "requirement"
        ].strip():
            failures.append(f"{owner}: requirement must be a non-empty string")
        dependencies = criterion.get("capability_dependencies")
        if (
            not isinstance(dependencies, list)
            or any(
                not isinstance(dependency, str) or not dependency.strip()
                for dependency in (dependencies or [])
            )
            or len(dependencies or []) != len(set(dependencies or []))
        ):
            failures.append(
                f"{owner}: capability_dependencies must be unique non-empty strings"
            )
        elif expected_verification == "attestation_required" and not dependencies:
            failures.append(f"{owner}: attestation has no capability dependency")
        elif expected_verification == "build_gated" and dependencies:
            failures.append(f"{owner}: build-gated criterion has external dependencies")
        check_evidence(owner, criterion.get("evidence", []))
        evidence_policy = REQUIRED_ACCEPTANCE_EVIDENCE.get(criterion_id)
        if evidence_policy is not None:
            policy_error = evidence_policy_error(
                criterion.get("evidence", []), evidence_policy
            )
            if policy_error is not None:
                failures.append(f"{owner}: {policy_error}")

    dimensions = model.get("required_dimensions", [])
    ids = [dimension.get("id") for dimension in dimensions]
    if len(ids) != len(set(ids)) or any(not item for item in ids):
        failures.append("dimension ids must be present and unique")
    missing_dimensions = REQUIRED_DIMENSIONS - set(ids)
    unexpected_dimensions = set(ids) - REQUIRED_DIMENSIONS
    if missing_dimensions:
        failures.append(f"missing required dimensions: {sorted(missing_dimensions)}")
    if unexpected_dimensions:
        failures.append(
            f"unexpected dimensions require verifier policy: {sorted(unexpected_dimensions)}"
        )

    for dimension in dimensions:
        evidence = dimension.get("evidence", [])
        check_evidence(dimension.get("id"), evidence)

    capabilities = model.get("capability_contracts", [])
    capability_ids = [capability.get("id") for capability in capabilities]
    if len(capability_ids) != len(set(capability_ids)) or any(
        not item for item in capability_ids
    ):
        failures.append("capability ids must be present and unique")
    missing_capabilities = REQUIRED_CAPABILITIES - set(capability_ids)
    unexpected_capabilities = set(capability_ids) - REQUIRED_CAPABILITIES
    if missing_capabilities:
        failures.append(
            f"missing required capability contracts: {sorted(missing_capabilities)}"
        )
    if unexpected_capabilities:
        failures.append(
            f"unexpected capability contracts require verifier policy: "
            f"{sorted(unexpected_capabilities)}"
        )
    for capability in capabilities:
        owner = f"capability {capability.get('id')}"
        if capability.get("state") not in CAPABILITY_STATES:
            failures.append(f"{owner}: invalid state {capability.get('state')!r}")
        profiles = capability.get("profiles")
        if (
            not isinstance(profiles, list)
            or not profiles
            or any(
                not isinstance(profile, str) or not profile.strip()
                for profile in profiles
            )
            or len(profiles) != len(set(profiles))
        ):
            failures.append(f"{owner}: profiles must be unique non-empty strings")
        for field in ("request", "behavior"):
            if not isinstance(capability.get(field), str) or not capability[field].strip():
                failures.append(f"{owner}: {field} must be a non-empty string")
        check_evidence(owner, capability.get("evidence", []))

    failures.extend(attestation_profile_errors(capabilities))

    capability_by_id = {item.get("id"): item for item in capabilities}
    for criterion in criteria:
        owner = f"acceptance criterion {criterion.get('id')}"
        declared_profiles = {
            profile
            for profile in criterion.get("profiles", [])
            if isinstance(profile, str)
        }
        dependency_profiles = set()
        for dependency in criterion.get("capability_dependencies", []):
            capability = capability_by_id.get(dependency)
            if capability is None:
                failures.append(f"{owner}: unknown capability dependency {dependency!r}")
            elif capability.get("state") != "attestation_required":
                failures.append(
                    f"{owner}: dependency {dependency!r} is not attestation_required"
                )
            else:
                dependency_profiles.update(
                    profile
                    for profile in capability.get("profiles", [])
                    if isinstance(profile, str)
                )
        binding_error = criterion_profile_binding_error(
            declared_profiles, dependency_profiles
        )
        if binding_error is not None:
            failures.append(f"{owner}: {binding_error}")

    if failures:
        for failure in failures:
            print(f"operating model: {failure}", file=sys.stderr)
        return 1
    print(
        f"operating model requires {len(criteria)} P/E acceptance criteria across "
        f"{len(tests)} governed source GoogleTests, covers {len(dimensions)} operating "
        f"dimensions, and defines {len(capabilities)} capability contracts with "
        "Mandatory evidence"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
