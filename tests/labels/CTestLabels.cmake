# GENERATED FILE - do not edit by hand.
# Regenerate with: python3 scripts/generate-ctest-labels.py
# Labels policy lives in that script; CI fails on a stale copy.

set_tests_properties(
    AccretionDiskTest.ISCO_Schwarzschild
    AccretionDiskTest.NonFiniteConfigurationFailsClosedWithoutRewritingTheRequest
    AccretionDiskTest.OuterEdgeInsideDerivedIscoFailsClosed
    AccretionDiskTest.ExplicitInnerEdgeCannotEnterTheUnstableOrbitDomain
    AccretionDiskTest.ShakuraSunyaevProfileHasZeroTorqueEdgeAndDeclaredScale
    AccretionDiskTest.ISCO_ExtremalKerr_Prograde
    AccretionDiskTest.ISCO_ExtremalKerr_Retrograde
    AccretionDiskTest.ISCO_ModerateSpin
    AccretionDiskTest.TemperatureProfileShape
    AccretionDiskTest.TemperatureScalingWithMass
    AccretionDiskTest.DiskBoundaries
    AccretionDiskTest.SpectralEmission
    AccretionDiskTest.FluxProfile
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    AlcubierreMetricTests.MetricSymmetry
    AlcubierreMetricTests.LorentzianSignature
    AlcubierreMetricTests.ReducesToMinkowskiOutside
    AlcubierreMetricTests.ShiftVectorAtCentre
    AlcubierreMetricTests.GttAtCentre
    AlcubierreMetricTests.SpatialComponentsFlat
    AlcubierreMetricTests.ShapeFunctionAtCentre
    AlcubierreMetricTests.ShapeFunctionFarField
    AlcubierreMetricTests.ShapeFunctionMonotoneDecrease
    AlcubierreMetricTests.ShapeFunctionDerivativeFinite
    AlcubierreMetricTests.SubluminalConstruction
    AlcubierreMetricTests.SuperluminalConstruction
    AlcubierreMetricTests.SetParameterMatchesTheOperatorVelocityDomain
    AlcubierreMetricTests.NoNaNAtBubbleWall
    AlcubierreMetricTests.NoNaNAtOrigin
    AlcubierreMetricTests.BubblePositionUpdate
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    AlignmentAuthority.CompiledReceiptMatchesTheStagedRuntimeAuthority
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    AnalyticValidationTest.PhotonSphereSchwarzschildExact
    AnalyticValidationTest.PhotonSphereKerrPrograde
    AnalyticValidationTest.PhotonSphereKerrRetrograde
    AnalyticValidationTest.ISCOSchwarzschildExact
    AnalyticValidationTest.ISCOKerrPrograde
    AnalyticValidationTest.ISCONearExtremal
    AnalyticValidationTest.PageThorneFluxMatchesIndependentQuadrature
    AnalyticValidationTest.PageThorneTemperatureHasZeroTorqueInnerEdge
    AnalyticValidationTest.TruncatedPageThorneDiskUsesDeclaredZeroTorqueEdge
    AnalyticValidationTest.PageThorneFluxApproachesNewtonianCubicFalloff
    AnalyticValidationTest.MalformedBeamEscapeDomainFailsClosed
    AnalyticValidationTest.MalformedIntegratorControlsFailClosedBeforeIntegration
    AnalyticValidationTest.SymplecticChartExitAndNullProjectionFailureRemainDistinct
    AnalyticValidationTest.NearExtremalKerrConservesEnergyAngularMomentumAndCarter
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    AnalyticValidationTests.SchwarzschildHorizonRadius
    AnalyticValidationTests.SchwarzschildISCO
    AnalyticValidationTests.KerrHorizonRadius
    AnalyticValidationTests.KerrISCOPrograde
    AnalyticValidationTests.KerrISCODecreaseWithSpin
    AnalyticValidationTests.KerrErgosphereAtEquator
    AnalyticValidationTests.KerrErgosphereAtPole
    AnalyticValidationTests.KerrReducesToSchwarzschildAtZeroSpin
    AnalyticValidationTests.AsymptoticFlatness
    AnalyticValidationTests.SchwarzschildKretschmannScalar
    AnalyticValidationTests.KretschmannMonotonicDecrease
    AnalyticValidationTests.HorizonMetricDegeneracy
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    BeamPropagationTest.BeamInitialisation
    BeamPropagationTest.JacobianDeterminantFlat
    BeamPropagationTest.BeamGeometryExtraction
    BeamPropagationTest.OrientationDescribesOutputEllipseRatherThanInputBasis
    BeamPropagationTest.CausticDetection
    BeamPropagationTest.BeamIntegrationStep
    BeamPropagationTest.SchwarzschildRadialCongruenceMatchesClosedFormToOnePartPerMillion
    BeamPropagationTest.SchwarzschildCircularPhotonCongruenceMatchesClosedFormToOnePartPerMillion
    BeamPropagationTest.HorizonTermination
    BeamPropagationTest.ConversionRoundTrip
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    BloomFilterTests.DisabledBloomLeavesBufferUnchanged
    BloomFilterTests.ZeroIntensityBloomLeavesUnchanged
    BloomFilterTests.BrightPixelsCauseBloom
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    BuildGateAuthority.ReleaseReceiptBindsEveryInstalledProductAtInitialisation
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    CameraFactoryTest.CreatePinhole
    CameraFactoryTest.CreateThinLens
    CameraFactoryTest.CreateFisheye
    CameraFactoryTest.MalformedLensValueFailsClosed
    CameraFactoryTest.ConfigPassthrough
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    CameraWorldlineTest.RestScreenRayAndWorldlineComposeOverLensModels
    CameraWorldlineTest.ZeroVelocityIsExactlyRepresented
    CameraWorldlineTest.InvalidInternalWorldlineFailsClosed
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    CaptureSurfaceTests.SchwarzschildHorizonIsCaptured
    CaptureSurfaceTests.KerrOblateHorizonUsesKerrRadius
    CaptureSurfaceTests.HorizonlessSpacetimesNeverCapture
    CaptureSurfaceTests.CheckTerminationUsesCartesianNormAndCapture
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ChristoffelTests.FlatSpaceChristoffelAllZero
    ChristoffelTests.TorsionFreeSymmetry
    ChristoffelTests.SphericalGammaRThetaTheta
    ChristoffelTests.SphericalGammaThetaRTheta
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ColourGradingTests.IdentityParametersPreserveValues
    ColourGradingTests.ZeroSaturationProducesGrey
    ColourGradingTests.OutputClamped
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    CommandRouter.UnknownCommandReturnsFailure
    CommandRouter.NoCommandStillPrintsHelpSuccessfully
    CommandRouter.ExplicitMissingConfigPropagatesFailure
    CommandRouter.InfoRejectsEmptyAndSurplusArguments
    CommandRouter.ConfigRejectsMalformedSubcommandArguments
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    ConfigEnvironment.IntegerOverridesApplied
    ConfigEnvironment.MetricOverridesApplyIdentityAwareMassDefault
    ConfigEnvironment.BooleanOverrideParsed
    ConfigEnvironment.ColorModeOverrideApplied
    ConfigEnvironment.MalformedIntegerDeclines
    ConfigEnvironment.TrailingNumericGarbageDeclines
    ConfigEnvironment.NonFiniteNumberDeclines
    ConfigEnvironment.MalformedBooleanDeclines
    ConfigEnvironment.FailedOverrideLeavesConfigurationUnchanged
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    ConfigLoading.ExplicitMissingFileDeclines
    ConfigLoading.MalformedJsonDeclinesWithoutDefaults
    ConfigLoading.UnknownNestedFieldDeclines
    ConfigLoading.UnknownTopLevelFieldDeclines
    ConfigLoading.DuplicateFieldsDeclineInsteadOfUsingParserOrder
    ConfigLoading.ValidateCommandUsesTheSameStrictParserAsStartup
    ConfigLoading.InvalidKnownValueDeclines
    ConfigLoading.ValidPartialFileMergesOverDefaults
    ConfigLoading.SaveDeclinesWhenParentCannotBeCreated
    ConfigLoading.SaveDeclinesInvalidConfiguration
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    ConfigSchema.DefaultsRoundTripThroughJson
    ConfigSchema.LegacyFieldSpellingsParse
    ConfigSchema.PartialJsonKeepsDefaultsForOmittedFields
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ConfigValidation.DefaultConfigurationIsValid
    ConfigValidation.WidthBelowMinimumRejected
    ConfigValidation.NonPowerOfTwoTileSizeRejected
    ConfigValidation.UnknownMetricNameRejected
    ConfigValidation.KnownMetricAliasAccepted
    ConfigValidation.SpinAboveNearExtremalRejected
    ConfigValidation.RotatingLambdaRejected
    ConfigValidation.ObserverInsidePoleBufferRejected
    ConfigValidation.UnknownTonemapperRejected
    ConfigValidation.VulkanBackendNameAccepted
    ConfigValidation.UnknownBackendNameRejected
    ConfigValidation.CpuBackendAccepted
    ConfigValidation.NonFiniteValuesAreRejected
    ConfigValidation.UnknownOutputExtensionIsRejected
    ConfigValidation.UnrepresentedNarrowbandEmissionDeclines
    ConfigValidation.TurbulenceRequiresVolumeAndSpectralCoronaDeclines
    ConfigValidation.UnimplementedDenoiserRequestIsRejected
    ConfigValidation.ObserverAzimuthRangeIsValidated
    ConfigValidation.CameraWorldlineAndLensAreValidated
    ConfigValidation.FeatureSpecificControlsRequireTheirOwningModels
    ConfigValidation.MetricMassAndObserverCoordinateRadiusAreIdentityAware
    ConfigValidation.DeSitterRequestsEnforcePositiveLambdaAndSubNariaiBlackHole
    ConfigValidation.PolarisationRequiresRepresentedThinBlackHoleDisk
    ConfigValidation.MotionBlurAndWormholeTopologyHaveExplicitOperatorBoundaries
    ConfigValidation.DiskRequestDeclinesForEveryMetricWithoutAnEmissionModel
    ConfigValidation.RayBundlesRejectMetricsWithoutStationaryCurvatureTransport
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ConservationLawTests.NullConditionPreservedSchwarzschild
    ConservationLawTests.NullConditionPreservedKerr
    ConservationLawTests.KillingEnergyConservedSchwarzschild
    ConservationLawTests.KillingEnergyConservedKerr
    ConservationLawTests.KillingAngularMomentumConservedSchwarzschild
    ConservationLawTests.KillingAngularMomentumConservedKerr
    ConservationLawTests.AllConservedQuantitiesSchwarzschild
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Contracts.SatisfiedContractsPassSilently
    Contracts.AxiomIsNeverEvaluated
    Contracts.LanguageCapabilitiesDistinguishNativeFeaturesFromSubstitutes
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ContractsDeathTest.EnforcedPreconditionViolationTerminates
    ContractsDeathTest.EnforcedAssertionViolationTerminates
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ConversionTests.ToVec4dFromFloats
    ConversionTests.ToFloat4RoundTrip
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    CoordinateTransformTests.RoundTripBLCartesian_Equator
    CoordinateTransformTests.RoundTripBLCartesian_NearPole
    CoordinateTransformTests.RoundTripBLCartesian_FarField
    CoordinateTransformTests.RoundTripBLCartesian_NearHorizon
    CoordinateTransformTests.RoundTripKerrSchild_LowSpin
    CoordinateTransformTests.RoundTripKerrSchild_ModerateSpin
    CoordinateTransformTests.RoundTripKerrSchild_HighSpin
    CoordinateTransformTests.KerrDiskRadiusIsSpheroidalNotCylindrical
    CoordinateTransformTests.KerrRadiusPreservesExactZeroAndScaleCovariance
    CoordinateTransformTests.RoundTripKerrSchild_NearPole
    CoordinateTransformTests.JacobianDeterminant_Equator
    CoordinateTransformTests.JacobianDeterminant_MidLatitude
    CoordinateTransformTests.VectorTransformRoundTrip
    CoordinateTransformTests.OriginHandling
    CoordinateTransformTests.PhiWrapping
    CoordinateTransformTests.NorthPoleCoordinates
    CoordinateTransformTests.SouthPoleCoordinates
    CoordinateTransformTests.KerrOblateSpheroidal
    CoordinateTransformTests.KerrSolveR_Accuracy
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    CpuGeodesicReferenceTests.CPUBaselineSchwarzschildEscaping
    CpuGeodesicReferenceTests.CPUBaselineSchwarzschildHorizon
    CpuGeodesicReferenceTests.CPUBaselineKerrPrograde
    CpuGeodesicReferenceTests.ParityToleranceSpecification
    CpuGeodesicReferenceTests.ReferenceRayGridGeneration
    CpuGeodesicReferenceTests.ReferenceRaysCoverImpactParameters
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    DNGRParityTest.ExtremalKerrConfiguration
    DNGRParityTest.CameraDistance
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    DeterminismTests.DualArithmeticDeterminism
    DeterminismTests.SchwarzschildMetricDeterminism
    DeterminismTests.KerrMetricDeterminism
    DeterminismTests.ChristoffelDeterminism
    DeterminismTests.GeodesicIntegrationDeterminism
    DeterminismTests.InnerProductDeterminism
    DeterminismTests.NullNormalizationDeterminism
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    DiskCoordinateTest.CylindricalHeightAndPolarAngleRoundTripOnTheirOwnedDomain
    DiskCoordinateTest.AxisAndNonFiniteCoordinatesDeclineInsteadOfBecomingEquatorial
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    DispatchGovernor.FirstBandIsTheInitialHeightClampedToTheTile
    DispatchGovernor.BandsNeverExceedRemainingRowsNorDropBelowOne
    DispatchGovernor.GrowthPerStepIsBoundedByTheCap
    DispatchGovernor.ZeroMeasurementTakesTheCappedGrowthStep
    DispatchGovernor.OvershootShrinksProportionallyInOneStep
    DispatchGovernor.TruncatedTailBandFeedsBackOnlyItsOwnWork
    DispatchGovernor.LearnedAreaNormalisesAcrossBandWidths
    DispatchGovernor.DisabledControllerDispatchesWholeTilesAndIgnoresFeedback
    DispatchGovernor.TargetDefaultsWhenTheEnvironmentIsUnset
    DispatchGovernor.TargetHonoursTheOverrideIncludingZero
    DispatchGovernor.TargetFailsLoudOnGarbageNegativesAndNonFinite
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    DisplayBuffer.NonFiniteRadianceIsIdentifiedBeforeEncoding
    DisplayBuffer.MalformedDimensionsAndTilesFailClosed
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    DopplerToggleTest.SuppressionCollapsesDiskAsymmetry
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    DualNumberTests.AdditiveIdentity
    DualNumberTests.MultiplicativeIdentity
    DualNumberTests.AdditiveInverse
    DualNumberTests.AdditionCommutativity
    DualNumberTests.MultiplicationCommutativity
    DualNumberTests.AdditionAssociativity
    DualNumberTests.Distributivity
    DualNumberTests.Nilpotency
    DualNumberTests.ProductFormula
    DualNumberTests.DivisionFormula
    DualNumberTests.SinDerivative
    DualNumberTests.SinDerivativeMultipleAngles
    DualNumberTests.CosDerivative
    DualNumberTests.SqrtDerivative
    DualNumberTests.SqrtDerivativeMultiplePoints
    DualNumberTests.ChainRuleSinSquare
    DualNumberTests.ChainRuleSqrtSin
    DualNumberTests.ScalarMultiplication
    DualNumberTests.ScalarDivision
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    EXRRoundTripTests.HDRGradientSurvivesWriteAndRead
    EXRRoundTripTests.WriteFailsCleanlyOnBadPath
    EXRRoundTripTests.NonFiniteRadianceIsRejected
    EXRRoundTripTests.PpmRgbaBoundaryAppliesExactlyOneSrgbEncode
    EXRRoundTripTests.MalformedBufferShapesAreRejectedByEveryPublicWriter
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    Error.ExpectedPropagatesValue
    Error.FailCarriesDomainOperationAndDetail
    Error.DescriptionNamesDomainOperationAndDetail
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    FilmSimulationTest.DefaultAndPresetConfigsUseOnlyRepresentedControls
    FilmSimulationTest.MalformedControlFailsClosedAtPipelineBoundary
    FilmSimulationTest.DisabledGrainRequiresNeutralControls
    FilmSimulationTest.DisabledHalationRequiresNeutralControls
    FilmSimulationTest.InterstellarPresetSetsImplementedLookControls
    FilmSimulationTest.SpaceOdysseyPresetUsesLowerGrain
    FilmSimulationTest.Grain_AddsNoise
    FilmSimulationTest.Grain_DifferentFrames_DifferentNoise
    FilmSimulationTest.Grain_NonNegativeResult
    FilmSimulationTest.Halation_AffectsBrightAreas
    FilmSimulationTest.Halation_RedBias
    FilmSimulationTest.Vignette_DarkensCorners
    FilmSimulationTest.Vignette_CenterUnchanged
    FilmSimulationTest.Exposure_Positive_Brightens
    FilmSimulationTest.Exposure_Negative_Darkens
    FilmSimulationTest.Saturation_Zero_Grayscale
    FilmSimulationTest.OutputClamped_0_1
    FilmSimulationTest.Interstellar_Preset
    FilmSimulationTest.DisabledVignetteAndBloomRequireNeutralControls
    FilmSimulationTest.FullPipeline_DoesNotCrash
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    FisheyeCameraTest.CentreRayPointsInward
    FisheyeCameraTest.EdgeRayPerpendicularAt180Fov
    FisheyeCameraTest.ProjectionMaskedRayIsInactiveAndRepresented
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    GeodesicPathTests.PhotonSphereInstability
    GeodesicPathTests.WeakFieldLightDeflection
    GeodesicPathTests.RadialInfallApproachesHorizon
    GeodesicPathTests.KerrFrameDragging
    GeodesicPathTests.PhiWrapContinuity
    GeodesicPathTests.MetricSignatureCorrect
    GeodesicPathTests.KerrMetricSignature
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    GeodesicStateDTests.DefaultConstruction
    GeodesicStateDTests.ConstructionFromPosAndMom
    GeodesicStateDTests.ConservedQuantitiesSchwarzschildLike
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    GeodesicTests.FlatSpaceZeroAcceleration
    GeodesicTests.NormalizeNullPreservesCondition
    GeodesicTests.RadialNullVectorSatisfiesSchwarzschildForm
    GeodesicTests.TimelikeVectorNormalization
    GeodesicTests.LorentzBoostedObserver
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    GeodesicTracerRedshift.NearExtremalInnerDiskEmissionRemainsFinite
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    GeodesicTracerTest.Construction
    GeodesicTracerTest.CameraRayGeneration
    GeodesicTracerTest.BasicTracing
    GeodesicTracerTest.HorizonCapture
    GeodesicTracerTest.EscapeToInfinity
    GeodesicTracerTest.DiskIntersection
    GeodesicTracerTest.LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile
    GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile
    GeodesicTracerTest.LiveDiskCrossingCarriesTransportedPhysicalStokesOrientation
    GeodesicTracerTest.NoNumericalFailures
    GeodesicTracerTest.KerrMetricTracing
    GeodesicTracerTest.InvariantFrequencyTransferSelectsTheCircularEmitterBranch
    GeodesicTracerTest.TracingPerformance
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    GeodesicTracerVolumetric.TransferAccumulatesAcrossEveryTraversedSegment
    GeodesicTracerVolumetric.RedshiftAndDopplerReachTheLiveVolumeSource
    GeodesicTracerVolumetric.ProceduralTurbulenceAltersLiveTransferDeterministically
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    HamiltonianStateDTests.DefaultConstruction
    HamiltonianStateDTests.RoundTripViaGeodesicState
    HamiltonianStateDTests.ConstructFromGeodesicState
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ICameraTest.GetPositionReturnsConfigCoordinates
    ICameraTest.SetConfigUpdatesRayGeneration
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    KernelBeam.BeamFlagWiresDeviationWithoutMovingDefault
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    KernelParity.KerrSchildMetricMatchesLegacyToOnePartInMillion
    KernelParity.RepresentedSubThresholdKerrMetricIsScaleCovariant
    KernelParity.UnrepresentedKerrStageShrinksBeforeMetricEvaluation
    KernelParity.WormholeAndWarpMetricMatchLegacy
    KernelParity.KerrSchildChristoffelMatchesLegacyToOnePartInMillion
    KernelParity.FullPageThorneDiskTemperatureMatchesIndependentCoreModel
    KernelParity.RepresentedSmallKerrSpinDoesNotAliasToSchwarzschildIsco
    KernelParity.ShakuraSunyaevZeroTorqueTemperatureMatchesCoreModel
    KernelParity.TruncatedGaussianOpacityMatchesFiniteColumnCoreClosure
    KernelParity.BlackbodyMatchesIntegratedCoreSpectrum
    KernelParity.DiskEmissionAppliesExactlyOneGFourthFactor
    KernelParity.NearExtremalKerrLiveRenderIntegratorConservesEnergyAngularMomentumAndCarter
    KernelParity.PrecisionProbeArtifactsCarryOnlyTheirDeclaredFloat64Capability
    KernelParity.PrecisionRungsConserveNearExtremalKerrWithoutImageComparison
    KernelParity.BeamEllipseRetainsBothAxesAndOutputOrientation
    KernelParity.GeodesicDeviationIsFiniteAndCurvedNearBlackHole
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    KernelPortability.CudaEmissionCarriesTheNativeComputeEntryPoint
    KernelPortability.MetalEmissionCarriesTheNativeComputeEntryPoint
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow
    KernelTrace.Fp64RungAgreesWithFp32OnKerrScene
    KernelTrace.CompensatedRungTracksFp64AtLeastAsWellAsFp32
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    KerrMetricDTest.MetricInverseIdentity_Schwarzschild
    KerrMetricDTest.MetricInverseIdentity_Kerr
    KerrMetricDTest.ChristoffelSymmetry
    KerrMetricDTest.HamiltonianNullGeodesic
    KerrMetricDTest.HorizonRadius
    KerrMetricDTest.ExactExtremalParametersAreNotSilentlyRewritten
    KerrMetricDTest.BoyerLindquistAxisAndNonFinitePhaseSpaceDeclineWithoutClamping
    KerrMetricDTest.ISCORadius
    KerrMetricDTest.PhotonSphereRadius
    KerrMetricDTest.ChristoffelMatchesFiniteDifferencesOfMetric
    KerrMetricDTest.SecondDerivativesMatchFiniteDifference
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    KerrTests.MetricMatchesIndependentCartesianKerrSchildForm
    KerrTests.ComputedRadiusSatisfiesTheDefiningOblateQuartic
    KerrTests.CartesianMetricIsScaleCovariantBelowTheFormerSpinFloor
    KerrTests.SingularKerrDiskDeclinesInsteadOfReceivingAnEpsilonRadius
    KerrTests.HorizonsAndCaptureFollowTheIndependentKerrPolynomial
    KerrTests.ExactExtremalIscoIsNotClampedAndUndefinedDomainsDecline
    KerrTests.StaticLimitAuthorityIsAZeroOfProductionGtt
    KerrTests.ZeroSpinEqualsProductionSchwarzschild
    KerrTests.ClosedFormInverseMultipliesToIdentity
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    KottlerHorizonTests.ExactRootsSeparateBlackHoleCaptureFromCosmologicalHorizon
    KottlerHorizonTests.PureDeSitterHasOnlyACosmologicalHorizon
    KottlerHorizonTests.RegularDeSitterOriginUsesTheExactAnalyticLimit
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    LivePathConservationTests.SchwarzschildEnergyAndAngularMomentum
    LivePathConservationTests.KerrEnergyAndAngularMomentum
    LivePathConservationTests.NearExtremalKerrEnergyAngularMomentumAndCarter
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MandatoryKillingTests.SchwarzschildEnergyConservation
    MandatoryKillingTests.SchwarzschildAngularMomentumConservation
    MandatoryKillingTests.KerrEnergyConservation
    MandatoryKillingTests.KerrAngularMomentumConservation
    MandatoryKillingTests.ReissnerNordstromEnergyConservation
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Mat4dTest.Identity
    Mat4dTest.MatrixVectorMultiply
    Mat4dTest.Determinant
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Mat4dTests.DefaultIsZero
    Mat4dTests.IdentityDiagonal
    Mat4dTests.Trace
    Mat4dTests.IdentityDeterminant
    Mat4dTests.ZeroDeterminant
    Mat4dTests.KnownDeterminant
    Mat4dTests.MatrixVectorMultiplication
    Mat4dTests.ScalarMultiplication
    Mat4dTests.Addition
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MemoryGovernor.TwoGigabyteBudgetSeatsAWorkableTile
    MemoryGovernor.SmallerBudgetYieldsSmallerTile
    MemoryGovernor.TinyBudgetDeclinesLoudly
    MemoryGovernor.OverheadLargerThanUsableBudgetDeclines
    MemoryGovernor.TileNeverExceedsImageExtent
    MemoryGovernor.WorkingSetMatchesTheDerivedTile
    MemoryGovernor.EnvironmentOverrideResolvesBudget
    MemoryGovernor.MalformedOverrideDeclinesInsteadOfBorrowingTheDeviceBudget
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    MetricDerivativeTests.KerrSchildDerivativeSymmetry
    MetricDerivativeTests.KerrSchildFiniteDifferenceAgreement
    MetricDerivativeTests.KerrSchildNonZeroDerivatives
    MetricDerivativeTests.EllisDrainholeDerivativeSymmetry
    MetricDerivativeTests.EllisDrainholeRadialDerivative
    MetricDerivativeTests.EllisDrainholeAngularDerivative
    MetricDerivativeTests.EllisDrainholeFiniteDifferenceAgreement
    MetricDerivativeTests.AlcubierreDerivativeSymmetry
    MetricDerivativeTests.AlcubierreTimeDerivative
    MetricDerivativeTests.AlcubirrreFiniteDifferenceAgreement
    MetricDerivativeTests.AlcubierreOutsideBubble
    MetricDerivativeTests.ChristoffelSymmetryFromDerivatives
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MetricLoaderChainTests.AllMetricsHaveName
    MetricLoaderChainTests.AllMetricsHaveParameters
    MetricLoaderChainTests.AllMetricsEvaluateAtStandardPosition
    MetricLoaderChainTests.AllMetricsHaveLorentzianSignature
    MetricLoaderChainTests.AllMetricsAreSymmetric
    MetricLoaderChainTests.SchwarzschildMassParameterWorks
    MetricLoaderChainTests.KerrSpinParameterWorks
    MetricLoaderChainTests.RNChargeParameterWorks
    MetricLoaderChainTests.KerrReducesToSchwarzschild
    MetricLoaderChainTests.RNReducesToSchwarzschild
    MetricLoaderChainTests.FarFieldAsymptoticFlatness
    MetricLoaderChainTests.MinkowskiIsExactlyFlat
    MetricLoaderChainTests.HandlesNearPoles
    MetricLoaderChainTests.HandlesNearHorizon
    MetricLoaderChainTests.MetricDerivativesFinite
    MetricLoaderChainTests.MinkowskiZeroDerivatives
    MetricLoaderChainTests.SequentialEvaluationsStable
    MetricLoaderChainTests.DifferentPositionsDifferentMetrics
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MetricRegistryTests.EveryCanonicalNameParsesToItsOwnId
    MetricRegistryTests.EveryAliasParsesToItsOwnId
    MetricRegistryTests.ParsingIsCaseInsensitive
    MetricRegistryTests.UnknownNamesFailInsteadOfDefaulting
    MetricRegistryTests.MetricInfoRoundTripsById
    MetricRegistryTests.PositiveLambdaObserverAndHorizonShareOneCausalDomain
    MetricRegistryTests.KnownMetricNamesListsEveryCanonicalName
    MetricRegistryTests.FamilyDisplayNamesRoundTrip
    MetricRegistryTests.BackendSupportMatchesImplementations
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MetricValidationTests.KerrSchild_MetricSymmetry
    MetricValidationTests.KerrSchild_LorentzianSignature
    MetricValidationTests.KerrSchild_NoNaNInf
    MetricValidationTests.KerrSchild_MinkowskiLimit
    MetricValidationTests.KerrSchild_SchwarzschildWeakField
    MetricValidationTests.KerrMetricD_MetricSymmetry
    MetricValidationTests.KerrMetricD_InverseAccuracy
    MetricValidationTests.KerrMetricD_ChristoffelSymmetry
    MetricValidationTests.KerrMetricD_RiemannAntisymmetry
    MetricValidationTests.KerrMetricD_KretschmannScalar
    MetricValidationTests.KerrMetricD_RiemannMatchesFiniteDifferenceChristoffel
    MetricValidationTests.KerrMetricD_RiemannVacuumRicci
    MetricValidationTests.KerrMetricD_RiemannLoweredSymmetries
    MetricValidationTests.KerrMetricD_KretschmannMatchesRiemannContraction
    MetricValidationTests.KerrMetricD_HorizonRadius
    MetricValidationTests.KerrMetricD_ISCORadius
    MetricValidationTests.KerrMetricD_NoNaNInf
    MetricValidationTests.MorrisThorne_MetricSymmetry
    MetricValidationTests.MorrisThorne_FlareOutCondition
    MetricValidationTests.MorrisThorne_ThroatCondition
    MetricValidationTests.MorrisThorne_NoNaNInf
    MetricValidationTests.MorrisThorne_LorentzianSignature
    MetricValidationTests.SchwarzschildConsistency
    MetricValidationTests.KerrMetricD_DeterminantNonZero
    MetricValidationTests.KerrMetricD_NearHorizonStability
    MetricValidationTests.KerrMetricD_NearPoleStability
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MorrisThorneCartesianTests.ChartAgreementWithSphericalFamily
    MorrisThorneCartesianTests.DerivativesMatchFiniteDifferencesOfMetric
    MorrisThorneCartesianTests.AnalyticInverseIsExact
    MorrisThorneCartesianTests.ThroatIsTheCaptureSurface
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MorrisThorneTracerTest.CentralRayCapturedAtThroat
    MorrisThorneTracerTest.EdgeRayEscapes
    MorrisThorneTracerTest.DeflectionFallsQuadraticallyWithImpactParameter
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    MuellerMatrixTests.IdentityPreservesStokes
    MuellerMatrixTests.HorizontalPolariserOnUnpolarised
    MuellerMatrixTests.VerticalPolariserOnUnpolarised
    MuellerMatrixTests.CrossedPolariersExtinguish
    MuellerMatrixTests.QuarterWavePlateConvertsToCircular
    MuellerMatrixTests.CompositionAssociativity
    MuellerMatrixTests.DepolariserReducesPolarisation
    MuellerMatrixTests.HalfWavePlateFlipsHandedness
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    NaNInfDetectionTests.DualDivisionByZero
    NaNInfDetectionTests.DualSqrtNegative
    NaNInfDetectionTests.DualLogNonPositive
    NaNInfDetectionTests.DualSqrtNearZeroDerivative
    NaNInfDetectionTests.SchwarzschildMetricFinite
    NaNInfDetectionTests.SchwarzschildNearHorizon
    NaNInfDetectionTests.KerrMetricFinite
    NaNInfDetectionTests.KerrNearPoles
    NaNInfDetectionTests.TensorOperationsFinite
    NaNInfDetectionTests.ChristoffelFinite
    NaNInfDetectionTests.LargeValuesDoNotOverflow
    NaNInfDetectionTests.SmallValuesDoNotUnderflow
    PROPERTIES LABELS "Mandatory;Stability"
)

set_tests_properties(
    NumericalStabilityTest.NoNaNInMetric
    NumericalStabilityTest.DeterministicIntegration
    PROPERTIES LABELS "Mandatory;Stability"
)

set_tests_properties(
    NumericalStabilityTests.DualZeroHandling
    NumericalStabilityTests.DualDivisionSmallDenominator
    NumericalStabilityTests.SqrtNearZero
    NumericalStabilityTests.TrigLargeAngles
    NumericalStabilityTests.VectorSmallMagnitude
    NumericalStabilityTests.VectorLargeComponents
    NumericalStabilityTests.InnerProductNearNull
    NumericalStabilityTests.ChristoffelNearSingular
    PROPERTIES LABELS "Mandatory;Stability"
)

set_tests_properties(
    ObserverFrameTests.MovingKerrSchildObserverProducesAnOrthonormalTetrad
    ObserverFrameTests.InvalidWorldlineCannotBecomeASilentStaticObserver
    ObserverFrameTests.MovingCameraLaunchMatchesIndependentLorentzTransform
    ObserverFrameTests.KerrCameraRayIsPastNullAndLaunchFrequencyIsOne
    ObserverFrameTests.EulerianReferenceRemainsTimelikeInsideTheKerrErgosphere
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    OracleConnection.AnalyticConnectionAgreesWithMetricDerivatives
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    PNGWriterTest.WriteImageBuffer
    PNGWriterTest.WriteImageBufferRGBA
    PNGWriterTest.GammaCorrection
    PNGWriterTest.LargeImage
    PNGWriterTest.EmptyBuffer
    PNGWriterTest.ZeroSizeBuffer
    PNGWriterTest.NullPixels
    PNGWriterTest.NonFiniteRadianceIsRejected
    PNGWriterTest.MalformedBufferShapesAreRejected
    PNGWriterTest.DecodeRoundTripMatchesSRGBEncoding
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    ParallelTransportTests.ApplyPreservesIntensity
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    PinholeCameraTest.CentreRayPointsInward
    PinholeCameraTest.RayDirectionIsNormalised
    PinholeCameraTest.OriginMatchesConfig
    PinholeCameraTest.RightPixelIncreasesAzimuth
    PinholeCameraTest.UpPixelDecreasesTheta
    PinholeCameraTest.RayIsActiveAndRepresented
    PinholeCameraTest.TimeComponentIsZero
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    PixelSampling.EmitsExactlyEveryRequestedNonSquareCount
    PixelSampling.NonPositiveInputFailsClosed
    PixelSampling.NonSquareCountsCoverBothAxesWithoutRemainderBias
    PixelSampling.PatternIsDeterministic
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    PlatformPaths.PlatformNameIsNonEmpty
    PlatformPaths.ConfigSearchPathsAreOrderedAndNonEmpty
    PlatformPaths.AbsentResourceResolvesToNullopt
    PlatformPaths.ExecutableDirectoryIsResolved
    PlatformPaths.ExplicitResourceRootDisablesFallbacks
    PlatformPaths.ResourceNamesAndSymlinksCannotEscapeTheSelectedVolume
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    PolarisedEmissionTests.ChandrasekharAtmosphereHasPhysicalEndpointPolarisation
    PolarisedEmissionTests.ChandrasekharAtmospherePreservesHemisphericFlux
    PolarisedEmissionTests.ChandrasekharAtmosphereRejectsInvalidDirectionCosines
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    RK45IntegratorTests.DefaultConfigReasonable
    RK45IntegratorTests.OptimalStepIncreasesForSmallError
    RK45IntegratorTests.OptimalStepDecreasesForLargeError
    RK45IntegratorTests.OptimalStepRespectsBounds
    RK45IntegratorTests.MinkowskiStraightLine
    RK45IntegratorTests.MinkowskiNullConstraint
    RK45IntegratorTests.SchwarzschildIntegration
    RK45IntegratorTests.SchwarzschildNullConstraint
    RK45IntegratorTests.StepAdaptsToCurvature
    RK45IntegratorTests.RK45MatchesRK4Approximately
    RK45IntegratorTests.StepRejectionWorks
    RK45IntegratorTests.UnrepresentedStageShrinksBeforeMetricEvaluation
    RK45IntegratorTests.NoNaNInResults
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    RayBundleTest.FlatSpaceBundleMagnificationIsUnity
    RayBundleTest.KretschmannMatchesOracleSchwarzschild
    RayBundleTest.KretschmannMatchesOracleKerrEquatorial
    RayBundleTest.KretschmannMatchesOracleKerrOffEquatorial
    RayBundleTest.BundleFiniteAndDeterministicKerr
    RayBundleTest.MagnificationComesOnlyFromJacobiMap
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ReissnerNordstromTests.MetricMatchesIndependentCartesianKerrSchildForm
    ReissnerNordstromTests.AnalyticDerivativesMatchIndependentFiniteDifferences
    ReissnerNordstromTests.ZeroChargeEqualsProductionSchwarzschild
    ReissnerNordstromTests.HorizonAuthoritiesSatisfyTheIndependentPolynomial
    ReissnerNordstromTests.SuperExtremalDomainHasNoHorizonOrCaptureSurface
    ReissnerNordstromTests.ClosedFormInverseMultipliesToIdentity
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    RenderCommandParse.BasicFlagsMapToConfig
    RenderCommandParse.RepresentedVolumetricAndFilmFlagsSetEnables
    RenderCommandParse.FeatureParametersSelectTheirOwningModels
    RenderCommandParse.MotionBlurAndWormholeTopologyReachTheValidatedSchema
    RenderCommandParse.ExoticMetricScalesReachTheValidatedSchema
    RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent
    RenderCommandParse.BackendVulkanDeclinesMetricOffTheRenderPath
    RenderCommandParse.ReusedCommandDoesNotRetainAnEarlierGpuRequest
    RenderCommandParse.CliCpuOverridesLowerLayerVulkanBackend
    RenderCommandParse.UnknownMetricFailsValidation
    RenderCommandParse.UnknownOptionRejected
    RenderCommandParse.TrailingNumericGarbageRejected
    RenderCommandParse.NonFiniteNumericValueRejected
    RenderCommandParse.UnexpectedPositionalArgumentRejected
    RenderCommandParse.MetricOverrideDefaultsMassWithoutDiscardingExplicitInput
    RenderCommandParse.RetiredBackendNamesAreNotRemapped
    RenderCommandParse.MalformedCameraBetaRejected
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    RenderPipelineTests.LightrayInitialization
    RenderPipelineTests.LightrayPositionFinite
    RenderPipelineTests.AccelerationIsFinite
    RenderPipelineTests.SchwarzschildAccelerationNearHorizon
    RenderPipelineTests.IntegrateStepProducesValidState
    RenderPipelineTests.MultipleStepsProgressRay
    RenderPipelineTests.RK45StepSucceeds
    RenderPipelineTests.RK45ProducesFiniteValues
    RenderPipelineTests.DifferentMetricsGiveDifferentPaths
    RenderPipelineTests.HandlesNearPoleTheta
    RenderPipelineTests.HandlesLargeRadius
    RenderPipelineTests.IntegrationIsDeterministic
    RenderPipelineTests.KerrFrameDragging
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    RenderSessionProbe.TraceDomainScalesWithMassAndEnclosesTheObserver
    RenderSessionProbe.GeometricMetadataNeverInventsPhysicalLengthUnits
    RenderSessionProbe.FeatureSpecificControlsRequireOwnersAtTypedBoundary
    RenderSessionProbe.LowLevelTracerRejectsUnownedOrPartialFeatureControls
    RenderSessionProbe.FilmFinishPresetsRetainUnspecifiedPresetControls
    RenderSessionProbe.BackendAutoResolvesByDeviceRegistryAndCapabilities
    RenderSessionProbe.ConfigurationConversionPreservesObserverAndDiskControls
    RenderSessionProbe.InMemoryPreviewRejectsInactiveOutputPath
    RenderSessionProbe.CpuKerrRenderProducesValidPngAndExr
    RenderSessionProbe.CpuKerrRenderProducesValidPpmThroughTheOwnedWriter
    RenderSessionProbe.FilmAffectsDisplayOutputButNeverLinearExr
    RenderSessionProbe.StartIsAsynchronousAndCancellationIsTerminalWithoutOutput
    RenderSessionProbe.CompletionCallbackCanReenterLifecycleWithoutDeadlock
    RenderSessionProbe.PointStarfieldRejectsValuesItsGeneratorWouldClamp
    RenderSessionProbe.SceneEvidenceBindsCanonicalTypedConfiguration
    RenderSessionProbe.TypedNumericBoundariesMatchTheExternalConfigurationBoundary
    RenderSessionProbe.PolarisedAndTwoSheetRequestsDeclineAtTheTypedBoundary
    RenderSessionProbe.CpuPolarisationModeConsumesTransportedDiskStokes
    RenderSessionProbe.CpuMorrisThorneRenderCompletes
    RenderSessionProbe.EveryRegisteredCpuMetricCompletesAFrame
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    SchwarzschildTests.MetricMatchesIndependentCartesianKerrSchildForm
    SchwarzschildTests.AnalyticDerivativesMatchIndependentFiniteDifferences
    SchwarzschildTests.HorizonAndCaptureUseTheExactArealRadius
    SchwarzschildTests.FarFieldPerturbationHasExactInverseRadiusScaling
    SchwarzschildTests.StationaryClockRateMatchesTheExactStaticPotential
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Sha256Tests.EmptyStringMatchesNistKnownAnswer
    Sha256Tests.AbcMatchesNistKnownAnswer
    Sha256Tests.MultiBlockPaddingMatchesNistKnownAnswer
    Sha256Tests.FileStreamingMatchesTheSameKnownAnswer
    Sha256Tests.MillionByteFileMatchesNistKnownAnswerAcrossReadChunks
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ShadowBoundary.KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p
    ShadowBoundary.SchwarzschildCriticalImpactParameterMatchesAnalyticAt1080p
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SpectralEmissionTest.BlackbodyColourDirection
    SpectralEmissionTest.DopplerShiftDirection
    SpectralEmissionTest.TrueColorAppliesExactlyOneGFourthIntensityFactor
    SpectralEmissionTest.BolometricDiskAuthorityAppliesExactlyOneGFourthFactor
    SpectralEmissionTest.MotionBlurAveragesNonlinearTemporalRadiance
    SpectralEmissionTest.SpinDisplayFormat
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SpectralRadianceTest.BlackbodyPeakWavelength
    SpectralRadianceTest.BlackbodyWhitePoint
    SpectralRadianceTest.BlackbodyBinsMatchPlanckAuthorityAndRejectInvalidTemperature
    SpectralRadianceTest.RedshiftRebinsILambdaWithGFiveAndGFourBolometricScaling
    SpectralRadianceTest.RedshiftWavelengthShift
    SpectralRadianceTest.SRGBConversionRange
    SpectralRadianceTest.ACESConversion
    SpectralRadianceTest.SpectralArithmetic
    SpectralRadianceTest.WavelengthBinIndexing
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SpectralUtilsTests.PhysicalConstantsValid
    SpectralUtilsTests.PlanckRadiancePositive
    SpectralUtilsTests.PlanckRadianceIncreasesWithTemp
    SpectralUtilsTests.PlanckRadianceHandlesInvalid
    SpectralUtilsTests.WienPeakWavelengthReasonable
    SpectralUtilsTests.ColorMatchingPositiveInVisible
    SpectralUtilsTests.YPeaksNearGreen
    SpectralUtilsTests.WavelengthToXYZOutOfRange
    SpectralUtilsTests.XYZToRGBValid
    SpectralUtilsTests.SRGBGammaCorrection
    SpectralUtilsTests.LinearToSRGBClamps
    SpectralUtilsTests.BlackbodyColorTemperature
    SpectralUtilsTests.SunColorReasonable
    SpectralUtilsTests.DopplerFactorZeroVelocity
    SpectralUtilsTests.DopplerFactorReceding
    SpectralUtilsTests.DopplerFactorApproaching
    SpectralUtilsTests.TotalRedshiftCombined
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SpectralValidationTests.PlanckRadianceMatchesIndependentCodataEquation
    SpectralValidationTests.WienAuthorityMatchesNumericalPlanckMaximum
    SpectralValidationTests.StefanBoltzmannAuthorityMatchesIntegratedPlanckSpectrum
    SpectralValidationTests.DopplerFactorMatchesIndependentLorentzFormula
    SpectralValidationTests.TotalRedshiftComposesGravitationalAndDopplerFactors
    SpectralValidationTests.StaticMetricRedshiftHasThePhysicalDirection
    SpectralValidationTests.KerrDiskTransferMatchesKillingFieldContraction
    SpectralValidationTests.ZamoBranchRemainsTimelikeInsideTheErgosphere
    SpectralValidationTests.ComovingOpacityUsesInvariantAffinePathLength
    SpectralValidationTests.ObserverToSourceTransferPreservesForegroundEmissionOrder
    SpectralValidationTests.BlackbodyColourProgressionConsumesIntegratedSpectrum
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    StarEntryTests.ComputeColorProducesValidRGB
    StarEntryTests.ComputeColorUsesTheIntegratedBlackbodyAuthority
    StarEntryTests.HotStarIsBluer
    StarEntryTests.InvalidTemperatureFailsClosedInsteadOfDefaultingToSolar
    StarEntryTests.IntensityFromMagnitude
    StarEntryTests.IntensityMagnitudeRelation
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    StarfieldConfigTests.PointCatalogueProjectionHasNoUnconsumedSamplingControls
    StarfieldConfigTests.OutOfRangeStarCountFailsClosed
    StarfieldConfigTests.NonPositiveMinimumDistanceFailsClosed
    StarfieldConfigTests.UnorderedDistanceDomainFailsClosed
    StarfieldConfigTests.NarrowOrderedDistanceDomainIsPreservedExactly
    StarfieldConfigTests.OutOfRangeMagnitudeLimitFailsClosed
    StarfieldConfigTests.OutOfRangeApertureFailsClosed
    StarfieldConfigTests.NonFiniteScalarFailsClosedWithoutRewritingTheRequest
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    StarfieldGeneratorTests.GeneratesNonEmptyCatalog
    StarfieldGeneratorTests.SpatialIndexOwnsValidatedCatalogueSnapshot
    StarfieldGeneratorTests.CatalogSizeBounded
    StarfieldGeneratorTests.DirectionVectorsNormalised
    StarfieldGeneratorTests.AllTemperaturesPositive
    StarfieldGeneratorTests.AllDistancesPositive
    StarfieldGeneratorTests.DeterministicWithSameSeed
    StarfieldGeneratorTests.DifferentSeedsDifferentCatalogs
    StarfieldGeneratorTests.NoNaNInCatalog
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    StarfieldPointTest.CatalogueMeetsSizeFloorAndIsFinite
    StarfieldPointTest.CataloguePreservesTheRequestedCount
    StarfieldPointTest.BeamAccumulationFiniteAndNonConstant
    StarfieldPointTest.SpatialIndexMatchesExhaustiveBeamOracle
    StarfieldPointTest.EllipticalFootprintUsesBothAxesAndOrientation
    StarfieldPointTest.ImaxCatalogueIndexFitsTheTwoGigabyteOperatingEnvelope
    StarfieldPointTest.BeamFootprintSuppressesStarFlicker
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    StokesVectorTests.UnpolarisedConstruction
    StokesVectorTests.HorizontalPolarisation
    StokesVectorTests.VerticalPolarisation
    StokesVectorTests.CircularPolarisation
    StokesVectorTests.PhysicalConstraint
    StokesVectorTests.PolarisationDegree
    StokesVectorTests.LinearPolarisationDegree
    StokesVectorTests.CircularPolarisationDegree
    StokesVectorTests.EVPAComputation
    StokesVectorTests.IsPhysicalRejectsUnphysical
    StokesVectorTests.NormalisationProjection
    StokesVectorTests.ZeroIntensityHandling
    StokesVectorTests.ArithmeticOperators
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SymplecticIntegratorTest.HamiltonianConservation_10kSteps
    SymplecticIntegratorTest.EnergyConservation
    SymplecticIntegratorTest.AngularMomentumConservation
    SymplecticIntegratorTest.SymplecticStructurePreserved
    SymplecticIntegratorTest.HorizonTermination
    SymplecticIntegratorTest.EscapeDetection
    SymplecticIntegratorTest.KerrMetricIntegration
    SymplecticIntegratorTest.SingleStepNullCondition
    SymplecticIntegratorTest.IntegratorOrderComparison
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    TensorInverseTests.DiagonalMetricInverseIsExact
    TensorInverseTests.OffDiagonalMetricInverseIsExact
    TensorInverseTests.KerrSchildDeterminantIsMinusOne
    TensorInverseTests.KerrSchildNumericInverseSatisfiesIdentity
    TensorInverseTests.AnalyticInverseMatchesNumericInverse
    TensorInverseTests.ChristoffelSymbolsAreSymmetricInLowerIndices
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    TensorTests.VectorZeroInitialization
    TensorTests.VectorAddition
    TensorTests.VectorSubtraction
    TensorTests.VectorScalarMultiplication
    TensorTests.VectorLength
    TensorTests.VectorUnaryMinus
    TensorTests.MatrixZeroInitialization
    TensorTests.MinkowskiSymmetry
    TensorTests.MinkowskiDiagonal
    TensorTests.MinkowskiSignature
    TensorTests.InnerProductMinkowski
    TensorTests.InnerProductSpatial
    TensorTests.NullVectorInnerProduct
    TensorTests.Tensor3ZeroInitialization
    TensorTests.FlatSpaceChristoffelZero
    TensorTests.ChristoffelSymmetry
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ThinLensCameraTest.CentreRayWithDefaultSampling
    ThinLensCameraTest.DifferentSamplesGiveDifferentRays
    ThinLensCameraTest.CentreApertureSharesPinholeProjectionAndFieldOfView
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    TileScheduler.ReinitialiseResetsCompletionLedger
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    TonemapTests.ACESAtZero
    TonemapTests.ACESMonotone
    TonemapTests.ACESBounded
    TonemapTests.ReinhardAnalytic
    TonemapTests.FilmicAtZero
    TonemapTests.FilmicMonotone
    TonemapTests.FilmicBounded
    TonemapTests.ApplyExposureScaling
    TonemapTests.ApplyACESDispatch
    TonemapTests.ApplyReinhardDispatch
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    TurbulenceTest.MalformedProceduralDomainFailsClosedWithoutRewritingTheRequest
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Vec4dTest.Arithmetic
    Vec4dTest.IndexedAccess
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Vec4dTests.DefaultConstructionAllZeros
    Vec4dTests.ParameterisedConstruction
    Vec4dTests.IndexedAccess
    Vec4dTests.IndexedMutation
    Vec4dTests.Addition
    Vec4dTests.Subtraction
    Vec4dTests.ScalarMultiplication
    Vec4dTests.ScalarMultiplicationCommutative
    Vec4dTests.ScalarDivision
    Vec4dTests.CompoundAssignment
    Vec4dTests.Negation
    Vec4dTests.Norm2
    Vec4dTests.IsZero
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ViewCommandOperational.StrictParsingAndSessionProjection
    ViewCommandOperational.RelativisticJetsDeclineBeforeViewerInitialisation
    ViewCommandOperational.HeadlessRefinementProducesASynchronisedFrame
    ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames
    ViewCommandOperational.InputStateHandlesPressRepeatReleaseMouseAndScroll
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    VolumetricDiskClosure.TruncatedGaussianColumnEqualsDeclaredOpticalDepth
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    VulkanBackend.EnumerationReportsInsteadOfThrowing
    VulkanBackend.SlangKernelMatchesCpuReference
    VulkanBackend.WorkerThreadDispatchTearsDownSafely
    VulkanBackend.DeviceSelectionIsStrictAndRangeChecked
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    VulkanRenderSession.CapabilityBoundaryAcceptsRepresentedSceneSemantics
    VulkanRenderSession.CapabilityBoundaryRejectsUnrepresentedSceneSemantics
    VulkanRenderSession.ProceduralVolumetricTurbulenceReachesLiveKernel
    VulkanRenderSession.ThinAndVolumetricDopplerSuppressionAffectLiveEmission
    VulkanRenderSession.NonSquareMultisamplingCameraAndLensReachLiveKernel
    VulkanRenderSession.IndexedPointCatalogueReachesLiveKernel
    VulkanRenderSession.CombinedParitySceneRetainsResolvedImageStructure
    VulkanRenderSession.CpuVulkanPointCatalogueAgreeOnFlatScene
    VulkanRenderSession.ConstrainedBudgetDeclinesRatherThanChangingBackground
    VulkanRenderSession.Fp64RungRendersOrDeclinesLoudly
    VulkanRenderSession.CompensatedRungRendersOnAnyDevice
    VulkanRenderSession.Kerr160x120CompletesAcrossMultipleGovernedTiles
    VulkanRenderSession.CpuVulkanAgreeOnKerrGeometryWithinStatisticalBounds
    VulkanRenderSession.KerrNearExtremalBardeenBoundaryAt1080p
    VulkanRenderSession.CpuVulkanAgreeOnMorrisThorneGeometryWithinStatisticalBounds
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    WalkerPenrose.MalformedIntegratorStepDomainFailsClosed
    WalkerPenrose.BoyerLindquistInitialDataAndAxisExitDeclineWithoutSubstitution
    WalkerPenrose.ImaginaryPartEqualsKillingYanoContraction
    WalkerPenrose.PolarisationStaysOrthogonalAndNormalisedUnderTransport
    WalkerPenrose.ConstantConservedAlongKerrGeodesics
    WalkerPenrose.SchwarzschildEquatorialPerpendicularTransportsRigidly
    WalkerPenrose.FlatSpacePolarisationIsConstantInCartesian
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    WalkerPenroseLivePath.ConservesConstantAndOrthonormality
    WalkerPenroseLivePath.AgreesWithOracleAcrossCharts
    WalkerPenroseLivePath.TransportedVectorRotatesStokesByEvpa
    PROPERTIES LABELS "Mandatory;Correctness"
)

