# GENERATED FILE - do not edit by hand.
# Regenerate with: python3 scripts/generate-ctest-labels.py
# Labels policy lives in that script; CI fails on a stale copy.

set_tests_properties(
    AccretionDiskTest.ISCO_Schwarzschild
    AccretionDiskTest.ISCO_ExtremalKerr_Prograde
    AccretionDiskTest.ISCO_ExtremalKerr_Retrograde
    AccretionDiskTest.ISCO_ModerateSpin
    AccretionDiskTest.TemperatureProfileShape
    AccretionDiskTest.TemperatureScalingWithMass
    AccretionDiskTest.LimbDarkening
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
    AnalyticValidationTest.PhotonSphereSchwarzschildExact
    AnalyticValidationTest.PhotonSphereKerrPrograde
    AnalyticValidationTest.PhotonSphereKerrRetrograde
    AnalyticValidationTest.ISCOSchwarzschildExact
    AnalyticValidationTest.ISCOKerrPrograde
    AnalyticValidationTest.ISCONearExtremal
    AnalyticValidationTest.EinsteinRingSchwarzschildWeak
    AnalyticValidationTest.CriticalImpactParameterExact
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    AnalyticValidationTests.SchwarzschildHorizonRadius
    AnalyticValidationTests.SchwarzschildPhotonSphere
    AnalyticValidationTests.SchwarzschildISCO
    AnalyticValidationTests.SchwarzschildOrbitalVelocity
    AnalyticValidationTests.KerrHorizonRadius
    AnalyticValidationTests.KerrISCOPrograde
    AnalyticValidationTests.KerrISCODecreaseWithSpin
    AnalyticValidationTests.KerrErgosphereAtEquator
    AnalyticValidationTests.KerrErgosphereAtPole
    AnalyticValidationTests.WeakFieldDeflectionFormula
    AnalyticValidationTests.SolarDeflectionOrderOfMagnitude
    AnalyticValidationTests.KerrReducesToSchwarzschildAtZeroSpin
    AnalyticValidationTests.AsymptoticFlatness
    AnalyticValidationTests.SchwarzschildKretschmannScalar
    AnalyticValidationTests.KretschmannMonotonicDecrease
    AnalyticValidationTests.HorizonMetricDegeneracy
    AnalyticValidationTests.SchwarzschildISCOEnergy
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
    BlackbodyParityTest.GPUMatchesCPUAtReferenceTemperatures
    BlackbodyParityTest.ColourTemperatureOrdering
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    BloomFilterTests.DisabledBloomLeavesBufferUnchanged
    BloomFilterTests.ZeroIntensityBloomLeavesUnchanged
    BloomFilterTests.BrightPixelsCauseBloom
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    CameraAberrationTest.MatchesAnalyticFormulaAlongViewAxis
    CameraAberrationTest.ZeroBetaIsExactNoOp
    CameraAberrationTest.ForwardMotionBeamsTowardAxis
    CameraAberrationTest.ComposesOverLensModels
    PROPERTIES LABELS "Mandatory;Correctness"
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
    CaptureSurfaceTests.SchwarzschildHorizonIsCaptured
    CaptureSurfaceTests.KerrOblateHorizonUsesKerrRadius
    CaptureSurfaceTests.HorizonlessSpacetimesNeverCapture
    CaptureSurfaceTests.CheckTerminationUsesCartesianNormAndCapture
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ChristoffelBenchmark.SphericalVsCartesianPerformance
    ChristoffelBenchmark.PoleHandlingComparison
    PROPERTIES LABELS "Performance"
)

set_tests_properties(
    ChristoffelTests.FlatSpaceChristoffelAllZero
    ChristoffelTests.TorsionFreeSymmetry
    ChristoffelTests.SphericalGammaRThetaTheta
    ChristoffelTests.SphericalGammaThetaRTheta
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ChristoffelTransformTests.FlatSpaceReproducesTextbookSphericalSymbols
    ChristoffelTransformTests.KerrTransformedConnectionSymmetricAndFinite
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ColourGradingTests.IdentityParametersPreserveValues
    ColourGradingTests.ZeroSaturationProducesGrey
    ColourGradingTests.OutputClamped
    PROPERTIES LABELS "Correctness"
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
    ConfigEnvironment.MetricNameAndSpinOverridesApplied
    ConfigEnvironment.BooleanOverrideParsed
    ConfigEnvironment.ColorModeOverrideApplied
    ConfigEnvironment.MalformedIntegerDeclines
    ConfigEnvironment.TrailingNumericGarbageDeclines
    ConfigEnvironment.NonFiniteNumberDeclines
    ConfigEnvironment.MalformedBooleanDeclines
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
    ConfigValidation.VolumetricExtensionsRequireTheLiveVolumePath
    ConfigValidation.UnimplementedDenoiserRequestIsRejected
    ConfigValidation.ObserverAzimuthRangeIsValidated
    ConfigValidation.CameraWorldlineAndLensAreValidated
    ConfigValidation.MasslessMetricUsesAUnitObserverDistanceScale
    ConfigValidation.PolarisationRequiresRepresentedThinBlackHoleDisk
    ConfigValidation.MotionBlurAndWormholeTopologyHaveExplicitOperatorBoundaries
    ConfigValidation.DiskRequestDeclinesForEveryMetricWithoutAnEmissionModel
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
    CoronaConfigTests.ValidateClampsTemperature
    CoronaConfigTests.ValidateClampsOpticalDepth
    CoronaConfigTests.ValidateEnsuresOuterGreaterThanInner
    CoronaConfigTests.ComptonizationParameter
    CoronaConfigTests.ComptonizationHighOpticalDepth
    CoronaConfigTests.SpectralIndexFinite
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    CoronaEmissivityTests.ZeroOutsideBounds
    CoronaEmissivityTests.DisabledReturnsZero
    CoronaEmissivityTests.DecreasesWithRadius
    CoronaEmissivityTests.GaussianVerticalFalloff
    CoronaEmissivityTests.OpticalDepthDisabledIsZero
    CoronaEmissivityTests.OpticalDepthPositiveInsideCorona
    CoronaEmissivityTests.ScatteredIntensityZeroForZeroTau
    CoronaEmissivityTests.ComptonizedSourceLeavesOpticalDepthToTheRayMarcher
    CoronaEmissivityTests.ScatteredIntensityIncreasesWithTau
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    CoronaGeometryTests.DisabledReturnsNoContainment
    CoronaGeometryTests.OutsideRadialBoundsRejected
    CoronaGeometryTests.SlabGeometryEquatorialInside
    CoronaGeometryTests.SlabGeometryPolarOutside
    CoronaGeometryTests.LamppostNearSourceInside
    CoronaGeometryTests.SphereContainment
    CoronaGeometryTests.ExtendedScalesWithRadius
    PROPERTIES LABELS "Correctness"
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
    DNGRParityTest.InclinationAngle
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
    DisplayBuffer.MalformedDimensionsTilesAndGammaFailClosed
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
    EXRRoundTripTests.MalformedBufferShapesAreRejectedByEveryPublicWriter
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    EinsteinRingTest.EinsteinRingRadius
    EinsteinRingTest.CriticalImpactParameter
    EinsteinRingTest.MagnificationDefinition
    EinsteinRingTest.DeflectionAngleFormula
    EinsteinRingTest.MagnificationScaling
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    Error.ExpectedPropagatesValue
    Error.FailCarriesDomainOperationAndDetail
    Error.DescriptionNamesDomainOperationAndDetail
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    FPSBenchmarks.MinkowskiFPS
    FPSBenchmarks.SchwarzschildFPS
    FPSBenchmarks.SchwarzschildNearHorizonFPS
    FPSBenchmarks.KerrFPS
    FPSBenchmarks.KerrSpinComparison
    FPSBenchmarks.MetricComparison
    FPSBenchmarks.StepCountScaling
    PROPERTIES LABELS "Performance"
)

set_tests_properties(
    FPSThresholdTests.FrameBudgetCalculations
    FPSThresholdTests.ResolutionCalculations
    FPSThresholdTests.TargetsExceedBounds
    FPSThresholdTests.MinkowskiEvaluationTime
    FPSThresholdTests.SchwarzschildEvaluationTime
    FPSThresholdTests.KerrEvaluationTime
    FPSThresholdTests.MetricComplexityOrdering
    FPSThresholdTests.MinkowskiPerRayBudget
    FPSThresholdTests.SchwarzschildPerRayBudget
    FPSThresholdTests.KerrPerRayBudget
    FPSThresholdTests.PerStepBudgetAnalysis
    FPSThresholdTests.GPUParallelismAssumptions
    FPSThresholdTests.EvaluationConsistency
    FPSThresholdTests.NoNaNInCriticalPath
    PROPERTIES LABELS "Performance"
)

set_tests_properties(
    FilmSimulationTest.IMAX_AspectRatio_143
    FilmSimulationTest.IMAX5perf_AspectRatio_220
    FilmSimulationTest.ComputeHeight_Correct
    FilmSimulationTest.ComputeHeight_EvenDimension
    FilmSimulationTest.KodakVision3_500T_Settings
    FilmSimulationTest.KodakVision3_50D_LowerGrain
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
    FilmSimulationTest.DigitalClean_NoEffects
    FilmSimulationTest.FullPipeline_DoesNotCrash
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    FisheyeCameraTest.CentreRayPointsInward
    FisheyeCameraTest.EdgeRayPerpendicularAt180Fov
    FisheyeCameraTest.OutOfBoundsRayHasZeroWeight
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    GeodesicBenchmarks.CircularOrbitPeriod
    GeodesicBenchmarks.ISCORadius
    GeodesicBenchmarks.PhotonSphereRadius
    GeodesicBenchmarks.WeakFieldLightDeflection
    GeodesicBenchmarks.SolarDeflectionOrder
    GeodesicBenchmarks.GeodesicConservedQuantities
    GeodesicBenchmarks.RadialFallProperTime
    GeodesicBenchmarks.PerihelionPrecession
    GeodesicBenchmarks.VerificationDataPoints
    PROPERTIES LABELS "Performance"
)

set_tests_properties(
    GeodesicDeviationTests.MinkowskiRiemannZero
    GeodesicDeviationTests.DeviationConstantInFlatSpace
    GeodesicDeviationTests.SchwarzschildRiemannNonzero
    GeodesicDeviationTests.RiemannDecreasesWithDistance
    GeodesicDeviationTests.EllipseSemiAxesPositive
    GeodesicDeviationTests.CircularBundleInitialization
    GeodesicDeviationTests.BundleInitializationCorrect
    GeodesicDeviationTests.PropagationRemainsFinite
    PROPERTIES LABELS "Correctness"
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
    GeodesicTests.InnerProductConservation
    GeodesicTests.TimelikeVectorNormalization
    GeodesicTests.LorentzBoostedObserver
    GeodesicTests.CircularOrbitISCO
    GeodesicTests.PhotonSphereOrbit
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    GeodesicTracerTest.Construction
    GeodesicTracerTest.CameraRayGeneration
    GeodesicTracerTest.BasicTracing
    GeodesicTracerTest.HorizonCapture
    GeodesicTracerTest.EscapeToInfinity
    GeodesicTracerTest.DiskIntersection
    GeodesicTracerTest.DiskTemperatureProfile
    GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile
    GeodesicTracerTest.LiveDiskCrossingCarriesTransportedPhysicalStokesOrientation
    GeodesicTracerTest.NoNumericalFailures
    GeodesicTracerTest.KerrMetricTracing
    GeodesicTracerTest.GFactorDecompositionConsistency
    GeodesicTracerTest.MotionBlurConvergence
    GeodesicTracerTest.GFactorCoefficientNormalization
    GeodesicTracerTest.TracingPerformance
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    GeodesicTracerVolumetric.TransferAccumulatesAcrossEveryTraversedSegment
    GeodesicTracerVolumetric.TurbulenceAndCoronaAlterLiveTransferDeterministically
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
    JetDopplerTests.HeadOnApproachBoosted
    JetDopplerTests.RecedingDeBoosted
    JetDopplerTests.TransverseDirection
    JetDopplerTests.AnalyticFormula
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    JetEmissionTests.MagneticFieldAtLaunch
    JetEmissionTests.MagneticFieldDecays
    JetEmissionTests.ElectronDensityAtLaunch
    JetEmissionTests.ElectronDensityDecays
    JetEmissionTests.SynchrotronEmissivityPositive
    JetEmissionTests.BeamingApproachingBrighterThanReceding
    JetEmissionTests.PolarisationDegreeFormula
    JetEmissionTests.VelocityDirection
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    JetGeometryTests.InsideJetOnAxis
    JetGeometryTests.OutsideJetBelowLaunch
    JetGeometryTests.OutsideJetFarOffAxis
    JetGeometryTests.SouthernJet
    JetGeometryTests.JetRadiusMonotone
    JetGeometryTests.JetRadiusZeroBelowLaunch
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    JetRayMarchTests.EmissionOutsideJetIsZero
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    KernelBeam.BeamFlagWiresDeviationWithoutMovingDefault
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    KernelParity.KerrSchildMetricMatchesLegacyToOnePartInMillion
    KernelParity.WormholeAndWarpMetricMatchLegacy
    KernelParity.KerrSchildChristoffelMatchesLegacyToOnePartInMillion
    KernelParity.SymplecticStepMatchesLegacy
    KernelParity.YoshidaS4StepMatchesLegacy
    KernelParity.FullPageThorneDiskTemperatureMatchesIndependentCoreModel
    KernelParity.ChebyshevBlackbodyMatchesLegacy
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
    KerrMetricDTest.ISCORadius
    KerrMetricDTest.PhotonSphereRadius
    KerrMetricDTest.ChristoffelMatchesFiniteDifferencesOfMetric
    KerrMetricDTest.SecondDerivativesMatchFiniteDifference
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    KerrTests.ReducesToSchwarzschildAtZeroSpin
    KerrTests.ApproachesMinkowskiAtZeroMass
    KerrTests.MetricIsSymmetric
    KerrTests.ErgosphereGttPositive
    KerrTests.GttNegativeOutsideErgosphere
    KerrTests.FrameDraggingPresent
    KerrTests.FrameDraggingSign
    KerrTests.FrameDraggingVanishesAtPoles
    KerrTests.FrameDraggingAngularVelocity
    KerrTests.DeterminantFormula
    KerrTests.ChristoffelSymmetry
    KerrTests.ChristoffelNonZeroForRotatingBH
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    LivePathConservationTests.SchwarzschildEnergyAndAngularMomentum
    LivePathConservationTests.KerrEnergyAndAngularMomentum
    LivePathConservationTests.NearExtremalKerrEnergyAndAngularMomentum
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
    MemoryUsageTests.FundamentalTypeSizes
    MemoryUsageTests.DualNumberSize
    MemoryUsageTests.Vec4Size
    MemoryUsageTests.Metric4DSize
    MemoryUsageTests.PerRayMemory
    MemoryUsageTests.PerPixelMetricStorage
    MemoryUsageTests.NumericalMetricSmallGrid
    MemoryUsageTests.NumericalMetricMediumGrid
    MemoryUsageTests.NumericalMetricLargeGrid
    MemoryUsageTests.ChristoffelTextureMemory
    MemoryUsageTests.TotalVRAMAnalytic1080p
    MemoryUsageTests.TotalVRAMNumerical1080p
    MemoryUsageTests.VectorAllocationSafe
    MemoryUsageTests.MetricEvaluationNoLeak
    MemoryUsageTests.TestMemoryOverhead
    PROPERTIES LABELS "Performance"
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
    PROPERTIES LABELS "Correctness"
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
    NumericalMetricTests.MinkowskiADMReturnsMinkowski
    NumericalMetricTests.MinkowskiInverseMetric
    NumericalMetricTests.MetricIsSymmetric
    NumericalMetricTests.InverseMetricIsSymmetric
    NumericalMetricTests.MetricHasNegativeGtt
    NumericalMetricTests.SpatialMetricPositiveDefinite
    NumericalMetricTests.InverseMetricIdentity
    NumericalMetricTests.InverseMetricIdentitySchwarzschild
    NumericalMetricTests.SchwarzschildAsymptoticFlatness
    NumericalMetricTests.SchwarzschildLapseDecreasesTowardHorizon
    NumericalMetricTests.SchwarzschildConformalFactorIncreasesTowardHorizon
    NumericalMetricTests.ZeroShiftGivesDiagonalTimeComponents
    NumericalMetricTests.NonzeroShiftCreatesOffDiagonal
    NumericalMetricTests.LapseSquaredInGtt
    NumericalMetricTests.FourMetricDeterminant
    NumericalMetricTests.SphericalToCartesian
    NumericalMetricTests.SmallLapseHandled
    NumericalMetricTests.NoNaNInMetric
    NumericalMetricTests.ProperTimeReal
    NumericalMetricTests.NullVectorsExist
    PROPERTIES LABELS "Mandatory;Correctness"
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
    NumericalStabilityTests.SchwarzschildNearHorizon
    NumericalStabilityTests.PolarAxisSingularity
    NumericalStabilityTests.VectorSmallMagnitude
    NumericalStabilityTests.VectorLargeComponents
    NumericalStabilityTests.InnerProductNearNull
    NumericalStabilityTests.ChristoffelNearSingular
    NumericalStabilityTests.MachineEpsilonAddition
    NumericalStabilityTests.AccumulationError
    PROPERTIES LABELS "Mandatory;Stability"
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
    ParallelTransportTests.ZeroSpinNoRotation
    ParallelTransportTests.RotationIncreasesWithSpin
    ParallelTransportTests.ApplyPreservesIntensity
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    PinholeCameraTest.CentreRayPointsInward
    PinholeCameraTest.RayDirectionIsNormalised
    PinholeCameraTest.OriginMatchesConfig
    PinholeCameraTest.RightPixelIncreasesAzimuth
    PinholeCameraTest.UpPixelDecreasesTheta
    PinholeCameraTest.WeightIsOne
    PinholeCameraTest.TimeComponentIsZero
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    PixelSampling.EmitsExactlyEveryRequestedNonSquareCount
    PixelSampling.NonPositiveInputStillEmitsOneDefensiveSample
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
    PolarisedEmissionTests.SynchrotronPolarisationDegree
    PolarisedEmissionTests.SynchrotronEmissionIsPhysical
    PolarisedEmissionTests.ThomsonScatteringAt90Degrees
    PolarisedEmissionTests.ThomsonScatteringForward
    PROPERTIES LABELS "Correctness"
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
    RK45IntegratorTests.NoNaNInResults
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    RayBundleTest.FlatSpaceBundleMagnificationIsUnity
    RayBundleTest.KretschmannMatchesOracleSchwarzschild
    RayBundleTest.KretschmannMatchesOracleKerrEquatorial
    RayBundleTest.KretschmannMatchesOracleKerrOffEquatorial
    RayBundleTest.BundleFiniteAndDeterministicKerr
    RayBundleTest.MagnificationConsistencyWithScalarJacobian
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ReissnerNordstromTests.MetricIsDiagonal
    ReissnerNordstromTests.MetricIsSymmetric
    ReissnerNordstromTests.MetricHasLorentzianSignature
    ReissnerNordstromTests.SignatureForVariousCharges
    ReissnerNordstromTests.ReducesToSchwarzschildZeroCharge
    ReissnerNordstromTests.SchwarzschildLimitMultipleRadii
    ReissnerNordstromTests.AsymptoticFlatness
    ReissnerNordstromTests.ChargeEffectsFallOffFaster
    ReissnerNordstromTests.OuterHorizonFormula
    ReissnerNordstromTests.InnerHorizonFormula
    ReissnerNordstromTests.SchwarzschildHorizon
    ReissnerNordstromTests.NearExtremalHorizons
    ReissnerNordstromTests.EventHorizonGtt
    ReissnerNordstromTests.EventHorizonGrr
    ReissnerNordstromTests.MetricThetaTheta
    ReissnerNordstromTests.MetricPhiPhi
    ReissnerNordstromTests.MetricDeterminant
    ReissnerNordstromTests.DeterminantIsNegative
    ReissnerNordstromTests.ProperTimeStationaryObserver
    ReissnerNordstromTests.ChargeIncreasesRedshift
    ReissnerNordstromTests.PhotonSphereExists
    ReissnerNordstromTests.ChargeShrinksOuterHorizon
    ReissnerNordstromTests.ChargeExpandsInnerHorizon
    ReissnerNordstromTests.FunctionFCorrectness
    ReissnerNordstromTests.FPositiveOutsideHorizon
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    RenderCommandParse.BasicFlagsMapToConfig
    RenderCommandParse.RepresentedVolumetricAndFilmFlagsSetEnables
    RenderCommandParse.MotionBlurAndWormholeTopologyReachTheValidatedSchema
    RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent
    RenderCommandParse.BackendVulkanDeclinesMetricOffTheRenderPath
    RenderCommandParse.ReusedCommandDoesNotRetainAnEarlierGpuRequest
    RenderCommandParse.CliCpuOverridesLowerLayerVulkanBackend
    RenderCommandParse.UnknownMetricFailsValidation
    RenderCommandParse.UnknownOptionRejected
    RenderCommandParse.TrailingNumericGarbageRejected
    RenderCommandParse.NonFiniteNumericValueRejected
    RenderCommandParse.UnexpectedPositionalArgumentRejected
    RenderCommandParse.ExplicitMassOnMasslessMetricIsNotSilentlyDiscarded
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
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    RenderSessionProbe.CpuKerrRenderProducesValidPngAndExr
    RenderSessionProbe.FilmAffectsDisplayOutputButNeverLinearExr
    RenderSessionProbe.BackendAutoResolvesByDeviceRegistryAndCapabilities
    RenderSessionProbe.ConfigurationConversionPreservesObserverAndDiskControls
    RenderSessionProbe.StartIsAsynchronousAndCancellationIsTerminalWithoutOutput
    RenderSessionProbe.CompletionCallbackCanReenterLifecycleWithoutDeadlock
    RenderSessionProbe.PointStarfieldRejectsValuesItsGeneratorWouldClamp
    RenderSessionProbe.TypedNumericBoundariesMatchTheExternalConfigurationBoundary
    RenderSessionProbe.PolarisedAndTwoSheetRequestsDeclineAtTheTypedBoundary
    RenderSessionProbe.CpuPolarisationModeConsumesTransportedDiskStokes
    RenderSessionProbe.CpuMorrisThorneRenderCompletes
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    SchwarzschildTests.MetricIsDiagonal
    SchwarzschildTests.MetricIsSymmetric
    SchwarzschildTests.MetricHasLorentzianSignature
    SchwarzschildTests.AsymptoticFlatness
    SchwarzschildTests.EventHorizonGtt
    SchwarzschildTests.EventHorizonGrr
    SchwarzschildTests.MetricThetaTheta
    SchwarzschildTests.MetricPhiPhi
    SchwarzschildTests.MetricPhiPhiAtPoles
    SchwarzschildTests.ProperTimeStationaryObserver
    SchwarzschildTests.GravitationalRedshift
    SchwarzschildTests.MetricDeterminant
    SchwarzschildTests.DeterminantIsNegative
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    ShadowBoundary.KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p
    ShadowBoundary.SchwarzschildCriticalImpactParameterMatchesAnalyticAt1080p
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SpectralEmissionTest.NT_QFactorZeroAtISCO_Schwarzschild
    SpectralEmissionTest.NT_QFactorZeroAtISCO_Kerr
    SpectralEmissionTest.NT_QFactorPositiveOutsideISCO
    SpectralEmissionTest.NT_QFactorBounded
    SpectralEmissionTest.TemperatureRange
    SpectralEmissionTest.BlackbodyColourDirection
    SpectralEmissionTest.DopplerShiftDirection
    SpectralEmissionTest.TrueColorAppliesExactlyOneGFourthIntensityFactor
    SpectralEmissionTest.MotionBlurAveragesNonlinearTemporalRadiance
    SpectralEmissionTest.SpinDisplayFormat
    PROPERTIES LABELS "Mandatory;Correctness"
)

set_tests_properties(
    SpectralRadianceTest.BlackbodyPeakWavelength
    SpectralRadianceTest.BlackbodyWhitePoint
    SpectralRadianceTest.RedshiftEnergyConservation
    SpectralRadianceTest.RedshiftWavelengthShift
    SpectralRadianceTest.SRGBConversionRange
    SpectralRadianceTest.ACESConversion
    SpectralRadianceTest.SpectralArithmetic
    SpectralRadianceTest.WavelengthBinIndexing
    PROPERTIES LABELS "Correctness"
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
    SpectralUtilsTests.ApplyRedshiftEffect
    SpectralUtilsTests.ApplyBlueshiftEffect
    SpectralUtilsTests.TotalRedshiftCombined
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    SpectralValidationTests.PlanckFunctionBasic
    SpectralValidationTests.PlanckFunctionMonotonicity
    SpectralValidationTests.WienDisplacementLaw
    SpectralValidationTests.SolarPeakWavelength
    SpectralValidationTests.CandelaColor
    SpectralValidationTests.WienKnownValues
    SpectralValidationTests.RedshiftIntensityScaling
    SpectralValidationTests.RedshiftWavelengthTransform
    SpectralValidationTests.ColorTemperatureOrdering
    SpectralValidationTests.LimbDarkeningCoefficient
    SpectralValidationTests.BlackbodyColorProgression
    SpectralValidationTests.StefanBoltzmannLaw
    SpectralValidationTests.SolarLuminosity
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    StarEntryTests.ComputeColorProducesValidRGB
    StarEntryTests.HotStarIsBluer
    StarEntryTests.ZeroTemperatureDefaultsToSolar
    StarEntryTests.IntensityFromMagnitude
    StarEntryTests.IntensityMagnitudeRelation
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    StarfieldConfigTests.ValidateClampsStarCount
    StarfieldConfigTests.ValidateClampsMinDistance
    StarfieldConfigTests.ValidateEnsuresMaxGreaterThanMin
    StarfieldConfigTests.ValidateClampsMagnitudeLimit
    StarfieldConfigTests.ValidateClampsAperture
    PROPERTIES LABELS "Correctness"
)

set_tests_properties(
    StarfieldGeneratorTests.GeneratesNonEmptyCatalog
    StarfieldGeneratorTests.CatalogSizeBounded
    StarfieldGeneratorTests.DirectionVectorsNormalised
    StarfieldGeneratorTests.AllTemperaturesPositive
    StarfieldGeneratorTests.AllDistancesPositive
    StarfieldGeneratorTests.DeterministicWithSameSeed
    StarfieldGeneratorTests.DifferentSeedsDifferentCatalogs
    StarfieldGeneratorTests.NoNaNInCatalog
    PROPERTIES LABELS "Correctness"
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
    PROPERTIES LABELS "Correctness"
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
    ViewCommandOperational.HeadlessRefinementProducesASynchronisedFrame
    ViewCommandOperational.InputStateHandlesPressRepeatReleaseMouseAndScroll
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    VulkanBackend.EnumerationReportsInsteadOfThrowing
    VulkanBackend.SlangKernelMatchesCpuReference
    VulkanBackend.DeviceSelectionIsStrictAndRangeChecked
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
    VulkanRenderSession.CapabilityBoundaryAcceptsRepresentedSceneSemantics
    VulkanRenderSession.CapabilityBoundaryRejectsUnrepresentedSceneSemantics
    VulkanRenderSession.VolumetricTurbulenceAndCoronaReachLiveKernel
    VulkanRenderSession.NonSquareMultisamplingCameraAndLensReachLiveKernel
    VulkanRenderSession.IndexedPointCatalogueReachesLiveKernel
    VulkanRenderSession.CpuVulkanPointCatalogueAgreeOnFlatScene
    VulkanRenderSession.ConstrainedBudgetDeclinesRatherThanChangingBackground
    VulkanRenderSession.Fp64RungRendersOrDeclinesLoudly
    VulkanRenderSession.CompensatedRungRendersOnAnyDevice
    VulkanRenderSession.Kerr160x120CompletesAcrossMultipleGovernedTiles
    VulkanRenderSession.CpuVulkanAgreeOnKerrGeometryWithinStatisticalBounds
    VulkanRenderSession.KerrNearExtremalBardeenBoundaryAt1080p
    VulkanRenderSession.CpuVulkanAgreeOnMorrisThorneGeometryWithinStatisticalBounds
    VulkanRenderSession.BackendCompiledOrSkipped
    PROPERTIES LABELS "Mandatory;Operational"
)

set_tests_properties(
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

