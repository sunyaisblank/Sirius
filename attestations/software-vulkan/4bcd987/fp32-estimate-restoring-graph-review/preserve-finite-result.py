"""Preserve a terminal finite packet with an isolated Git index; never change HEAD."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

root = Path.cwd()
kind = sys.argv[1]
parent = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert parent == ('4498cc41168ad29fc4c111819a4cb7ccaa54a8bf' if kind in ('minimum-step-floor-alignment-stop', 'minimum-step-floor-focused-controls', 'divide-fp64-estimate-compiler-reader-stop', 'fp64-estimate-float-controls-properties', 'divide-fp64-estimate-reader-recovery', 'divide-fp64-estimate-ordinary-jit-screen', 'fp64-estimate-original-arithmetic', 'fp64-estimate-original-controls', 'fp64-estimate-original-benefit', 'fp64-estimate-final-retirement', 'original-cold-kerr-source-review', 'fp64-estimate-controls-profile-reader-stop', 'fp64-estimate-profile-recovery', 'divide-fp32-estimate-compiler-policy-stop', 'divide-fp32-estimate-policy-recovery', 'divide-fp32-estimate-graph-review', 'divide-fp32-estimate-ordinary-jit-screen', 'divide-fp32-estimate-restoring-compiler', 'divide-fp32-estimate-restoring-graph-review') else '27f5d9cdc6f9268b598a2b954febca8ec3468394' if kind == 'minimum-step-floor-configure-stop' else '46a9599a0cb30b0667e09838eac7a87690e982e4' if kind in ('multiply64-original-benefit', 'sanitizer-debug-o1-focused', 'multiply64-compiler-screen', 'multiply64-ordinary-jit-screen', 'multiply64-original-arithmetic', 'multiply64-profile-graph-stop', 'multiply64-o0-parser-fixture-stop', 'native-windows-system-5b5fad9', 'linux-gate-incomplete-38577826600', 'macos-full-ci-incomplete-38577826600', 'multiply64-original-controls', 'p1-kerr-observer-rss-stop', 'macos-current-provider-recovery', 'windows-runner-loss-38577826600', 'p2-cpu-radial-congruence', 'p2-source-sky-differential', 'e3-dispatch-continuations') else 'f5d649bf5067d1d68dc9215462f37830122f0540')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()

def pin(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}

def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

runner_source = parent
if kind in ('p4', 'kerr', 'volume'):
    label = {'p4': 'page-thorne', 'kerr': 'kerr-authority', 'volume': 'disk-volume'}[kind]
    case = {'p4': 'original-three-page-thorne-controls', 'kerr': 'original-four-kerr-authority-controls',
            'volume': 'original-four-disk-volume-controls'}[kind]
    work = root / ('out/science-p4-' + label)
    bundle = root / ('attestations/scientific-core/4bcd987/p4-' + label + '-current-source')
    tag = 'evidence/issue-53-p4-' + label + '-current-source'
    result = json.loads((work / 'result.json').read_text())
    assert result['passed'] and not any(result[k] for k in ('failures', 'errors', 'skips', 'disabled'))
    owner = json.loads((work / case / 'owner.json').read_text())
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
    mapping = [(work, 'files/out/science-p4-' + label),
               (root / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests', 'producer/sirius_oracle_tests'),
               (root / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    if kind in ('kerr', 'volume'):
        mapping += [(root / 'bin/linux-gcc/tests/core/sirius_core_tests', 'producer/sirius_core_tests'),
                    (root / 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake', 'producer/core-registrations.cmake')]
    if kind == 'volume':
        mapping += [(root / 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'producer/sirius_render_tests'),
                    (root / 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake', 'producer/render-registrations.cmake')]
    summary = result
    message = 'Preserve original current-source P4 ' + label + ' witnesses (#53)\n'
elif kind in ('e1', 'p1-e4', 'p2-schwarzschild-congruences', 'p2-cpu-boundary', 'p2-cpu-radial-congruence', 'p2-source-sky-differential', 'e3-dispatch-continuations', 'point-source-transfer', 'point-source-restricted-cells', 'point-source-catalogue-composition', 'point-source-response-remainder'):
    import xml.etree.ElementTree as ET
    if kind == 'e1':
        label, work, case = ('e1-field-equations', root / 'out/science-e1-field-equations', 'original-two-field-equation-controls')
        expected = {'KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell',
                    'KerrSchildFieldEquations.KerrNewmanPotentialIsSourceFreeOutsideTheRing'}
        producer, registrations = 'bin/linux-gcc/tests/core/sirius_core_tests', 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
    elif kind == 'p2-schwarzschild-congruences':
        label, work, case = ('p2-schwarzschild-congruences', root / 'out/science-p2-schwarzschild-congruences-current-source', 'original-two-p2-schwarzschild-congruence-controls')
        expected = {'BeamPropagationTest.SchwarzschildRadialCongruenceMatchesClosedFormToOnePartPerMillion',
                    'BeamPropagationTest.SchwarzschildCircularPhotonCongruenceMatchesClosedFormToOnePartPerMillion'}
        producer, registrations = 'bin/linux-gcc/tests/oracle/sirius_oracle_tests', 'bin/linux-gcc/tests/oracle/sirius_oracle_tests[1]_tests.cmake'
    elif kind == 'p2-source-sky-differential':
        label, work, case = ('p2-source-sky-differential', root / 'out/science-p2-source-sky-differential-current-source', 'original-six-source-sky-differential-controls')
        expected = {'SourceSkyDifferential.' + name for name in
                    ['KerrPositionAndTangentVariationsMatchIndependentDifferences',
                     'MovingFlatCameraMatchesAnalyticLorentzDirectionDerivative',
                     'NormalizesNonNullSpatialDirectionAndItsDerivative',
                     'DirectionIsInvariantUnderFrequencyAndLongitudinalVariation',
                     'ReseedsMetricVariationWithoutInheritedDualMetadata',
                     'DeclinesUnrepresentedFrequencyGeometryAndVariation']}
        producer, registrations = 'bin/linux-gcc/tests/core/sirius_core_tests', 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
    elif kind == 'e3-dispatch-continuations':
        label, work, case = ('e3-dispatch-continuations', root / 'out/science-e3-dispatch-continuations-current-source', 'original-eight-dispatch-continuation-controls')
        expected = {'DispatchGovernor.' + name for name in
                    ['ContinuationReadbackRejectsStaleOrUncommittedState',
                     'OnlyPhysicalTerminalReadbackCanFinalizeWithoutChangingTrajectory',
                     'ContinuationsPreserveFourCameraSamplesAndPublishOnce',
                     'LateContinuationOvershootReplaysChildrenFromBothZeroOrdinals',
                     'EveryContinuationTimingEnforcesThePhysicalSafetyBoundary',
                     'CancellationStopsBetweenContinuationsAndBeforePublication',
                     'CheapContinuationPhasesCannotHideThePeakSoftSizingFeedback',
                     'ContinuationExhaustionNeverImplicitlyCompletesOrAdvancesSample']}
        producer, registrations = 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake'
    elif kind == 'p2-cpu-radial-congruence':
        label, work, case = ('p2-cpu-radial-congruence', root / 'out/science-p2-cpu-radial-congruence-current-source', 'original-live-cpu-radial-congruence')
        expected = {'CpuJacobiOracle.RadialPointSourceCongruenceMatchesClosedForm'}
        producer, registrations = 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake'
    elif kind == 'p2-cpu-boundary':
        label, work, case = ('p2-cpu-boundary', root / 'out/science-p2-cpu-boundary-current-source', 'original-three-live-cpu-boundary-controls')
        expected = {'CpuTraceBoundary.CentralEventIsInvariantUnderBundleFeatureToggle',
                    'CpuTraceBoundary.JacobiBundleTerminatesAtTheSameCausalEvent',
                    'CpuTraceBoundary.OrdinaryEscapeClipsCentralAndJacobiAffineIntervals'}
        producer, registrations = 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake'
    elif kind == 'point-source-response-remainder':
        label, work, case = ('point-source-response-remainder', root / 'out/science-point-source-response-remainder-current-source', 'original-nine-cpu-point-source-response-controls')
        expected = {'PointSourceResponseTest.' + name for name in
                    ['TracedImageDensityKeepsTheOriginalGaussianAndSignedMeasure',
                     'IntegratedDetectorFluxFollowsSignedAffineLensMeasure',
                     'UnresolvedImageShapeStaysInDetectorCoordinates',
                     'EllipticalSupportRejectsInsideMajorCircleOutsideEllipse',
                     'SingularAndBroadMapsRequestRefinementWithoutFloors',
                     'InvalidArithmeticIsDistinctFromEmptyContribution',
                     'CatalogueAccumulationPreservesFluxColorAndExactSupport',
                     'CompositionCancellationAndUnrepresentedParityDecline',
                     'MovingThinLensAndSourceMapComposeIntoReceivedFlux']}
        producer, registrations = 'bin/linux-gcc/tests/core/sirius_core_tests', 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
    elif kind == 'point-source-catalogue-composition':
        label, work, case = ('point-source-catalogue-composition', root / 'out/science-point-source-catalogue-composition-current-source', 'original-three-cpu-catalogue-composition-controls')
        expected = {'PointSourceResponseTest.CatalogueBandTransferComposesDensityAndMovingFrequency',
                    'PointSourceResponseTest.CatalogueBandTransferKeepsFluxAndDensityLinear',
                    'PointSourceResponseTest.CatalogueBandTransferDeclinesInvalidFrequencyAndPartialFailure'}
        producer, registrations = 'bin/linux-gcc/tests/core/sirius_core_tests', 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
    elif kind == 'point-source-restricted-cells':
        label, work, case = ('point-source-restricted-cells', root / 'out/science-point-source-restricted-cells-current-source', 'original-three-cpu-restricted-cell-controls')
        expected = {'PointSourceResponseTest.RestrictedCellsPreserveOriginalGaussianAndNormalization',
                    'PointSourceResponseTest.RestrictedCellsOwnSharedEdgesAndKeepParentCutoff',
                    'PointSourceResponseTest.RestrictedCellsDeclineInvalidOrUnrepresentedGeometry'}
        producer, registrations = 'bin/linux-gcc/tests/core/sirius_core_tests', 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
    elif kind == 'point-source-transfer':
        label, work, case = ('point-source-transfer', root / 'out/science-point-source-transfer-current-source', 'original-six-cpu-point-source-transfer-controls')
        expected = {'PointSourceTransfer.ReferenceFrequencyPreservesCatalogueColourAndFlux',
                    'PointSourceTransfer.ShiftedSpectrumMatchesIndependentFixedBandPlanckLaw',
                    'PointSourceTransfer.RejectsReciprocalExtraGainAndShiftedNormalizerSubstitutions',
                    'PointSourceTransfer.MovingObserverUsesBandGainAndAngularAreaExactlyOnce',
                    'PointSourceTransfer.ShiftedTemperatureHasNoCatalogueLookupClamp',
                    'PointSourceTransfer.InvalidOrOverflowingTransferHasNoPublishedValue'}
        producer, registrations = 'bin/linux-gcc/tests/core/sirius_core_tests', 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
    else:
        label, work, case = ('p1-e4-oracle', root / 'out/science-p1-e4-oracle-current-source', 'original-three-p1-e4-oracle-controls')
        expected = {'AnalyticValidationTest.NearExtremalKerrConservesEnergyAngularMomentumAndCarter',
                    'AnalyticValidationTest.SymplecticChartExitAndNullProjectionFailureRemainDistinct',
                    'OracleConnection.AnalyticConnectionAgreesWithMetricDerivatives'}
        producer, registrations = 'bin/linux-gcc/tests/oracle/sirius_oracle_tests', 'bin/linux-gcc/tests/oracle/sirius_oracle_tests[1]_tests.cmake'
    bundle = root / ('attestations/scientific-core/4bcd987/' + label + '-current-source')
    tag = 'evidence/issue-53-' + label + '-current-source'
    result = json.loads((work / 'result.json').read_text())
    owner = json.loads((work / case / 'owner.json').read_text())
    assert result['passed'] and not any(result[k] for k in ('failures', 'errors', 'skips', 'disabled'))
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
    document = ET.parse(work / case / 'ctest.xml').getroot()
    cases = document.findall('.//testcase')
    assert len(cases) == len(expected) and {row.get('name') for row in cases} == expected
    assert all(row.get('status', 'run') == 'run' and not any(row.findall(k) for k in ('failure', 'error', 'skipped')) for row in cases)
    for suite in [document, *document.findall('.//testsuite')]:
        assert all(int(suite.get(k, '0')) == 0 for k in ('failures', 'errors', 'skipped', 'disabled'))
    mapping = [(work, 'files/out/' + work.name),
               (root / producer, 'producer/' + Path(producer).name),
               (root / registrations, 'producer/registrations.cmake')]
    summary = {'schema': 'sirius-original-finite-science-controls-v1', 'result': result, 'owner': owner,
               'full_criterion_complete': False, 'qualification_claimed': False}
    message = 'Preserve original finite ' + label + ' controls (#53)\n'
elif kind == 'sanitizer-debug-o1-source':
    work = root / 'out/sanitizer-debug-o1'
    source = json.loads((work / 'source-commit.json').read_text())
    runner_source = source['commit']
    assert runner_source == '46a9599a0cb30b0667e09838eac7a87690e982e4'
    assert source['base'] == '5b5fad91c08c0533e29c4f5c7ae249725162e405'
    assert subprocess.check_output(['git', 'show', runner_source + ':CMakePresets.json']) == (work / 'CMakePresets.json').read_bytes()
    assert subprocess.check_output(['git', 'diff-tree', '--no-commit-id', '--name-only', '-r', source['base'], runner_source], text=True).splitlines() == ['CMakePresets.json']
    old, new = [json.loads((work / name).read_text()) for name in ['CMakePresets.base.json', 'CMakePresets.json']]
    preset = next(row for row in new['configurePresets'] if row['name'] == 'linux-gcc-sanitize')
    for key in ['CMAKE_C_FLAGS_DEBUG', 'CMAKE_CXX_FLAGS_DEBUG']:
        assert preset['cacheVariables'].pop(key) == '-g -O1'
    assert new == old
    bundle = root / 'attestations/qualification/46a9599/sanitizer-debug-o1-source'
    tag = 'evidence/issue-48-sanitizer-debug-o1-source'
    mapping = [(work, 'files/out/sanitizer-debug-o1')]
    summary = {'schema': 'sirius-sanitizer-debug-o1-reviewed-source-v1', 'source': source,
               'static_checks': json.loads((work / 'checks.json').read_text()),
               'change': 'Only existing sanitizer Debug C/CXX flags add -O1. Debug assertions, enforce contracts, all instrumentation/environments, workloads, deadlines and other presets remain unchanged.',
               'configure_build_runtime_benefit_deadline_or_qualification_validated': False,
               'original_failed_Vulkan_retry_performed': False,
               'qualification_claimed': False}
    message = 'Preserve reviewed sanitizer Debug O1 source candidate without runtime claims (#48)\n'
elif kind == 'windows-runner-loss-38577826600':
    runner_source = '5b5fad91c08c0533e29c4f5c7ae249725162e405'
    work = root / 'out/native-provider-integration/full-ci-38577826600/windows-build'
    bundle = root / 'attestations/qualification/5b5fad9/windows-runner-loss-38577826600'
    tag = 'evidence/issue-48-windows-runner-loss-38577826600'
    job = json.loads((work / 'job-api.json').read_text())
    check = json.loads((work / 'check-run-api.json').read_text())
    annotations = json.loads((work / 'check-run-annotations-api.json').read_text())
    summary = json.loads((work / 'analysis.json').read_text())
    assert job['id'] == check['id'] == 115899500048 and job['run_id'] == 38577826600
    assert job['head_sha'] == check['head_sha'] == runner_source
    assert job['status'] == check['status'] == 'completed' and job['conclusion'] == check['conclusion'] == 'failure'
    assert len(annotations) == 1 and 'hosted runner lost communication with the server' in annotations[0]['message']
    assert json.loads((work / 'independent-review.json').read_text())['passed']
    assert not summary['underlying_cause_established'] and not summary['raw_complete_job_log_obtained']
    assert not summary['full_job_diagnostics_or_attestation_uploaded'] and not summary['qualification_claimed']
    assert pin(work / 'job-log-endpoint-response.txt') == {'bytes': 215, 'sha256': 'cdea14bcd09918756094706dc7e1371da88ede6fc214b0c57ececf5456bc3de8'}
    mapping = [(work, 'files/out/native-provider-integration/full-ci-38577826600/windows-build')]
    message = 'Preserve Windows hosted-runner communication loss and unavailable-log evidence (#48)\n'
elif kind == 'linux-sanitizer-38577826600':
    import zipfile
    runner_source = '5b5fad91c08c0533e29c4f5c7ae249725162e405'
    work = root / 'out/native-provider-integration/full-ci-38577826600/linux-sanitizers'
    bundle = root / 'attestations/qualification/5b5fad9/linux-sanitizer-incomplete-38577826600'
    tag = 'evidence/issue-48-linux-sanitizer-incomplete-38577826600'
    job = json.loads((work / 'job-api.json').read_text())
    result = json.loads((work / 'analysis.json').read_text())
    assert job['head_sha'] == runner_source and job['status'] == 'completed' and job['conclusion'] == 'cancelled'
    assert result['passed_cases'] == 23 and len(result['timed_out_cases']) == 2
    assert not result['final_JUnit_or_Mandatory_pass_present'] and not result['sanitizer_or_runtime_or_release_qualification_claimed']
    assert pin(work / 'diagnostics.zip') == {'bytes': 46116, 'sha256': '38cb874752bb03f23f2ae876756b899cb7478cb0e83fc136aaa0cb21e90bd9a2'}
    metadata = json.loads((work / 'run-artifacts-api.json').read_text())
    artifact = next(row for row in metadata['artifacts'] if row['id'] == 11915007461)
    assert artifact['digest'] == 'sha256:38cb874752bb03f23f2ae876756b899cb7478cb0e83fc136aaa0cb21e90bd9a2'
    with zipfile.ZipFile(work / 'diagnostics.zip') as archive:
        assert archive.testzip() is None
        assert {row.filename: {'bytes': row.file_size, 'sha256': hashlib.sha256(archive.read(row)).hexdigest()}
                for row in archive.infolist()} == result['artifact']['complete_raw_ZIP_members']
        assert archive.read('bin/linux-gcc-sanitize/generated/sirius/mandatory_gate_ctest.log') == (work / 'mandatory-ctest.log').read_bytes()
    mapping = [(work, 'files/out/native-provider-integration/full-ci-38577826600/linux-sanitizers')]
    summary = result
    message = 'Preserve incomplete sanitizer estate, original CPU and cold Vulkan timeouts, and full raw diagnostics (#48)\n'
elif kind == 'native-failure':
    work = root / 'out/macos-moltenvk-integration'
    bundle = root / 'attestations/native-vulkan/4bcd987/macos-provider-discovery-f5d649b'
    tag = 'evidence/issue-48-provider-discovery-f5d649b'
    archive = work / 'native-38556777741.zip'
    assert pin(archive) == {'bytes': 75835745, 'sha256': '2b0dbcd768b646a9d258c1942c797e3a642a96ece71a4f2c2d51e2da00acf030'}
    import zipfile
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        prefix = 'out/macos-moltenvk-integration/observation/'
        observation = json.loads(zipped.read(prefix + 'integration.json'))
        assert observation['runner_source']['runner_revision'] == parent
        assert observation['provider_source_inputs_unchanged'] and not observation['cases']
        assert not any(name.startswith(prefix + 'provider-build/') for name in zipped.namelist())
        owner = json.loads(zipped.read(prefix + 'toolchain/owner.json'))
        assert owner['status'] == 'terminal' and owner['child_exit'] == 1 and owner['stop_reason'] is None
        assert owner['whole_inputs_unchanged'] and owner['containment_complete'] and owner['observation_complete']
        assert not owner['cleanup_errors'] and all(row['absent'] for row in owner['observed_births'])
        assert b'CMake module directory is ambiguous' in zipped.read(prefix + 'toolchain/stderr.log')
    mapping = [(archive, 'artifact.zip'), (work / 'run-38556777741-failure.log', 'workflow-failure.log')]
    summary = {'schema': 'sirius-failed-native-provider-discovery-v1', 'run': 38556777741,
               'job': 115827784818, 'artifact': 11890208325, 'archive': pin(archive),
               'runner_revision': parent, 'provider_source_commit': 'e2422dd92017879a452377a1122f1838e6610b3c',
               'failure': 'CMake module-directory layout guessed incorrectly; query its authoritative CMAKE_ROOT instead.',
               'observation': observation, 'toolchain_owner': owner, 'provider_build_started': False,
               'qualification_claimed': False, 'raw_scope': 'Complete immutable native artifact ZIP plus source snapshots and workflow failure log. Remote-only system input pins are observations.'}
    message = 'Preserve failed native CMake discovery and complete raw artifact (#48)\n'
elif kind == 'native-four-controls':
    import xml.etree.ElementTree as ET
    import zipfile
    runner_source = '5b22bee3f15d628998624aceb9ede6e9088a6401'
    work = root / 'out/macos-moltenvk-integration'
    archive = work / 'native-38567408322.zip'
    assert pin(archive) == {'bytes': 105870928, 'sha256': '3e235b95061eb4e534401b66d24865b3b8799b87003faa2a3075be0d3d5ceb68'}
    bundle = root / 'attestations/native-vulkan/4bcd987/macos-custom-provider-four-controls'
    tag = 'evidence/issue-48-custom-provider-four-controls'
    names = ['PhysicalInitializationRetainsTheHamiltonianResidual',
             'ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation',
             'JointRkStagesRetainCriticalIncrementsAndEmbeddedError',
             'DenseSegmentsPreserveSmallCovariantArrivalDerivatives']
    phases = ['toolchain', 'stock-inventory', 'stock-uuid', 'provider-build', 'diagnostic-inventory', 'diagnostic-uuid']
    phases += [str(i + 1) + '-' + name for i, name in enumerate(names)]
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        prefix = 'out/macos-moltenvk-integration/observation/'
        observation = json.loads(zipped.read(prefix + 'integration.json'))
        assert observation['runner_source']['runner_revision'] == runner_source
        assert observation['status'] == 'four_original_controls_passed_custom_provider'
        assert all(observation[k] is True for k in ('provider_identity_accepted', 'all_four_controls_passed', 'numerics_executed', 'provider_source_inputs_unchanged'))
        assert not observation['qualification_claimed'] and not observation['production_adopted']
        assert [row['case'] for row in observation['cases']] == ['RetainedComputeTest.' + name for name in names]
        assert all(row['owner_complete'] and row['test']['passed'] and row['exit'] == 0 and row['shader_capture_bound'] for row in observation['cases'])
        identities = json.loads(zipped.read(prefix + 'provider-identities.json'))
        stock, diagnostic = identities['stock'], identities['diagnostic']
        assert stock['device']['pipeline_cache_uuid'] == 'db445ff20f0702050000000100000000'
        assert diagnostic['device']['pipeline_cache_uuid'] == '76cb8b210f0702050000000100000000'
        assert diagnostic['device']['driver_info'] == '1.4.1+sirius-issue48-unforced.76cb8b212bc049261146b526d75ebddafe017b18fe1f4c2c7065c3024c4ffe9e'
        assert diagnostic['device']['name'] == stock['device']['name'] == 'Apple Paravirtual device'
        assert diagnostic['device']['api_version'] == stock['device']['api_version'] == 4206926
        assert diagnostic['loader'] == stock['loader']
        owners = {}
        for phase in phases:
            folder = prefix + phase
            owner = json.loads(zipped.read(folder + '/owner.json'))
            assert owner['status'] == 'terminal' and owner['child_exit'] == 0 and owner['stop_reason'] is None
            assert owner['whole_inputs_unchanged'] and owner['containment_complete'] and owner['observation_complete']
            assert not owner['cleanup_errors'] and not owner['remaining_observed']
            assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
            assert json.loads(zipped.read(folder + '/inputs-before.json')) == json.loads(zipped.read(folder + '/inputs-after.json'))
            owners[phase] = owner
        for i, name in enumerate(names):
            folder = prefix + str(i + 1) + '-' + name
            document = ET.fromstring(zipped.read(folder + '/gtest.xml'))
            cases = document.findall('.//testcase')
            assert len(cases) == 1 and cases[0].get('classname') + '.' + cases[0].get('name') == 'RetainedComputeTest.' + name
            assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
            assert not any(cases[0].findall(k) for k in ('failure', 'error', 'skipped'))
            for suite in [document, *document.findall('.//testsuite')]:
                assert all(int(suite.get(k, '0')) == 0 for k in ('failures', 'errors', 'skipped', 'disabled'))
            joined = json.loads(zipped.read(folder + '/provider-shader-joins.json'))
            assert joined['loaded']['libMoltenVK'] == diagnostic['driver'] and joined['loaded']['libvulkan'] == diagnostic['loader']
            assert joined['shaders']['generic']['all_dumped_spirv_bound']
            for filename, row in joined['shaders']['generic']['files'].items():
                raw = zipped.read(folder + '/shaders/' + filename)
                assert {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()} == {k: row[k] for k in ('bytes', 'sha256')}
            assert any(joined['shaders']['expected_module'] in row['embedded_aliases'] for row in joined['shaders']['pairs'])
            if i == 0:
                assert joined['shaders']['pairs'][0]['exact_initialize_control']
                assert joined['shaders']['pairs'][0]['msl']['sha256'] == '80b33059934cc77a685c5738b9413c94379e45f1b18d0e6635a6432ff0407921'
        binary = zipped.read('bin/macos/macos-moltenvk-integration/MoltenVK/libMoltenVK.1.4.1.dylib')
        assert {'bytes': len(binary), 'sha256': hashlib.sha256(binary).hexdigest()} == {k: diagnostic['driver'][k] for k in ('bytes', 'sha256')}
    mapping = [(work / 'run-38567408322-success.log', 'workflow-success.log')]
    summary = {'schema': 'sirius-original-native-custom-provider-controls-v1', 'run': 38567408322, 'job': 115863521359,
               'artifact': 11895738765, 'archive': pin(archive), 'runner_revision': runner_source,
               'provider_source_commit': '21f87e4806c36cd1ae9e45f912a5817d6878b774',
               'observation': observation, 'identities': identities, 'owners': owners,
               'scope': 'Ten complete actual action owners and four original queued numerical controls on hosted Apple Paravirtual device with custom diagnostic MoltenVK. Three expected owner-self-test observations remain separately preserved in raw archive.',
               'raw_scope': 'Complete native ZIP split into consecutive32MiB parts, including provider source tar, build output, all owners/input seals, live loaded-image records and shader captures.',
               'physical_apple_gpu_claimed': False, 'stock_sdk_runtime_qualified': False,
               'external_runtime_domain_added': False, 'production_adopted': False, 'qualification_claimed': False}
    message = 'Preserve four original queued numerical controls on custom native MoltenVK (#48)\n'
elif kind in ('native-escape', 'native-side-effects', 'native-metadata-compile', 'native-header-boundary'):
    work = root / 'out/macos-moltenvk-integration'
    if kind == 'native-escape':
        label, runner_source, run, job, artifact = ('ninja-escape-a7a0c9d', 'a7a0c9df594da6a7944ba97c00823643a7de7075', 38557811215, 115831201718, 11890771513)
        archive_pin = {'bytes': 85241436, 'sha256': '7cd11604400b28bbe3c39e55cf1f4d86ec8f1cd64dc2aa67bcb79c0f6ca74c6d'}
    elif kind == 'native-side-effects':
        label, runner_source, run, job, artifact = ('build-side-effects-bfefab7', 'bfefab73d571b3cbeafb9f0f657924169f067239', 38559236199, 115835910513, 11891043528)
        archive_pin = {'bytes': 108476597, 'sha256': '56d46acefdbffb7e84fdcf6835eb48608d286e04d1038f99c10a47d92d5fcb6d'}
    elif kind == 'native-metadata-compile':
        label, runner_source, run, job, artifact = ('metadata-compile-a75d4b5', 'a75d4b547e4456b02b802ee2cee0f9e53974ef9d', 38560779867, 115841015706, 11893520167)
        archive_pin = {'bytes': 90481242, 'sha256': 'a86193c41a5f8fd5d5f431f31e8d4a9a3427169e25b92f2456386ef4e4e1f9cb'}
    else:
        label, runner_source, run, job, artifact = ('header-boundary-00be3ee', '00be3ee9040027bf8b0627d751cc3c97d8759e76', 38564864901, 115854858358, 11895486380)
        archive_pin = {'bytes': 97777079, 'sha256': '82e76c9d090586a909326a6bfdaf25b2ccbce0a2c84175c19a320062a15c2102'}
    bundle = root / ('attestations/native-vulkan/4bcd987/macos-provider-' + label)
    tag = 'evidence/issue-48-provider-' + label
    archive = work / ('native-' + str(run) + '.zip')
    assert pin(archive) == archive_pin
    import zipfile
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        prefix = 'out/macos-moltenvk-integration/observation/'
        observation = json.loads(zipped.read(prefix + 'integration.json'))
        assert observation['runner_source']['runner_revision'] == runner_source
        assert not observation['cases']
        owner = json.loads(zipped.read(prefix + 'provider-build/owner.json'))
        assert owner['status'] == 'terminal'
        if kind == 'native-escape':
            assert observation['provider_source_inputs_unchanged']
            assert owner['child_exit'] == -15 and owner['stop_reason'] == 'observed descendant escaped launch group'
            assert not owner['containment_complete'] and not owner['observation_complete']
        elif kind == 'native-side-effects':
            assert not observation['provider_source_inputs_unchanged']
            assert owner['child_exit'] == 120 and owner['stop_reason'] is None
            assert owner['containment_complete'] and owner['observation_complete']
            assert len(zipped.read(prefix + 'provider-build/stderr.log')) == 64 * 1024**2
        elif kind == 'native-metadata-compile':
            assert observation['provider_source_inputs_unchanged']
            assert owner['child_exit'] == 1 and owner['stop_reason'] is None
            assert owner['containment_complete'] and owner['observation_complete']
            stderr = zipped.read(prefix + 'provider-build/stderr.log')
            assert len(stderr) < 64 * 1024**2 and b"MVKDevice.mm:829:57: error: expected ')'" in stderr
        else:
            assert observation['provider_source_inputs_unchanged']
            assert owner['child_exit'] == 0 and owner['stop_reason'] is None
            assert owner['containment_complete'] and owner['observation_complete']
            assert json.loads(zipped.read(prefix + 'provider-build/inputs-before.json')) == json.loads(zipped.read(prefix + 'provider-build/inputs-after.json'))
            generated = zipped.read('bin/macos/macos-moltenvk-integration/MoltenVK/mvkGitRevDerived.h')
            template = (root / 'out/macos-moltenvk-identity-fixed/payload/sources/moltenvk-patched/Templates/cmake/mvkGitRevDerived.h.in').read_bytes()
            expected = template.replace(b'@MVK_GIT_REV@', observation['provider_content_discriminator'].encode())
            assert generated == expected + b'\n' and generated != expected
            assert observation['error'] == 'ValueError: generated legacy revision field does not contain the declared diagnostic discriminator'
            binary = zipped.read('bin/macos/macos-moltenvk-integration/MoltenVK/libMoltenVK.dylib')
            assert len(binary) == 5038480 and hashlib.sha256(binary).hexdigest() == 'cb81f199aead4872d04fedd80ffc8cb9b6c6cd26b75b058fde66ec8b29fcd26a'
        assert owner['whole_inputs_unchanged'] and not owner['cleanup_errors']
        assert not owner['remaining_observed'] and all(row['absent'] for row in owner['observed_births'])
        for phase in ('toolchain', 'stock-inventory', 'stock-uuid'):
            prior = json.loads(zipped.read(prefix + phase + '/owner.json'))
            assert prior['child_exit'] == 0 and prior['containment_complete'] and prior['observation_complete']
            assert prior['whole_inputs_unchanged'] and not prior['cleanup_errors']
    mapping = [(work / ('run-' + str(run) + '-failure.log'), 'workflow-failure.log')]
    summary = {'schema': 'sirius-rejected-native-build-v1', 'run': run,
               'job': job, 'artifact': artifact, 'archive': pin(archive),
               'runner_revision': runner_source, 'provider_source_commit': observation['provider_source_evidence_commit'],
               'failure': ('Ninja compiler ABI command creates its own process group; original strict owner rejects it.' if kind == 'native-escape' else
                           'Five child-Python bytecode additions violate frozen tree; stderr exactly64MiB ends mid-dyld-line and child exits120. Truncated output does not establish a compiler semantic failure.' if kind == 'native-side-effects' else
                           'Actual MVK_VERSION_STRING is a parenthesized expression; diagnostic adjacent-literal concatenation fails in MVKDevice.mm:829. Correct only identity formatting with bounded std::snprintf.' if kind == 'native-metadata-compile' else
                           'Provider build completed. Parent generated-header whole-byte check forgot CMake final LF; actual content discriminator is correct. No diagnostic UUID or original numerical control ran.'),
               'observation': observation, 'provider_build_owner': owner, 'qualification_claimed': False,
               'original_controls_started': False,
               'raw_scope': 'Complete native artifact ZIP split into consecutive 32MiB chunks; concatenation has the recorded archive bytes/hash. Contains complete provider source tar and all remote observations.'}
    message = 'Preserve rejected native provider build and complete raw artifact (#48)\n'
    if kind == 'native-header-boundary':
        for name in ('local-checks.json', 'header-final-lf-check.json', 'source-commit.json'):
            mapping.append((root / 'out/macos-moltenvk-identity-fixed-integration' / name, 'local-controls/' + name))
        summary['provider_build_completed'] = True
        summary['native_provider_binary'] = {'bytes': 5038480, 'sha256': 'cb81f199aead4872d04fedd80ffc8cb9b6c6cd26b75b058fde66ec8b29fcd26a'}
        message = 'Preserve completed native provider build and exact header-boundary refusal (#48)\n'
elif kind == 'p4-interrupted':
    work = root / 'out/science-p4-cpu-remainder'
    case = work / 'original-four-non-frame-cpu-controls'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-cpu-interrupted-batch'
    tag = 'evidence/issue-53-p4-cpu-interrupted-batch'
    result = json.loads((work / 'result.json').read_text())
    owner = json.loads((case / 'owner.json').read_text())
    assert result['passed'] is False and owner['status'] == 'terminal'
    assert owner['stop_reason'] == 'diagnostic_timeout' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert all(row['absent'] for row in owner['observed_births'])
    log = (case / 'stdout.log').read_text()
    name = 'CpuTraceBoundary.TruncatedPageThorneLiveProfileUsesDeclaredZeroTorqueEdge'
    assert '[       OK ] ' + name + ' (296730 ms)' in log
    assert '[  PASSED  ] 1 test.' in log and 'Passed  296.75 sec' in log
    assert not (case / 'ctest.xml').exists() and not (case / 'ctest.log').exists()
    assert (case / 'ctest.log.tmp7d2d4').is_file()
    mapping = [(work, 'files/out/science-p4-cpu-remainder'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'producer/sirius_backend_tests'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake', 'producer/backend-registrations.cmake'),
               (root / 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'producer/sirius_render_tests'),
               (root / 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake', 'producer/render-registrations.cmake')]
    summary = {'schema': 'sirius-completed-component-in-incomplete-batch-v1', 'original_batch_result': result,
               'owner': owner, 'completed_component': name, 'component_milliseconds': 296730,
               'component_verdict': 'Completed GTest and CTest pass observed in raw stdout of incomplete batch; no passing standalone owner or final XML exists.',
               'reference_scope': 'Original edge-change/cache-consumption fixture: prime at6, change to10, shared production temperature authority within2e-3, at least6 comparisons. Actual count not logged; no independent quadrature claim.',
               'second_case': 'Shakura-Sunyaev started without verdict', 'remaining_cases': 'Full Page-Thorne and near-extremal did not start',
               'batch_passed': False, 'p4_complete': False, 'qualification_claimed': False}
    message = 'Preserve completed disk-edge component within interrupted P4 batch (#53)\n'
elif kind == 'p4-parity-incomplete':
    work = root / 'out/science-p4-ordinary-fp32-parity'
    folder = work / 'parity-1'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-ordinary-fp32-setup-incomplete'
    tag = 'evidence/issue-53-p4-ordinary-fp32-setup-incomplete'
    result = json.loads((work / 'parity-1-result.json').read_text())
    prepared = json.loads((work / 'preparation.json').read_text())
    owner = json.loads((folder / 'owner.json').read_text())
    assert not result['passed'] and not owner['passed']
    assert owner['status'] == 'terminal' and owner['stop_reason'] == 'diagnostic_timeout'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
    assert json.loads((folder / 'inputs-before.json').read_text()) == json.loads((folder / 'inputs-after.json').read_text())
    assert '[ParityProbe] phase=' not in (folder / 'stdout.log').read_text()
    assert not (folder / 'gtest.xml').exists() and not (folder / 'ctest.xml').exists()
    assert all(not (work / ('parity-' + str(i))).exists() for i in range(2, 10))
    mapped = {}
    for observation in owner['provider_mappings']:
        assert observation['owned_session'] == owner['child_pid']
        for line in observation['maps']:
            fields = line.split(None, 5)
            assert len(fields) == 6 and not fields[5].endswith(' (deleted)')
            path = Path(fields[5])
            assert str(path) in prepared['required_provider_mappings']
            info = path.stat()
            assert tuple(int(piece, 16) for piece in fields[3].split(':')) == (os.major(info.st_dev), os.minor(info.st_dev))
            assert int(fields[4]) == info.st_ino
            mapped[str(path)] = {**pin(path), 'resolved_path': str(path.resolve(strict=True))}
    assert set(mapped) == set(prepared['required_provider_mappings'])
    assert all(row == prepared['pins'][name] for name, row in mapped.items())
    mapping = [(work, 'files/out/science-p4-ordinary-fp32-parity'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'producer/sirius_backend_tests'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake', 'producer/backend-registrations.cmake')]
    summary = {'schema': 'sirius-incomplete-original-p4-software-setup-v1', 'result': result, 'owner': owner,
               'classification': 'INCOMPLETE_DIAGNOSTIC_LIMIT', 'completed_initialization_logged': False,
               'actual_dispatch_or_numerical_verdict': 'unavailable', 'other_eight_actions_launched': False,
               'mapped_provider_modules': mapped,
               'interpretation': 'No completed initialization phase is logged. This boundary covers pipeline setup and the original zero-record initialization dispatch; the archive alone does not locate the exact blocking call.',
               'original_actual_dispatch_limit_ms': 1000, 'p4_complete': False, 'qualification_claimed': False}
    message = 'Preserve incomplete original ordinary-FP32 P4 software setup (#53)\n'
elif kind in ('p4-separated', 'p4-doppler'):
    import xml.etree.ElementTree as ET
    if kind == 'p4-separated':
        label = 'cpu-separated'
        work = root / 'out/science-p4-cpu-remaining-separated'
        selection = {'shakura-sunyaev-live': 'GeodesicTracerTest.LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile',
                     'page-thorne-live': 'GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile',
                     'near-extremal-live': 'GeodesicTracerRedshift.NearExtremalInnerDiskEmissionRemainsFinite'}
    else:
        label = 'doppler-raster'
        work = root / 'out/science-p4-doppler-raster'
        selection = {'original-doppler-raster-control': 'DopplerToggleTest.SuppressionCollapsesDiskAsymmetry'}
    bundle = root / ('attestations/scientific-core/4bcd987/p4-' + label + '-current-source')
    tag = 'evidence/issue-53-p4-' + label + '-current-source'
    observations, passed = {}, []
    for case, name in selection.items():
        result_path = work / (case + '-result.json') if kind == 'p4-separated' else work / 'result.json'
        result = json.loads(result_path.read_text())
        owner = json.loads((work / case / 'owner.json').read_text())
        assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert all(row['absent'] for row in owner['observed_births'])
        if result['passed']:
            assert owner['passed'] and owner['child_exit'] == 0 and owner['stop_reason'] is None
            document = ET.parse(work / case / 'ctest.xml').getroot()
            cases = document.findall('.//testcase')
            assert len(cases) == 1 and cases[0].get('name') == name
            assert cases[0].get('status', 'run') == 'run'
            assert not any(cases[0].findall(k) for k in ('failure', 'error', 'skipped'))
            for suite in [document, *document.findall('.//testsuite')]:
                assert all(int(suite.get(k, '0')) == 0 for k in ('failures', 'errors', 'skipped', 'disabled'))
            assert (work / case / 'ctest.log').is_file()
            passed.append(name)
        observations[case] = {'result': result, 'owner': owner}
    mapping = [(work, 'files/out/' + work.name),
               (root / 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'producer/sirius_render_tests'),
               (root / 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake', 'producer/render-registrations.cmake')]
    if kind == 'p4-separated':
        mapping += [(root / 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'producer/sirius_backend_tests'),
                    (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake', 'producer/backend-registrations.cmake')]
    summary = {'schema': 'sirius-original-finite-p4-controls-v1', 'runner_revision': parent,
               'product_source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc', 'observations': observations,
               'passed_cases': passed, 'all_selected_passed': len(passed) == len(selection),
               'original_assertions_changed': False, 'p4_complete': False, 'qualification_claimed': False}
    message = 'Preserve original finite current-source P4 ' + label + ' observations (#53)\n'
elif kind in ('p4-hint-clear', 'p4-dont-unroll'):
    dont = kind == 'p4-dont-unroll'
    work = root / ('out/p4-parity-dont-unroll-setup' if dont else 'out/p4-parity-unroll-hint-setup')
    label = 'p4-dont-unroll-setup' if dont else 'p4-unroll-hint-clear-setup'
    bundle = root / ('attestations/scientific-core/4bcd987/' + label)
    tag = 'evidence/issue-53-' + label
    results, owners = {}, {}
    for phase in ('validate', 'run'):
        result = json.loads((work / (phase + '-result.json')).read_text())
        owner = json.loads((work / phase / 'owner.json').read_text())
        assert not result['numerical_acceptance'] and not result['performance_or_adoption_claimed'] and not result['qualification_claimed']
        assert owner['status'] == 'terminal'
        if dont and phase == 'run':
            assert not result['passed'] and not owner['passed'] and owner['child_exit'] == -15
            assert owner['stop_reason'] == 'diagnostic_sampled_rss_guard'
            assert owner['peak_sampled_rss_kib'] > owner['limits']['rss_kib'] == 4194304
            markers = result['partial_setup_observation']['markers']
            assert markers[-1]['phase'] == 'vkQueueWaitIdle' and markers[-1]['event'] == 'begin'
            assert not any(m['phase'] == 'ReadBuffer' for m in markers)
            assert not (work / 'run/setup-observation.json').exists()
        else:
            assert result['passed'] and owner['passed'] and owner['child_exit'] == 0 and owner['stop_reason'] is None
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes'] and all(row['absent'] for row in owner['observed_births'])
        results[phase] = result; owners[phase] = owner
    proof = json.loads((work / 'unroll-bit-diff.json').read_text())
    assert proof['changed_word_count'] == 2494 and proof['all_other_words_identical'] and proof['exact_reversal_equals_stock'] and proof['dont_unroll_added'] == dont
    assert proof['candidate']['sha256'] == ('2ca5fd3a5eb1b23089dfe69a6b7520d0055d89e057625ae2c3d8b6bc2790c602' if dont else '2fe2dc5d2bae1a3b8c015278a7338a02f229a9490856286640938d942f1c8bac')
    review = json.loads((work / 'independent-review.json').read_text());assert not review['blocking_findings' if dont else 'material_blockers']
    mapping = [(work, 'files/out/' + work.name),
               (root / ('bin/linux-gcc/' + work.name), 'files/bin/linux-gcc/' + work.name)]
    summary = {'schema': 'sirius-' + label + '-diagnostic-v1', 'results': results, 'owners': owners,
               'exact_hint_bit_proof': proof, 'independent_review': review,
               'unchanged_baseline_elf_and_compile_source_evidence': 'evidence/issue-53-p4-original-setup-diagnostic at4632e77b0cca32fc4402c7b1ce613cd24101a61b',
               'all_observed_costs_including_owner_cleanup_retained': True,
               'normal_setup_readback_and_teardown_completed': not dont,
               'incomplete_reason': '4 GiB sampled RSS guard stopped submission; no readback or complete setup' if dont else None,
               'actual_driver_internal_unrolling_path_observed': False,
               'numerical_performance_adoption_scientific_or_release_qualification_claimed': False}
    message = ('Preserve exact DontUnroll setup diagnosis without a benefit claim (#53)\n' if dont else 'Preserve exact Unroll-hint removal setup diagnosis without a benefit claim (#53)\n')
elif kind == 'native-loader-metadata':
    work = root / 'out'
    bundle = root / 'attestations/native-vulkan/macos-loader-metadata-5b5fad9'
    tag = 'evidence/issue-48-macos-loader-metadata-5b5fad9'
    metadata = work / 'native-macos-loader-metadata.json'
    assert pin(metadata)['sha256'] == '9ef33586bea9eee5a04d815e56a49e739fc1d2b3bee53da8d5f7c049809e69ad'
    summary = json.loads(metadata.read_text())
    runner_source = summary['source_revision']
    assert runner_source == '5b5fad91c08c0533e29c4f5c7ae249725162e405'
    assert summary['scope_decision']['external_configured_vulkan_dependency_permitted']
    assert not summary['scope_decision']['sdk_free_self_contained_install_required_by_reviewed_contract']
    assert not summary['scope_decision']['native_installed_gpu_readiness_render_view_witness_executed_by_this_review']
    assert not summary['tracked_changes'] and not summary['jobs_launched']
    mapping = [(metadata, 'native-macos-loader-metadata.json')]
    message = 'Preserve actual macOS linkage metadata and external-dependency contract boundary (#48)\n'
elif kind == 'native-build-domains':
    work = root / 'out/native-provider-integration/build-domains'
    bundle = root / 'attestations/native-vulkan/build-domains-5b5fad9'
    tag = 'evidence/issue-48-native-build-domains-5b5fad9'
    summary = json.loads((work / 'verification.json').read_text())
    runner_source = summary['source_revision']
    assert runner_source == '5b5fad91c08c0533e29c4f5c7ae249725162e405' and summary['passed']
    assert summary['external_domains_admitted'] == 2 and not summary['new_runtime_domain_admitted'] and not summary['scientific_or_release_qualified']
    jobs = json.loads((work / 'jobs-api.json').read_text())['jobs']
    for role in ('macos', 'windows'):
        record = summary['results'][role]
        assert record['verified'] and record['domains'] == [role + '-native-build']
        assert record['claims']['test_report'] == {'cases': 9, 'failures': 0, 'errors': 0, 'skipped': 0}
        assert record['claims']['runtime_ready'] is False
        assert pin(work / (role + '.zip')) == record['archive']
        matching = [j for j in jobs if j['name'] == 'integration-' + role + '-no-render']
        assert len(matching) == 1 and matching[0]['status'] == 'completed' and matching[0]['conclusion'] == 'success'
    mapping = [(work / name, 'verification/' + name) for name in
               ('verify_downloaded.py', 'verification.json', 'run-api.json', 'jobs-api.json',
                'macos-artifact-api.json', 'windows-artifact-api.json', 'source-authority')]
    message = 'Preserve independently verified exact-source Windows and macOS build-only bundles (#48)\n'
elif kind == 'multiply64-compiler-screen':
    import importlib.util
    work = root / 'out/portable-multiply64-compiler'
    products = root / 'bin/linux-gcc/portable-multiply64-compiler'
    source_packet = root / 'out/portable-multiply64-screen'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-compiler-screen'
    tag = 'evidence/issue-44-multiply64-compiler-screen'
    spec = importlib.util.spec_from_file_location('preserved_multiply64_compiler', work / 'screen.py')
    screen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(screen)
    screen.verify_frozen()
    stock_reuse = screen.verify_stock_reuse()
    load = lambda path: json.loads(path.read_text())
    acceptance = load(work / 'compiler-acceptance.json')
    result = load(work / 'compiler-result.json')
    owner = load(work / 'compile-endpoint-normal/owner.json')
    action = load(work / 'actions.json')['compile-endpoint-normal']
    assert acceptance['passed'] and load(work / 'independent-terminal-review.json')['passed']
    assert acceptance['compiler_owner'] == screen.record(work / 'compile-endpoint-normal/owner.json')
    assert acceptance['compiler_result'] == screen.record(work / 'compiler-result.json')
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source'] == {'revision': parent, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert owner['argv'] == action['argv'] and owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    folder = work / 'compile-endpoint-normal'
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    for name, expected in load(folder / 'inputs-before.json').items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    for row in result['modules'].values():
        for expected in [row['module'], row['assembly'], *row['bodies'].values()]:
            assert screen.record(expected['resolved_path']) == expected
    assert result['stock_equals_current_canonical'] and result['original_divide_unpack_round_normalized_bodies_equal']
    assert result['complete_descriptor_and_workgroup_type_decorations_exact'] and result['execution_modes_exact']
    assert result['ordered_barrier_operands_exact'] and not result['runtime_or_numerical_admission'] and not result['benefit_established']
    mapping = [(source_packet, 'files/out/' + source_packet.name),
               (work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-multiply64-complete-endpoint-compiler-evidence-v1', 'result': result,
               'acceptance': acceptance, 'owner': owner, 'stock_reuse': stock_reuse,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'scope': 'One full normal Endpoint candidate compiled; exact stock module/source/tool and12 policy controls reused. Direct unsigned64 product through20/21 shift/jam; original Divide/Unpack/Round/interfaces/modes/barriers retained. No ordinary driver lowering, arithmetic, benefit, native-wide/global profile or production adoption evidence.',
               'qualification_claimed': False}
    message = 'Preserve distinct multiply64 source proof and complete Endpoint compiler screen (#44)\n'
elif kind == 'multiply64-ordinary-jit-screen':
    work = root / 'out/portable-multiply64-jit'
    products = root / 'bin/linux-gcc/portable-multiply64-jit'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-ordinary-jit-screen'
    tag = 'evidence/issue-44-multiply64-ordinary-jit-screen'
    load = lambda path: json.loads(path.read_text())
    acceptance = load(work / 'capture-acceptance.json')
    independent = load(work / 'independent-terminal-review.json')
    cost = load(work / 'independent-cost-review.json')
    assert cost['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN'
    assert cost['conclusion']['duplicate_product_removal_observed']
    assert not cost['conclusion']['net_cost_determined'] and not cost['conclusion']['production_adopted']
    for expected in cost['inputs'].values():
        assert pin(root / expected['path']) == {k: expected[k] for k in ('bytes', 'sha256')}
    owner = load(work / 'capture-candidate-endpoint-normal/owner.json')
    assert acceptance['passed'] and independent['passed']
    assert acceptance['complete_55_feature_queue_extension_order_and_chain_equal']
    assert acceptance['complete_durable_stock_report_reproduced']
    assert acceptance['fresh_inferior_and_called_JIT_maps_join_all_three_provider_files']
    for key, name in (('candidate_owner', 'owner.json'), ('candidate_code_review', 'actual-code-review.json'), ('candidate_feature_review', 'feature-review.json')):
        actual = work / 'capture-candidate-endpoint-normal' / name
        assert pin(actual) == {k: acceptance[key][k] for k in ('bytes', 'sha256')}
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source'] == {'revision': parent, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    folder = work / 'capture-candidate-endpoint-normal'
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    for name, expected in load(folder / 'inputs-before.json').items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    for name, expected in load(work / 'source-freeze.json')['pins'].items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-multiply64-ordinary-called-jit-screen-v1',
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'compiler_evidence': {'tag': 'evidence/issue-44-multiply64-compiler-screen', 'commit': '9d75e3bb8cb67ef0e323ed1d0800197b5f64887a'},
               'acceptance': acceptance, 'owner': owner, 'independent_raw_review': independent,
               'ordinary_lowering_cost_review': cost,
               'scope': 'One fresh candidate full normal Endpoint ordinary optimizer/called JIT capture and exact historical stock reuse. Complete features, live provider maps, code sections and relocation joins accepted. Dispatch is one-word short-input refusal canary; retained multiplication has not executed.',
               'retained_arithmetic_or_performance_benefit_or_adoption_claimed': False,
               'production_Int64_capability_admitted': False, 'qualification_claimed': False}
    message = 'Preserve distinct multiply64 ordinary called-code review and countercosts (#44)\n'
elif kind == 'sanitizer-debug-o1-focused':
    import xml.etree.ElementTree as ET
    work = root / 'out/sanitizer-debug-o1-verification'
    bundle = root / 'attestations/qualification/46a9599/sanitizer-debug-o1-focused'
    tag = 'evidence/issue-48-sanitizer-debug-o1-focused'
    load = lambda path: json.loads(path.read_text())
    result = load(work / 'focused-result.json')
    build = load(work / 'build-review.json')
    actions = load(work / 'actions.json')
    assert result['passed'] and result['source_revision'] == parent
    assert not result['full_sanitizer_or_Mandatory_estate_qualified']
    assert not result['scientific_or_external_runtime_domain_qualified'] and not result['installed_release_qualified']
    assert load(work / 'independent-configure-review.json')['passed']
    assert load(work / 'independent-terminal-review.json')['passed']
    owners = {}
    for label in ('configure', 'build-focused', 'contract-deaths', 'relocated-cpu-volume'):
        folder = work / label
        owner = load(folder / 'owner.json')
        assert result['owners'][label]['owner'] == pin(folder / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['argv'] == [arg.replace('{folder}', str(folder)) for arg in actions[label]['argv']]
        assert owner['limits'] == actions[label]['limits']
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in load(folder / 'inputs-before.json').items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
        owners[label] = owner
    for label, expected in [('contract-deaths', set(build['original_contract_cases'])),
                            ('relocated-cpu-volume', {'OperationalVolume.InstallRelocateInitialiseAndFailClosed'})]:
        xml = ET.parse(work / label / 'ctest.xml').getroot()
        cases = xml.findall('.//testcase')
        assert {case.get('name') for case in cases} == expected and len(cases) == len(expected)
        assert all(case.get('status') == 'run' for case in cases)
        assert not any(xml.findall('.//' + name) for name in ('failure', 'error', 'skipped'))
        assert all(int(xml.get(key, '0')) == 0 for key in ('failures', 'errors', 'skipped', 'disabled'))
    assert len(build['original_contract_cases']) == 5 and result['actual_volume_ctest_deadline_seconds'] == 600
    assert not (root / 'bin/linux-gcc-sanitize/operational-volume-audit').exists()
    mapping = [(work, 'files/out/' + work.name)]
    for name, expected in build['test_input_seals'].items():
        source = root / name
        assert pin(source) == expected, name
        if source.is_relative_to(root / 'bin/linux-gcc-sanitize'):
            mapping.append((source, 'producer-inputs/' + str(source.relative_to(root))))
    for name in ('CMakePresets.json', 'cmake/sirius_toolchain.cmake', 'cmake/sirius_testing.cmake',
                 'src/sirius/core/CMakeLists.txt', 'tests/CMakeLists.txt', 'tests/base/contracts_test.cpp',
                 'tests/operational/volume_initialization.cmake.in', 'tests/sanitizers/lsan-vulkan.supp'):
        mapping.append((root / name, 'producer-source/' + name))
    mapping.append((root / 'bin/linux-gcc-sanitize/.ninja_log', 'producer-inputs/bin/linux-gcc-sanitize/.ninja_log'))
    summary = {'schema': 'sirius-sanitizer-debug-o1-focused-evidence-v1', 'result': result,
               'owners': owners, 'actual_products': build['actual_products'], 'qualification_claimed': False,
               'scope': 'Fresh local WSL2 GCC14 Debug O1 contracts plus original relocated CPU-volume check. No remote O0 timing-fix causal claim, full sanitizer/Mandatory estate, decoded retained frame, scientific/external domain or installed release qualification.'}
    message = 'Preserve fresh O1 sanitizer contract and original relocated CPU volume checks (#48)\n'
elif kind == 'fp64-estimate-original-benefit':
    import xml.etree.ElementTree as ET
    work = root / 'out/portable-divide-fp64-estimate-benefit'
    products = root / 'bin/linux-gcc/portable-divide-fp64-estimate-benefit'
    trial = work / 'benefit'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp64-estimate-original-benefit'
    tag = 'evidence/issue-44-fp64-estimate-original-benefit'
    load = lambda path: json.loads(path.read_text())
    preparation = load(work / 'preparation.json')
    acceptance = load(work / 'trial-acceptance.json')
    disposition = load(work / 'disposition.json')
    comparison = load(trial / 'software-comparison.json')
    freeze = load(trial / 'retention-gate-frozen.json')
    sequence = load(trial / 'software-sequence.json')
    actions = load(trial / 'actions.json')
    source = load(work / 'source-freeze.json')
    for name, expected in source['pins'].items():
        assert pin(Path(name)) == expected, name
    assert acceptance['passed'] and acceptance['original_trial_complete']
    assert sequence['completed'] and freeze['sequence'] == sequence['labels'] == ['a1', 'b1', 'b2', 'a2']
    assert freeze['within_slot_repeats'] == 3 and freeze['all20_gates_required']
    assert (freeze['complete_call_minimum_improvement_percent'], freeze['every_inclusive_stage_minimum_improvement_percent'],
            freeze['observed_instance_api_minimum_improvement_percent']) == (5, -5, 0)
    assert len(comparison['checks']) == 20 and len(comparison['accounting_checks']) == 2
    assert sum(row['scope'] == 'complete_call_ms' for row in comparison['checks']) == 4
    assert sum(row['scope'] != 'complete_call_ms' for row in comparison['checks']) == 16
    assert comparison['readback_identity_module_and_completion_prerequisites_pass']
    assert comparison['pass'] == acceptance['performance_acceptance'] == all(row['pass'] for row in comparison['checks'] + comparison['accounting_checks'])
    assert not disposition['production_adopted'] and disposition['original_benefit_trial_completed']
    assert disposition['decision'] == ('PROCEED_TO_PRODUCTION_IMPLEMENTATION' if comparison['pass'] else 'REJECT_GUARDED_FP64_ESTIMATE_COST')
    for key in ('reader', 'comparator', 'launcher'):
        name = {'reader': 'software_contracts.py', 'comparator': 'compare_software.py', 'launcher': 'run_sequence.py'}[key]
        assert pin(trial / name)['sha256'] == freeze[key + '_sha256']
        assert (trial / name).read_bytes() == (work / 'original-apparatus' / name).read_bytes()
    assert freeze['actions_sha256'] == pin(trial / 'actions.json')['sha256']
    for key in ('baseline', 'candidate'):
        executable = preparation['executables'][key]
        assert pin(root / executable['path']) == {k: executable[k] for k in ('bytes', 'sha256')}
    owners = {}
    for label in freeze['sequence']:
        folder = trial / label
        owner = load(folder / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['provider_mappings'] and owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['argv'] == actions[label]['argv'] and owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in load(folder / 'inputs-before.json').items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
        xml = ET.parse(folder / 'gtest.xml').getroot()
        assert int(xml.get('tests', '-1')) == 1 and not any(xml.findall('.//' + name) for name in ('failure', 'error', 'skipped'))
        assert all(int(xml.get(key, '0')) == 0 for key in ('failures', 'errors', 'skipped', 'disabled'))
        owners[label] = owner
    protected = load(work / 'protected-before.json')
    assert len(protected) == acceptance['protected_payloads_unchanged'] == 114
    assert all(pin(root / name) == expected for name, expected in protected.items())
    assert load(work / 'producer-reuse.json')['passed']
    inherited = preparation['producer_prerequisites']['inherited_prerequisites']
    assert len(preparation['producer_prerequisites']['fifteen_fresh_terminal_owners']) == 15
    assert len(inherited['historical_three_producers']) == 3
    assert inherited['original_failed_wrapper_passed'] is False and inherited['original_failed_wrapper_child_exit'] == 1
    assert inherited['arithmetic']['all_original_static_and10owner_input_joins']
    assert load(work / 'root-source-admission.json')['passed']
    assert load(work / 'root-runtime-admission.json')['passed']
    assert load(work / 'independent-source-review.json')['passed']
    assert load(work / 'independent-terminal-review.json')['passed']
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-fp64-estimate-original-fixed-benefit-v1', 'acceptance': acceptance,
               'disposition': disposition, 'comparison': comparison, 'owners': owners,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'prerequisite_evidence_commit': preparation['prerequisite_evidence_commit'],
               'prerequisite_manifest_sha256': preparation['prior_manifest_sha256'],
               'timing_scope': 'Original Create/all11 numerical API windows; fixture device negotiation excluded as original. Full owned wall separately retained.',
               'provider_scope': 'Original benefit owner observes birth/session-bound LVP. Loader/LLVM byte inputs frozen; broader mapped provenance belongs to earlier integrated prerequisites.',
               'production_adopted': False, 'qualification_claimed': False,
               'prerequisite_topology': '15fresh owners plus historical3 producers and separate1 fullprofile; originalfailedprofile owner remains false/1.',
               'candidate_scope': 'Private guarded Float64 estimate, exact integer verification and original fallback; no shaderInt64. Historical diagnostic device also enables unusedInt64; public optional-provider admission/fallback remains separate.',
               'sole_trial_policy': 'No retries, pooled orders, cost subtraction, substituted samples or weakened gates.'}
    message = 'Preserve sole original guarded FP64 estimate fixed benefit trial and verdict (#44)\n'
elif kind == 'divide64-original-benefit':
    import xml.etree.ElementTree as ET
    work = root / 'out/portable-divide64-benefit'
    products = root / 'bin/linux-gcc/portable-divide64-benefit'
    trial = work / 'benefit'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide64-original-benefit'
    tag = 'evidence/issue-44-divide64-original-benefit'
    load = lambda path: json.loads(path.read_text())
    preparation = load(work / 'preparation.json')
    acceptance = load(work / 'trial-acceptance.json')
    disposition = load(work / 'disposition.json')
    comparison = load(trial / 'software-comparison.json')
    freeze = load(trial / 'retention-gate-frozen.json')
    sequence = load(trial / 'software-sequence.json')
    actions = load(trial / 'actions.json')
    source = load(work / 'source-freeze.json')
    for name, expected in source['pins'].items():
        assert pin(Path(name)) == expected, name
    assert acceptance['passed'] and acceptance['original_trial_complete']
    assert sequence['completed'] and freeze['sequence'] == sequence['labels'] == ['a1', 'b1', 'b2', 'a2']
    assert freeze['within_slot_repeats'] == 3 and freeze['all20_gates_required']
    assert (freeze['complete_call_minimum_improvement_percent'], freeze['every_inclusive_stage_minimum_improvement_percent'],
            freeze['observed_instance_api_minimum_improvement_percent']) == (5, -5, 0)
    assert len(comparison['checks']) == 20 and len(comparison['accounting_checks']) == 2
    assert sum(row['scope'] == 'complete_call_ms' for row in comparison['checks']) == 4
    assert comparison['readback_identity_module_and_completion_prerequisites_pass']
    assert comparison['pass'] == acceptance['performance_acceptance'] == all(row['pass'] for row in comparison['checks'] + comparison['accounting_checks'])
    assert not disposition['production_adopted'] and disposition['original_benefit_trial_completed']
    assert disposition['decision'] == ('PROCEED_TO_PRODUCTION_IMPLEMENTATION' if comparison['pass'] else 'REJECT_OPTIONAL64_DIVIDE_COST')
    for key in ('reader', 'comparator', 'launcher'):
        name = {'reader': 'software_contracts.py', 'comparator': 'compare_software.py', 'launcher': 'run_sequence.py'}[key]
        assert pin(trial / name)['sha256'] == freeze[key + '_sha256']
        assert (trial / name).read_bytes() == (work / 'original-apparatus' / name).read_bytes()
    assert freeze['actions_sha256'] == pin(trial / 'actions.json')['sha256']
    for key in ('baseline', 'candidate'):
        executable = preparation['executables'][key]
        assert pin(root / executable['path']) == {k: executable[k] for k in ('bytes', 'sha256')}
    owners = {}
    for label in freeze['sequence']:
        folder = trial / label
        owner = load(folder / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['provider_mappings'] and owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['argv'] == actions[label]['argv'] and owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in load(folder / 'inputs-before.json').items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
        xml = ET.parse(folder / 'gtest.xml').getroot()
        assert int(xml.get('tests', '-1')) == 1 and not any(xml.findall('.//' + name) for name in ('failure', 'error', 'skipped'))
        assert all(int(xml.get(key, '0')) == 0 for key in ('failures', 'errors', 'skipped', 'disabled'))
        owners[label] = owner
    protected = load(work / 'protected-before.json')
    assert len(protected) == acceptance['protected_payloads_unchanged'] == 114
    assert all(pin(root / name) == expected for name, expected in protected.items())
    assert load(work / 'producer-reuse.json')['passed']
    assert load(work / 'independent-source-review.json')['passed']
    assert load(work / 'independent-terminal-review.json')['passed']
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-optional64-original-fixed-benefit-v1', 'acceptance': acceptance,
               'disposition': disposition, 'comparison': comparison, 'owners': owners,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'prerequisite_evidence_commit': preparation['prerequisite_evidence_commit'],
               'prerequisite_manifest_sha256': preparation['prior_manifest_sha256'],
               'timing_scope': 'Original Create/all11 numerical API windows; fixture device negotiation excluded as original. Full owned wall separately retained.',
               'provider_scope': 'Original benefit owner observes birth/session-bound LVP. Loader/LLVM byte inputs frozen; broader mapped provenance belongs to earlier integrated prerequisites.',
               'production_adopted': False, 'qualification_claimed': False,
               'sole_trial_policy': 'No retries, pooled orders, cost subtraction, substituted samples or weakened gates.'}
    message = 'Preserve sole original optional64 divide fixed benefit trial and verdict (#44)\n'
elif kind == 'multiply64-original-benefit':
    import xml.etree.ElementTree as ET
    work = root / 'out/portable-multiply64-benefit'
    products = root / 'bin/linux-gcc/portable-multiply64-benefit'
    trial = work / 'benefit'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-original-benefit'
    tag = 'evidence/issue-44-multiply64-original-benefit'
    load = lambda path: json.loads(path.read_text())
    preparation = load(work / 'preparation.json')
    acceptance = load(work / 'trial-acceptance.json')
    disposition = load(work / 'disposition.json')
    comparison = load(trial / 'software-comparison.json')
    freeze = load(trial / 'retention-gate-frozen.json')
    sequence = load(trial / 'software-sequence.json')
    actions = load(trial / 'actions.json')
    source = load(work / 'source-freeze.json')
    for name, expected in source['pins'].items():
        assert pin(Path(name)) == expected, name
    assert acceptance['passed'] and acceptance['original_trial_complete']
    assert sequence['completed'] and freeze['sequence'] == sequence['labels'] == ['a1', 'b1', 'b2', 'a2']
    assert freeze['within_slot_repeats'] == 3 and freeze['all20_gates_required']
    assert (freeze['complete_call_minimum_improvement_percent'], freeze['every_inclusive_stage_minimum_improvement_percent'],
            freeze['observed_instance_api_minimum_improvement_percent']) == (5, -5, 0)
    assert len(comparison['checks']) == 20 and len(comparison['accounting_checks']) == 2
    assert sum(row['scope'] == 'complete_call_ms' for row in comparison['checks']) == 4
    assert comparison['readback_identity_module_and_completion_prerequisites_pass']
    assert comparison['pass'] == acceptance['performance_acceptance'] == all(row['pass'] for row in comparison['checks'] + comparison['accounting_checks'])
    assert not disposition['production_adopted'] and disposition['original_benefit_trial_completed']
    assert disposition['decision'] == ('PROCEED_TO_PRODUCTION_IMPLEMENTATION' if comparison['pass'] else 'REJECT_OPTIONAL64_MULTIPLY_COST')
    for key in ('reader', 'comparator', 'launcher'):
        name = {'reader': 'software_contracts.py', 'comparator': 'compare_software.py', 'launcher': 'run_sequence.py'}[key]
        assert pin(trial / name)['sha256'] == freeze[key + '_sha256']
        assert (trial / name).read_bytes() == (work / 'original-apparatus' / name).read_bytes()
    assert freeze['actions_sha256'] == pin(trial / 'actions.json')['sha256']
    for key in ('baseline', 'candidate'):
        executable = preparation['executables'][key]
        assert pin(root / executable['path']) == {k: executable[k] for k in ('bytes', 'sha256')}
    owners = {}
    for label in freeze['sequence']:
        folder = trial / label
        owner = load(folder / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['provider_mappings'] and owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['argv'] == actions[label]['argv'] and owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in load(folder / 'inputs-before.json').items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
        xml = ET.parse(folder / 'gtest.xml').getroot()
        assert int(xml.get('tests', '-1')) == 1 and not any(xml.findall('.//' + name) for name in ('failure', 'error', 'skipped'))
        assert all(int(xml.get(key, '0')) == 0 for key in ('failures', 'errors', 'skipped', 'disabled'))
        owners[label] = owner
    protected = load(work / 'protected-before.json')
    assert len(protected) == acceptance['protected_payloads_unchanged'] == 114
    assert all(pin(root / name) == expected for name, expected in protected.items())
    assert load(work / 'producer-reuse.json')['passed']
    assert load(work / 'root-source-admission.json')['passed']
    assert load(work / 'root-runtime-admission.json')['passed']
    assert load(work / 'independent-trial-review.json')['passed']
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-optional64-original-fixed-benefit-v1', 'acceptance': acceptance,
               'disposition': disposition, 'comparison': comparison, 'owners': owners,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'prerequisite_evidence_commit': preparation['prerequisite_evidence_commit'],
               'prerequisite_manifest_sha256': preparation['prior_manifest_sha256'],
               'timing_scope': 'Original Create/all11 numerical API windows; fixture device negotiation excluded as original. Full owned wall separately retained.',
               'provider_scope': 'Original benefit owner observes birth/session-bound LVP. Loader/LLVM byte inputs frozen; broader mapped provenance belongs to earlier integrated prerequisites.',
               'production_adopted': False, 'qualification_claimed': False,
               'sole_trial_policy': 'No retries, pooled orders, cost subtraction, substituted samples or weakened gates.'}
    message = 'Preserve sole original optional64 multiply fixed benefit trial and verdict (#44)\n'
elif kind == 'divide64-original-controls':
    import importlib.util
    work = root / 'out/portable-divide64-original-controls'
    products = root / 'bin/linux-gcc/portable-divide64-original-controls'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide64-original-controls'
    tag = 'evidence/issue-44-divide64-original-controls'
    spec = importlib.util.spec_from_file_location('preserved_divide64_original_controls', work / 'gate.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    gate.verify_static()
    packet = json.loads((work / 'packet.json').read_text())
    result = json.loads((work / 'prerequisites-result.json').read_text())
    assert result['passed'] and result['original_eight_controls_passed']
    assert not result['production_adopted'] and not result['benefit_executed']
    assert result['arithmetic'] == gate.verify_arithmetic()
    labels = [*packet['outputs'], *['control-' + str(index) for index in range(1, 9)]]
    assert len(labels) == len(set(labels)) == 19
    owners = {label: gate.accepted(label) for label in labels}
    runtime = json.loads((work / 'runtime-context.json').read_text())
    assert runtime['producer_output_seals'] == {label: gate.verify_output_seal(label) for label in packet['producer_labels']}
    for key, filename in [('profile', 'profile.json'), ('producer', 'producer-acceptance.json'),
                          ('feature', 'producer-feature-acceptance.json'), ('products', 'products-acceptance.json')]:
        assert runtime[key] == gate.record(work / filename)
    assert runtime['executable'] == gate.record(products / 'controls/original_candidate_controls')
    for index in range(1, 9):
        label = 'control-' + str(index)
        gate.join_runtime_inputs(label, owners[label])
        assert owners[label]['limits'] == packet['control_limits'][label]
    for label in ('verify-controls', 'verify-endpoint-routes'):
        gate.verify_output_seal(label)
    control = json.loads((work / 'control-acceptance.json').read_text())
    routes = json.loads((work / 'endpoint-route-acceptance.json').read_text())
    profile = json.loads((work / 'products-acceptance.json').read_text())
    assert control['passed'] and len(control['original_controls']) == 8 and routes['passed']
    assert control['timestamp_schema']['complete_readbacks'] == 32 and control['timestamp_schema']['ordered_markers'] == 51
    assert control['timestamp_schema']['actual_module_command_joins'] == 54
    assert len(routes['routes']) == 4 and routes['complete_readbacks'] == 56
    assert profile['passed'] and profile['changed_modules'] == 17 and profile['unchanged_modules'] == 18
    assert profile['main_arrays_changed'] == 17 and profile['main_arrays_unchanged'] == 23
    assert profile['owner_arrays_unchanged'] == 2 and profile['all_seven_tables_and_nonarray_ABI_constants_equal']
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-optional64-original-integrated-prerequisites-v1', 'result': result,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc', 'owners': owners,
               'generated_profile': profile, 'control_acceptance': control, 'endpoint_routes': routes,
               'source_feature_scope': 'Private explicitly admitted Int64 diagnostic candidate; exact original stock producer preserved. Actual new ELF/header/module joins and numerical execution are distinct from historical public feature-argument capture reuse.',
               'stock_reuse_scope': '39 exact original stock artifacts with source/builder/five-tool/owner/manifest/canonical-byte joins; full35-module candidate compiled fresh.',
               'performance_benefit_or_adoption_claimed': False, 'qualification_claimed': False,
               'next_dependency': 'Sole original fixed A1/B1/B2/A2 with three repeats per family and unchanged20+2 benefit gates. Production optional-capability/fallback admission and full qualification remain separate.'}
    message = 'Preserve optional64 complete retained profile and original eight integrated controls (#44)\n'
elif kind == 'fp64-estimate-original-controls':
    import importlib.util
    work = root / 'out/portable-divide-fp64-estimate-controls-continuation'
    products = root / 'bin/linux-gcc/portable-divide-fp64-estimate-controls-continuation'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp64-estimate-original-controls'
    tag = 'evidence/issue-44-fp64-estimate-original-controls'
    sys.path.insert(0,str(work))
    spec = importlib.util.spec_from_file_location('preserved_fp64_estimate_controls_continuation',work/'gate.py')
    gate = importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
    gate.verify_static()
    load = lambda path:json.loads(path.read_text())
    packet = load(work/'packet.json');result = load(work/'prerequisites-result.json')
    assert result['schema']=='sirius-fp64-estimate-controls-continuation-result-v1'
    assert result['passed'] and result['original_eight_controls_passed'] and not result['production_adopted'] and not result['benefit_executed']
    assert result['fresh_owner_count']==15 and result['fresh_producer_count']==5
    assert load(work/'root-source-admission.json')['passed'] and load(work/'root-runtime-admission.json')['passed']
    assert load(work/'independent-source-review.json')['passed']
    independent=load(work/'independent-terminal-review.json');assert independent['passed']
    prerequisites=gate.verify_prerequisites()
    assert result['inherited_prerequisites']==prerequisites and result['arithmetic']==prerequisites['arithmetic']
    assert prerequisites['original_failed_wrapper_passed'] is False and prerequisites['original_failed_wrapper_child_exit']==1
    assert len(prerequisites['historical_three_producers'])==3
    labels=[*packet['outputs'],*['control-'+str(i) for i in range(1,9)]]
    assert len(labels)==len(set(labels))==15 and len(packet['producer_labels'])==5
    owners={label:gate.accepted(label) for label in labels}
    for label,owner in owners.items():
        assert owner['source']=={'revision':parent,'status':''}
        assert owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
        for name,expected in load(work/label/'inputs-before.json').items():assert gate.record(root/name)==expected,name
    runtime=load(work/'runtime-context.json')
    assert runtime['inherited_full_profile']==prerequisites and runtime['arithmetic']==prerequisites['arithmetic']
    assert runtime['producer_output_seals']=={label:gate.verify_output_seal(label) for label in packet['producer_labels']}
    for key,filename in [('profile','profile.json'),('producer','producer-acceptance.json'),('feature','producer-feature-acceptance.json')]:assert runtime[key]==gate.record(work/filename)
    assert runtime['executable']==gate.record(products/'controls/original_candidate_controls')
    for i in range(1,9):
        label='control-'+str(i);gate.join_runtime_inputs(label,owners[label]);assert owners[label]['limits']==packet['control_limits'][label]
    for label in ('verify-controls','verify-endpoint-routes'):gate.verify_output_seal(label)
    controls=load(work/'control-acceptance.json');routes=load(work/'endpoint-route-acceptance.json')
    assert controls['passed'] and len(controls['original_controls'])==8 and routes['passed']
    assert (controls['timestamp_schema']['complete_readbacks'],controls['timestamp_schema']['ordered_markers'],controls['timestamp_schema']['actual_module_command_joins'])==(32,51,54)
    assert len(routes['routes'])==4 and routes['complete_readbacks']==56
    profile=load(Path(prerequisites['profile_result']['resolved_path']))['profile']
    assert profile['passed'] and (profile['changed_modules'],profile['unchanged_modules'])==(17,18)
    assert (profile['main_arrays_changed'],profile['main_arrays_unchanged'],profile['owner_arrays_unchanged'])==(17,23,2)
    assert profile['all_seven_tables_and_nonarray_ABI_constants_equal']
    mapping=[(work,'files/out/'+work.name),(products,'files/bin/linux-gcc/'+products.name)]
    ancillary = root/'out/original-cold-kerr-source-review-continuation'
    finding = load(ancillary/'source-findings.json')
    assert finding['source_only'] and finding['no_launch_prepared'] and finding['source_revision']==parent
    assert load(ancillary/'root-review.json')['passed'] and load(ancillary/'root-review.json')['finding']==pin(ancillary/'source-findings.json')
    for name,row in finding['source_pins'].items():assert pin(root/name)=={k:row[k] for k in ('bytes','sha256')}
    mapping.append((ancillary,'files/out/'+ancillary.name))
    summary={'schema':'sirius-fp64-estimate-original-integrated-controls-continuation-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc','result':result,'fresh15_owners':owners,
        'explicit_inherited_prerequisites':prerequisites,'independent_terminal_review':independent,
        'control_acceptance':controls,'endpoint_routes':routes,'full_profile':profile,
        'topology':'3historical successful producers +1separate accepted fullprofile +5fresh producer/8original controls/2readers; originalfailedprofile owner remains false/1.',
        'scope':'Private explicit17-module Float64-estimate profile. Exact producer/header/ELF/load/readback/runtime evidence; historical feature-enable arguments and extra unused hostInt64 remain explicit. No new shader compiler/validator/provider-canary rerun.',
        'ancillary_cold_frame_source_attribution':{'finding':pin(ancillary/'source-findings.json'),'scope':'Original CI last snapshot, no new frame or deadline repair claim.'},
        'performance_benefit_or_production_adoption_claimed':False,'qualification_claimed':False,
        'next_dependency':'Sole original fixed ABBA/3repeats with unchanged4whole+16stage+2inclusiveAPI benefit gates; public optional-provider fallback/admission and final science/release qualification remain separate.'}
    message='Preserve guarded FP64 estimate matched producer and original eight controls after profile recovery (#44)\n'
elif kind == 'fp64-estimate-profile-recovery':
    import importlib.util
    work = root / 'out/portable-divide-fp64-estimate-profile-recovery'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp64-estimate-profile-recovery'
    tag = 'evidence/issue-44-fp64-estimate-profile-recovery'
    sys.path.insert(0,str(work))
    spec = importlib.util.spec_from_file_location('preserved_estimate_profile_recovery',work/'recover.py')
    gate = importlib.util.module_from_spec(spec); spec.loader.exec_module(gate)
    gate.verify_static()
    load = lambda path:json.loads(path.read_text())
    assert load(work/'root-source-admission.json')['passed'] and load(work/'independent-source-review.json')['passed']
    independent = load(work/'independent-terminal-review.json');assert independent['passed']
    proof = gate.original_completed()
    owner = gate.terminal(work/'read-full-profile')
    for name,pin0 in load(work/'read-full-profile/inputs-before.json').items():assert gate.record(root/name)==pin0,name
    result = load(work/'profile-result.json');acceptance = load(work/'profile-acceptance.json')
    assert acceptance['passed'] and result['passed'] and result['original_failed_wrapper_passed'] is False and result['original_failed_wrapper_child_exit']==1
    assert result['signature_controls']['passed'] and len(result['signature_controls']['negative_controls'])==6
    assert all(x['rejected'] for x in result['signature_controls']['negative_controls'])
    profile = result['profile']
    assert profile['passed'] and len(profile['modules'])==35 and (profile['changed_modules'],profile['unchanged_modules'])==(17,18)
    assert (profile['main_arrays_changed'],profile['main_arrays_unchanged'],profile['owner_arrays_unchanged'])==(17,23,2)
    assert profile['all_seven_tables_and_nonarray_ABI_constants_equal']
    assert result['completed_producers']==proof
    seal = load(work/'read-full-profile/output-seal.json')
    assert acceptance['profile_result']==seal['profile_result']==gate.record(work/'profile-result.json')
    assert acceptance['reader_owner']==gate.record(work/'read-full-profile/owner.json')
    assert acceptance['output_seal']==gate.record(work/'read-full-profile/output-seal.json')
    assert len(result['new64_assembly'])==len(seal['new64_assembly'])==64
    assert result['new64_assembly']==seal['new64_assembly']=={str(path):gate.record(path) for path in (work/'analysis').rglob('*.spvasm')}
    assert len(result['reused6_assembly'])==6
    assert acceptance['original_failed_wrapper_passed'] is False and acceptance['original_failed_wrapper_child_exit']==1
    assert acceptance['original35_modules_40main_2owner_arrays_7tables_full_profile_accepted']
    mapping = [(work,'files/out/'+work.name)]
    summary = {'schema':'sirius-fp64-estimate-full-profile-recovery-v1','source_revision':parent,
        'result':result,'acceptance':acceptance,'owner':owner,'independent_terminal_review':independent,
        'original_false_reader_relabelled':False,'compiler_validator_policy_provider_rerun':False,
        'new64_disassemblies_and_reused6_dumps':True,'integrated_numerical_or_benefit_adoption_qualification_claimed':False,
        'next_dependency':'Separate matched producer and original eight controls continuation, joining historical3successful producers and this1successfulprofile owner while retaining oldfalse/1 reader.'}
    message = 'Preserve guarded FP64 estimate full35 profile reader recovery without recompilation (#44)\n'
elif kind == 'fp64-estimate-controls-profile-reader-stop':
    import importlib.util
    work = root / 'out/portable-divide-fp64-estimate-controls'
    products = root / 'bin/linux-gcc/portable-divide-fp64-estimate-controls'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp64-estimate-controls-profile-reader-stop'
    tag = 'evidence/issue-44-fp64-estimate-controls-profile-reader-stop'
    sys.path.insert(0,str(work))
    spec = importlib.util.spec_from_file_location('preserved_estimate_controls_stop',work/'gate.py')
    gate = importlib.util.module_from_spec(spec); spec.loader.exec_module(gate)
    gate.verify_static()
    load = lambda path: json.loads(path.read_text())
    assert len(load(work/'source-freeze.json')['pins']) == 2324
    assert load(work/'root-source-admission.json')['passed'] and load(work/'independent-source-review.json')['passed']
    independent = load(work/'independent-terminal-stop-review.json')
    assert independent['passed']
    accepted_labels = ['reuse-products-stock','compile-products-candidate','native-policy-selfcheck']
    owners = {label:gate.accepted(label) for label in accepted_labels}
    for label in accepted_labels: gate.verify_output_seal(label)
    failed = load(work/'verify-products/owner.json')
    assert failed['status']=='terminal' and failed['passed'] is False and failed['child_exit']==1
    assert failed['stop_reason'] is None and not failed['cleanup_errors'] and not failed['remaining_owned_processes']
    assert failed['source_unchanged'] and failed['whole_inputs_unchanged'] and all(x['absent'] for x in failed['observed_births'])
    owners['verify-products'] = failed
    for label,owner in owners.items():
        assert owner['source']=={'revision':parent,'status':''}
        assert owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
        assert (work/label/'inputs-before.json').read_bytes()==(work/label/'inputs-after.json').read_bytes()
        for name,pin0 in load(work/label/'inputs-before.json').items(): assert gate.record(root/name)==pin0,name
    assert 'camera_profile.py' in (work/'verify-products/stderr.log').read_text()
    assert not (work/'products-acceptance.json').exists() and not (work/'verify-products/output-seal.json').exists()
    assert not (products/'controls').exists() and not (work/'runtime-context.json').exists()
    assert all(not (work/('control-'+str(i))).exists() for i in range(1,9))
    assert len(list((products/'candidate').glob('*.spv')))==len(list((products/'stock').glob('*.spv')))==35
    assert len(list((products/'candidate').iterdir()))==len(list((products/'stock').iterdir()))==39
    assert len(list((products/'analysis').rglob('*.spvasm')))==6
    mapping = [(work,'files/out/'+work.name),(products,'files/bin/linux-gcc/'+products.name)]
    summary = {'schema':'sirius-fp64-estimate-controls-profile-reader-stop-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc','owners':owners,
        'root_finding':load(work/'root-profile-reader-stop.json'),'independent_terminal_stop_review':independent,
        'compiler_candidate35_complete':True,'stock39_reuse_complete':True,'native_policy_selfcheck_complete':True,
        'full_product_profile_accepted':False,'downstream_producer_or_numerical_controls_executed':False,
        'failed_owner_remains_failed':True,'performance_adoption_or_qualification_claimed':False,
        'next_dependency':'Separate exact32-bit signed-op reader recovery/full35 profile from completed products; reuse6 completed disassemblies and disassemble only64 missing. No compiler/validator/provider rerun, then separately admitted producer/original8 continuation.'}
    message = 'Preserve completed FP64-estimate retained compiler products and signed-type reader stop (#44)\n'
elif kind == 'fp64-estimate-final-retirement':
    work=root/'out/portable-divide-fp64-estimate-final-retirement'
    bundle=root/'attestations/software-vulkan/4498cc4/fp64-estimate-final-retirement'
    tag='evidence/issue-44-fp64-estimate-final-retirement'
    result=json.loads((work/'retirement.json').read_text())
    assert result['completed'] and len(result['removed'])==len(result['roots'])==10
    assert result['source_revision']==parent and result['tracked_status_after']==''
    assert all(not (root/name).exists() for name in result['removed'])
    for row in result['roots']:
        assert pin(root/row['bundle']/'manifest.json')==row['manifest']
    mapping=[(work,'files/out/'+work.name)]
    summary={'schema':'sirius-fp64-estimate-completed-private-retirement-v1','result':result,
        'preserved_fixed_decision':'REJECT_GUARDED_FP64_ESTIMATE_COST',
        'scope':'Ten completed private arithmetic/compiler/profile/control/benefit/source-review roots removed only after exact pushed archive/Git/manifest checks and recorded-birth retirement. No qualification or public adoption.'}
    message='Preserve retirement of completed guarded FP64 estimate private workspaces (#44)\n'
elif kind == 'original-cold-kerr-source-review':
    import zipfile
    work = root / 'out/minimum-step-floor-original-cold-kerr'
    bundle = root / 'attestations/software-vulkan/4498cc4/original-cold-kerr-source-review'
    tag = 'evidence/issue-44-original-cold-kerr-source-review'
    load = lambda path: json.loads(path.read_text())
    finding = load(work / 'source-findings.json')
    assert not finding['admitted'] and not finding['launch_prepared'] and finding['source_revision'] == parent
    assert pin(work / 'source-findings.json') == {'bytes':17129,'sha256':'19b87984d03d6d5dc92f7a88dcf9c3a24c679bd969202e01ae2210979b8fef10'}
    assert finding['exact_original_gate']['ctest'] == 'OperationalVulkan.RequiredRuntimeDispatch'
    assert finding['material_change_blocker']['finding'] == 'The exact-floor repair does not reach this Kerr retained renderer; this repair alone does not justify an unchanged failed-frame retry.'
    for name, expected in finding['material_change_blocker']['exact_original_CI_source_comparisons'].items():
        assert pin(root / name) == {k:expected[k] for k in ('bytes','sha256')}, name
        assert expected['original_CI_equal'] and (root / name).read_bytes() == subprocess.check_output(['git','show',finding['original_CI_revision']+':'+name]), name
    for name, expected in finding['freshness_blockers']['products'].items():
        assert pin(root / name) == {k:expected[k] for k in ('bytes','sha256')}, name
    prior = finding['original_CI_evidence']
    for role in ('whole_zip','manifest'):
        expected = prior[role]
        assert pin(Path(expected['resolved_path'])) == {k:expected[k] for k in ('bytes','sha256')}
    with zipfile.ZipFile(prior['whole_zip']['resolved_path']) as archive:
        raw = archive.read(prior['log_member']['name'])
        assert {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()} == {k:prior['log_member'][k] for k in ('bytes','sha256')}
    mapping = [(work,'files/out/'+work.name)]
    summary = {'schema':'sirius-original-cold-kerr-bounded-source-review-v1','source_revision':parent,'finding':finding,
               'root_rechecked_current_products_and23_CI_source_bytes':True,
               'actual_cold_gate_reexecuted':False,'render_build_or_resource_admitted':False,
               'qualification_claimed':False,'next_dependency':'Grounded repair in the actual retained Kerr path, then current matched render producer/resources and original300s full frame.'}
    message = 'Preserve cold Kerr retained-path source scope and render freshness gaps (#44 #53)\n'
elif kind == 'fp64-estimate-original-arithmetic':
    import importlib.util
    work = root / 'out/portable-divide-fp64-estimate-arithmetic'
    products = root / 'bin/linux-gcc/portable-divide-fp64-estimate-arithmetic'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp64-estimate-original-arithmetic'
    tag = 'evidence/issue-44-fp64-estimate-original-arithmetic'
    sys.path.insert(0, str(work))
    spec = importlib.util.spec_from_file_location('preserved_fp64_estimate_arithmetic', work / 'gate.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    gate.verify_static()
    load = lambda path: json.loads(path.read_text())
    result = load(work / 'arithmetic-result.json')
    packet = load(work / 'packet.json')
    runtime = load(work / 'runtime-context.json')
    independent = load(work / 'independent-terminal-review.json')
    assert load(work / 'root-source-admission.json')['passed']
    assert load(work / 'root-build-admission.json')['passed']
    assert load(work / 'independent-source-review.json')['passed']
    assert independent['passed'] and independent['source_revision'] == parent
    assert result['passed'] and result['original_records'] == 174268 and result['original_operation_count'] == 14
    assert result['complete_output_words_per_provider'] == 697072
    assert result['stock_result_reused'] and result['fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected']
    assert result['candidate_macro'] == 'SIRIUS_PB32_DIVIDE_FP64_ESTIMATE=1'
    assert len(load(work / 'source-freeze.json')['pins']) == 1886
    labels = [*packet['outputs'], 'execute-host', 'execute-gpu-candidate']
    owners = {label:gate.accepted(label) for label in labels}
    assert len(owners) == 10 and set(result['owners']) == set(owners)
    for label,owner in owners.items():
        assert owner['source'] == {'revision':parent,'status':''}
        assert owner['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
        assert owner['wall_seconds'] < 300 and owner['peak_sampled_rss_kib'] <= 4194304 and owner['minimum_sampled_available_kib'] >= 2097152
        assert result['owners'][label] == gate.record(work / label / 'owner.json')
        for name,expected in load(work / label / 'inputs-before.json').items():
            assert gate.record(root / name) == expected,name
    assert result['runtime_context'] == gate.record(work / 'runtime-context.json')
    assert result['build_review'] == runtime['build_review'] == gate.record(work / 'build-review.json')
    assert runtime['all_build_output_seals'] == {label:gate.verify_output_seal(label) for label in packet['outputs']}
    modules = gate.module_checks()
    assert gate.canonical_json(result['modules']) == gate.canonical_json(runtime['modules']) == gate.canonical_json(modules)
    assert result['stock_arithmetic_receipt'] == gate.verify_stock_reuse()
    assert result['compiler_and_jit_receipts'] == gate.verify_reference_acceptances()
    assert result['original_corpus'] == gate.corpus_checks()
    assert runtime['same_admitted_host_reuse'] == gate.verify_host_reuse()
    for label in ('execute-host','execute-gpu-candidate'):
        gate.join_runtime_inputs(label,owners[label])
    gate.provider_mappings(owners['execute-gpu-candidate'])
    assert (work / 'execute-gpu-candidate/actual.bin').read_bytes() == (work / 'expected.bin').read_bytes()
    assert pin(work / 'expected.bin') == {'bytes':2788288,'sha256':'8ff8c685140f5443bf989b281cdbf744b698ef30b1bc3b063464f0a6e3165414'}
    assert modules['candidate']['module']['sha256'] == '7a4290be1d6d52e8dea78dbe57f9d1465cb50c3ae4b62d2b9b7c9c1ef42930c4'
    mapping = [(work,'files/out/' + work.name),(products,'files/bin/linux-gcc/' + products.name)]
    summary = {'schema':'sirius-fp64-estimate-original-arithmetic-corpus-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc','result':result,'owners':owners,
        'independent_terminal_review':independent,'complete_output_bytes_per_gpu':pin(work / 'expected.bin'),
        'same_admitted_host_reuse':runtime['same_admitted_host_reuse'],
        'prelaunch_correction':'Preserved unexecuted checked draft and two evidence-join repairs: complete durable-JSON comparison plus actual prior owner/before/after receipt inputs. No compiler/test/provider retry.',
        'scope':'Unchanged original174268-record14-operation finite primitive corpus. Fresh candidate host checks all original record fields; fresh candidate GPU and exact reused stock complete697072-word outputs equal original expected. Native estimate accuracy/hit-rate remains unobserved, with no such correctness premise. No fresh numeric-process feature-argument capture claimed.',
        'integrated_retained_programs_stage_profiles_or_native_precision_admitted':False,
        'performance_benefit_or_production_adoption_claimed':False,'qualification_claimed':False,
        'next_dependency':'Original eight complete integrated controls, then unchanged original20+2 benefit/adoption gates.'}
    message = 'Preserve guarded FP64 estimate passing unchanged original arithmetic corpus (#44)\n'
elif kind == 'divide-fp64-estimate-ordinary-jit-screen':
    work = root / 'out/portable-divide-fp64-estimate-screen/driver'
    products = root / 'bin/linux-gcc/portable-divide-fp64-estimate-screen/driver'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide-fp64-estimate-ordinary-jit-screen'
    tag = 'evidence/issue-44-divide-fp64-estimate-ordinary-jit-screen'
    load = lambda path: json.loads(path.read_text())
    acceptance = load(work / 'capture-acceptance.json')
    independent = load(work / 'independent-terminal-review.json')
    cost = load(work / 'independent-cost-review.json')
    assert load(work / 'root-source-admission.json')['passed']
    assert load(work / 'independent-source-review.json')['passed']
    assert independent['passed'] and acceptance['passed']
    assert cost['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN' and cost['passed']
    assert cost['arithmetic_performance_adoption_or_qualification_claimed'] is False
    assert len(cost['input_pins']) == 46
    for expected in cost['input_pins']:
        assert pin(Path(expected['path'])) == {key:expected[key] for key in ('bytes','sha256')}
    assert len(cost['clones']) == 6 and all(row['complete_original_divide_block_is_skipped_on_no_failed_active_lane'] for row in cost['clones'])
    assert cost['lane_and_fallback_limits']['global_scope']['module_LLVM_udiv_sites_still_present'] == 59
    assert cost['lane_and_fallback_limits']['global_scope']['all_native_div_instructions_still_present'] == 191
    assert load(work / 'root-terminal-admission.json')['passed'] and load(work / 'root-cost-premise-admission.json')['passed']
    assert acceptance['original_failed_compiler_wrapper_preserved']
    assert all(acceptance[key] for key in ('complete_55_feature_queue_extension_order_and_chain_equal',
        'complete_durable_stock_report_reproduced', 'fresh_inferior_and_called_JIT_maps_join_all_three_provider_files'))
    folder = work / 'capture-candidate-endpoint-normal'
    owner = load(folder / 'owner.json')
    action = load(work / 'capture-actions.json')['capture-candidate-endpoint-normal']
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert owner['limits'] == action['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert owner['argv'] == [value.replace('{folder}',str(folder)) for value in action['argv']]
    assert owner['wall_seconds'] < 300 and owner['peak_sampled_rss_kib'] <= 4194304 and owner['minimum_sampled_available_kib'] >= 2097152
    for relative in ('inputs-before.json','inputs-after.json'):
        for name,expected in load(folder / relative).items():
            path = root / name
            assert dict(pin(path),resolved_path=str(path.resolve(strict=True))) == expected,name
    assert load(folder / 'inputs-before.json') == load(folder / 'inputs-after.json')
    frozen = load(work / 'source-freeze.json')['pins']
    assert len(frozen) == 1528 and len(load(folder / 'inputs-before.json')) == 1534
    for name,expected in frozen.items():
        path = root / name
        assert dict(pin(path),resolved_path=str(path.resolve(strict=True))) == expected,name
    for key,name in (('candidate_owner','owner.json'),('candidate_code_review','actual-code-review.json'),('candidate_feature_review','feature-review.json')):
        path = folder / name
        assert dict(pin(path),resolved_path=str(path.resolve(strict=True))) == acceptance[key]
    code = load(folder / 'actual-code-review.json')
    assert code['actual_loaded_module']['sha256'] == '87fc489831e49e8e1bdc7f60e3cf3f778c5a1297f21c4bc5d026968418c9fe7f'
    assert code['actual_loaded_module']['bytes'] == 153168
    expected = code['ordinary_ir_export']
    assert pin(folder / expected['path']) == {key:expected[key] for key in ('bytes','sha256')}
    old = root / 'attestations/software-vulkan/4bcd987/divide-fp64-estimate-compiler-reader-stop/files/out/portable-divide-fp64-estimate-screen/compiler'
    assert load(old / 'compile-endpoint-normal/owner.json')['passed'] is False
    assert load(old / 'compile-endpoint-normal/owner.json')['child_exit'] == 1
    assert not (old / 'compiler-result.json').exists() and not (old / 'compiler-acceptance.json').exists()
    mapping = [(work,'files/out/portable-divide-fp64-estimate-screen/driver'),
               (products,'files/bin/linux-gcc/portable-divide-fp64-estimate-screen/driver'),
               (work.parent / 'prelaunch-cache-pin-stop','files/out/portable-divide-fp64-estimate-screen/prelaunch-cache-pin-stop')]
    mapping.extend((path,'files/out/portable-divide-fp64-estimate-screen/' + path.name)
                   for path in sorted(work.parent.iterdir()) if path.is_file())
    summary = {'schema':'sirius-fp64-estimate-ordinary-called-jit-screen-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc','acceptance':acceptance,'owner':owner,
        'independent_terminal_review':independent,'ordinary_lowering_cost_review':cost,
        'prior_immutable_compiler_commit':'d547d32330a6f4cabb2a8f4f8345fa1ac3c79a93',
        'prior_immutable_reader_commit':'56076d8dff4b19c0bb40ffecb1f7432ce83bb1e6',
        'prior_immutable_property_commit':'3dccbf65373cd13e1dacca2aa57db195b4712eef',
        'original_failed_compiler_wrapper_passed':False,'original_failed_compiler_wrapper_child_exit':1,
        'static_input_count':1528,'actual_consumed_input_count':1534,'stock_recaptured':False,'host_recompiled':False,
        'retained_arithmetic_executed':False,'estimate_accuracy_or_hit_rate_established':False,
        'performance_benefit_or_production_adoption_claimed':False,'qualification_claimed':False,
        'scope':'One candidate-only original short-input ordinary driver code/callable capture, complete historical host and stock reuse, exact known-hash cache recovery and immutable failed-wrapper/new-reader/property joins. Six guarded estimate sites plus full fallback remain; static code/frame countercosts preserved. No retained arithmetic, hit-rate, benefit, cold-frame or qualification verdict.'}
    message = 'Preserve guarded FP64 estimate ordinary called-code screen and countercosts (#44)\n'
elif kind == 'divide-fp64-estimate-reader-recovery':
    work = root / 'out/portable-divide-fp64-estimate-screen/reader-recovery'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide-fp64-estimate-reader-recovery'
    tag = 'evidence/issue-44-divide-fp64-estimate-reader-recovery'
    load = lambda path: json.loads(path.read_text())
    frozen = load(work / 'source-freeze.json')['pins']
    assert len(frozen) == 126
    for name, expected in frozen.items():
        path = root / name
        assert dict(pin(path),resolved_path=str(path.resolve(strict=True))) == expected, name
    assert load(work / 'root-source-admission.json')['passed']
    assert load(work / 'independent-source-review.json')['passed']
    review = load(work / 'independent-terminal-review.json')
    assert review['passed'] and review['source_revision'] == parent
    folder = work / 'read-completed-endpoint'
    owner = load(folder / 'owner.json')
    action = load(work / 'actions.json')['read-completed-endpoint']
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert owner['argv'] == action['argv'] and owner['limits'] == action['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert owner['wall_seconds'] < 300 and owner['peak_sampled_rss_kib'] <= 4194304 and owner['minimum_sampled_available_kib'] >= 2097152
    before = load(folder / 'inputs-before.json')
    assert before == load(folder / 'inputs-after.json')
    for name, expected in before.items():
        path = root / name
        assert dict(pin(path),resolved_path=str(path.resolve(strict=True))) == expected, name
    approval = load(work / 'compiler-acceptance.json')
    result = load(work / 'compiler-result.json')
    controls = load(work / 'parser-control-result.json')
    assert approval['passed'] and result['passed'] and controls['passed']
    assert approval['original_wrapper_passed'] is False and approval['original_wrapper_child_exit'] == 1
    assert result['producer']['original_passed'] is False and result['producer']['original_child_exit'] == 1
    assert result['producer']['no_compiler_validator_or_policy_rerun'] is True
    assert approval['reader_owner'] == dict(pin(folder / 'owner.json'),resolved_path=str((folder / 'owner.json').resolve()))
    assert approval['compiler_result'] == dict(pin(work / 'compiler-result.json'),resolved_path=str((work / 'compiler-result.json').resolve()))
    assert controls['negative_count'] == 23 and len(controls['controls']) == 24
    assert sum(row['actual_acceptance'] for row in controls['controls']) == 1
    assert all(row['actual_acceptance'] == row['expected_acceptance'] for row in controls['controls'])
    shape = result['candidate_estimate_graph']
    assert shape['exact_product_remainder_inequalities'] and shape['original_complete_fallback_suffix_normalized_equal']
    assert shape['exact_remainder_sticky'] and shape['accepted_path_returns_same_round_sign_exponent']
    assert shape['native_estimate_accuracy_assumed'] is False and result['runtime_or_numerical_admission'] is False
    assert result['candidate_has_no_Int64_or_FToU'] and result['original_multiply_residual_unpack_round_normalized_bodies_equal']
    assert result['complete_descriptor_and_workgroup_type_decorations_exact'] and result['execution_modes_exact'] and result['ordered_barrier_operands_exact']
    for values in result['modules'].values():
        for expected in [values['module'],values['assembly'],*values['bodies'].values()]:
            path = Path(expected['resolved_path'])
            assert dict(pin(path),resolved_path=str(path.resolve(strict=True))) == expected
    old = root / 'attestations/software-vulkan/4bcd987/divide-fp64-estimate-compiler-reader-stop/files/out/portable-divide-fp64-estimate-screen/compiler'
    old_owner = load(old / 'compile-endpoint-normal/owner.json')
    assert old_owner['passed'] is False and old_owner['child_exit'] == 1
    assert not (old / 'compiler-result.json').exists() and not (old / 'compiler-acceptance.json').exists()
    assert not (root / 'out/portable-divide-fp64-estimate-screen/driver/runtime-context.json').exists()
    mapping = [(work,'files/out/portable-divide-fp64-estimate-screen/reader-recovery')]
    summary = {'schema':'sirius-fp64-estimate-reader-only-recovery-v1','source_revision':parent,
        'passed':True,'owner':owner,'result':result,'parser_controls':controls,'static_input_count':126,
        'prior_immutable_producer_commit':'d547d32330a6f4cabb2a8f4f8345fa1ac3c79a93',
        'original_wrapper_passed':False,'original_wrapper_child_exit':1,'compiler_retried':False,
        'driver_executed':False,'numerical_execution_started':False,'qualification_claimed':False,
        'scope':'Fresh pure-reader acceptance of exact previously compiled/validated complete module, typed local memory/closed Phi/UInt2 graph, full original fallback and profile.23 mutated text negatives refused. No compiler/validator/provider rerun, native estimate accuracy/hit-rate, original arithmetic, benefit, frame or qualification verdict.'}
    message = 'Preserve strict reader-only recovery of completed FP64 estimate module (#44)\n'
elif kind == 'fp64-estimate-float-controls-properties':
    work = root / 'out/portable-divide-fp64-estimate-property-query'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp64-estimate-float-controls-properties'
    tag = 'evidence/issue-44-fp64-estimate-float-controls-properties'
    load = lambda path: json.loads(path.read_text())
    frozen = load(work / 'source-freeze.json')['pins']
    assert len(frozen) == 53
    for name, expected in frozen.items():
        path = root / name
        assert dict(pin(path), resolved_path=str(path.resolve(strict=True))) == expected, name
    assert load(work / 'root-source-admission.json')['passed']
    assert load(work / 'independent-source-review.json')['passed']
    review = load(work / 'independent-terminal-review.json')
    assert review['passed'] and review['source_revision'] == parent
    actions = load(work / 'actions.json')
    assert set(actions) == {'query', 'read-properties'}
    owners = {}
    for name in ['query', 'read-properties']:
        folder = work / name
        owner = load(folder / 'owner.json')
        action = actions[name]
        assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['argv'] == [arg.replace('{folder}', str(folder)) for arg in action['argv']]
        assert owner['limits'] == action['limits'] == {'seconds':60,'rss_kib':4194304,'reserve_kib':2097152}
        assert owner['wall_seconds'] < 60 and owner['peak_sampled_rss_kib'] <= 4194304
        assert owner['minimum_sampled_available_kib'] >= 2097152
        before = load(folder / 'inputs-before.json')
        assert before == load(folder / 'inputs-after.json')
        for path_name, expected in before.items():
            path = root / path_name
            assert dict(pin(path), resolved_path=str(path.resolve(strict=True))) == expected, path_name
        if name == 'query':
            assert owner['provider_mappings'] and action['disable_software_shader_cache'] is False
        owners[name] = {key:owner[key] for key in ['wall_seconds','peak_sampled_rss_kib','limits']}
    result = load(work / 'read-properties/property-accepted.json')
    query = load(work / 'query/property-query.json')
    assert result['passed'] and result['source_revision'] == query['source_revision'] == parent
    assert result['qualification_claimed'] is False
    assert query['observation_complete'] and query['profile_admissible'] and query['instance_created'] and query['instance_destroyed']
    assert all(query[key] is False for key in ['logical_device_created','shader_pipeline_dispatch_or_WSI','estimate_accuracy_established','logical_device_float64_enablement_observed'])
    assert result['estimate_accuracy_established'] is False and result['logical_device_float64_enablement_observed'] is False
    assert all(result['profile_checks'].values())
    assert result['float_controls']['roundingModeIndependence'] in [0,1]
    mapping = [(work, 'files/out/' + work.name)]
    summary = {'schema':'sirius-fp64-estimate-float-controls-properties-v1','source_revision':parent,
        'passed':True,'owners':owners,'result':result,'static_input_count':53,'qualification_claimed':False,
        'scope':'Public instance/physical-device properties only. Eligible mixed RTE32/Float64 profile, with exact alive loader/LVP/LLVM map identities and raw structures; no logical-device enablement, shader/pipeline/dispatch/JIT, numerical precision/equivalence, performance, hardware-domain or release qualification.',
        'future_capture_rule':'Use exact archived owners/results/raw/inputs and selected provider/image identities; separately observe successful Float64-enabled VkDevice arguments and accepted exact candidate structure before dispatch.'}
    message = 'Preserve selected LVP mixed-width float-control property observation (#44)\n'
elif kind == 'divide-fp32-estimate-restoring-graph-review':
    import stat
    work = root / 'out/portable-divide-fp32-estimate-restoring-graph-review'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp32-estimate-restoring-graph-review'
    tag = 'evidence/issue-44-fp32-estimate-restoring-graph-review'
    load = lambda path: json.loads(path.read_text())
    def restoring_graph_record(path):
        return dict(pin(path), resolved_path=str(path.resolve(strict=True)),
                    mode=stat.S_IMODE(path.stat().st_mode),
                    symlink_target=os.readlink(path) if path.is_symlink() else None)
    def restoring_graph_matches(path, expected):
        actual = restoring_graph_record(path)
        assert {key:actual[key] for key in expected if key in actual} == expected, str(path)
    reviews = {}
    for name in ('independent-source-review.json','root-source-admission.json',
                 'independent-terminal-review.json','root-terminal-admission.json'):
        assert load(work/name)['passed'] is True, name
        reviews[name] = pin(work/name)
    assert load(work/'independent-terminal-review.json')['source_revision'] == parent
    assert reviews['independent-terminal-review.json'] == {'bytes':14851,'sha256':'6b6c89e87a5694569aa21bca0ede9c39fd8b1da964671ad409d7b998ffa85f07'}
    packet = load(work/'packet.json')
    assert packet['operator_revision'] == parent
    frozen = load(work/'source-freeze.json')['pins']
    assert len(frozen) == 699
    for name, expected in frozen.items(): restoring_graph_matches(Path(name),expected)
    folder = work/'read-graph'
    owner = load(folder/'owner.json')
    actions = load(work/'actions.json')
    assert list(actions) == ['read-graph'] and actions == load(work/'actions.unpinned.json')
    action = actions['read-graph']
    assert action['revision'] == parent and len(action['frozen_inputs']) == 702
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''} and owner['argv'] == action['argv']
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert not owner['provider_mappings']
    assert owner['limits'] == action['limits'] == packet['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    before = load(folder/'inputs-before.json')
    assert before == load(folder/'inputs-after.json') and len(before) == 703
    assert set(action['frozen_inputs']) <= {str(root/name) for name in before}
    for name, expected in before.items(): restoring_graph_matches(root/name,expected)
    context = load(work/'operator-context.json')
    assert context['revision'] == parent
    for key, name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:
        assert context[key] == restoring_graph_record(work/name)
    seal = load(folder/'output-seal.json')
    assert seal['label'] == 'read-graph' and seal['arithmetic_execution'] is False and seal['cost_or_runtime_or_adoption_admission'] is False
    for key, path in [('source_freeze',work/'source-freeze.json'),('operator_context',work/'operator-context.json'),('inputs_before',folder/'inputs-before.json')]:
        assert seal[key] == restoring_graph_record(path)
    assert seal['outputs'] == {str(work/name):restoring_graph_record(work/name) for name in ('graph-result.json','parser-control-result.json')}
    result = load(work/'graph-result.json')
    controls = load(work/'parser-control-result.json')
    acceptance = load(work/'graph-acceptance.json')
    assert result['schema'] == 'sirius-fp32-estimate-restoring-full-endpoint-graph-result-v1'
    assert acceptance['schema'] == 'sirius-fp32-estimate-restoring-full-endpoint-graph-acceptance-v1'
    assert result['passed'] and controls['passed'] and acceptance['passed'] and acceptance['structural_graph_accepted']
    assert not result['arithmetic_execution'] and not result['cost_or_hit_rate_or_adoption_established'] and not acceptance['arithmetic_or_cost_or_adoption_admitted']
    assert controls['positive_count'] == len(controls['accepted']) == 3
    assert controls['negative_count'] == len(controls['rejected']) == packet['required_graph_negative_controls'] == 75
    assert len({row['name'] for row in controls['rejected']}) == 75 and controls['all_five_restoring_positions_mutated']
    assert controls['helper_normalization_controls']['passed'] and controls['helper_normalization_controls']['negative_count'] == 2
    assert controls['mutants_are_complete_actual_module'] and controls['mutated_modules_assembled_or_validated'] is False
    for key, path in [('graph_owner',folder/'owner.json'),('graph_result',work/'graph-result.json'),
                      ('parser_controls',work/'parser-control-result.json'),('output_seal',folder/'output-seal.json')]:
        assert acceptance[key] == restoring_graph_record(path)
    chain = result['completed_compiler_chain']
    assert chain == acceptance['immutable_compiler_chain'] and chain['compiler_commit'] == packet['compiler_commit'] == '8078cfdda9afe5d7478be2c3727858746f22748a'
    assert chain['compiler102_static_and636_consumed_inputs_verified'] and chain['explicit_retired_root_relocation_preserves_historical_resolved_identity']
    assert chain['original_failed_passed'] is False and chain['original_failed_child_exit'] == 1
    assert acceptance['original_failed_passed'] is False and acceptance['original_failed_child_exit'] == 1
    archive = Path(packet['compiler_archive'])
    assert pin(archive/'manifest.json') == packet['compiler_manifest']
    members = load(archive/'manifest.json')['files']
    assert len(members) == packet['compiler_manifest_members'] == 79
    for name, expected in members.items(): assert not (archive/name).is_symlink() and pin(archive/name) == expected, name
    for key in ('compiler_manifest','compiler_owner','compiler_seal','compiler_acceptance','compiler_products','compiler_terminal_review',
                'compiler_independent_source_review','compiler_root_source_admission','compiler_root_terminal_admission','original_failed_owner'):
        restoring_graph_matches(Path(chain[key]['resolved_path']),chain[key])
    failed = load(Path(chain['original_failed_owner']['resolved_path']))
    assert failed['passed'] is False and failed['child_exit'] == 1
    assert result['products'] == chain['candidate_products']
    for value in result['products'].values(): restoring_graph_matches(Path(value['resolved_path']),value)
    assert result['products']['modules/candidate-endpoint-normal.spv']['sha256'] == packet['candidate_module_sha256'] == '93b2623d95fd98b68a60f254d3d9888cf884521f920396d7758ddfd9f922c40f'
    assert result['assembler_and_readback_complete_instructions_equal'] == packet['complete_module_instruction_count'] == 10418
    assert result['unchanged_normalized_PB32_helper_instances']['unchanged_non_Divide_helper_bodies'] == 23
    assert result['unchanged_normalized_PB32_helper_instances']['PB32_function_bodies'] == 24
    actual = result['actual_Divide']
    assert [row['shift'] for row in actual['restoring_steps']] == [4,3,2,1,0]
    assert all(actual[key] for key in ('corrective_UDiv_absent','only_original_fallback_UDiv_remains','original_complete_prefix_normalized_equal',
        'original_complete_fallback_suffix_normalized_equal','original_common_round_sign_exponent_join','exact_new_region_instruction_and_block_coverage'))
    retirement = load(work/'root-prerequisite-retirement.json')
    assert retirement['passed'] and pin(work/'root-prerequisite-retirement.json') == {'bytes':28098,'sha256':'4165f674597e567f56b01c8a54572062232a14abc0857db760365be2b76e14a2'}
    mapping = [(work,'files/out/'+work.name),
               (root/'out/preserve-finite-result-fp32-restoring-graph.patch','collector-source.patch')]
    summary = {'schema':'sirius-fp32-estimate-restoring-graph-evidence-v1','passed':True,'source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc','reviews':reviews,
        'graph_owner':restoring_graph_record(folder/'owner.json'),'result':restoring_graph_record(work/'graph-result.json'),
        'controls':restoring_graph_record(work/'parser-control-result.json'),'acceptance':restoring_graph_record(work/'graph-acceptance.json'),
        'output_seal':restoring_graph_record(folder/'output-seal.json'),
        'static_rich_inputs':699,'declared_inputs':702,'actual_consumed_inputs':703,'positive_count':3,'graph_negative_count':75,'helper_negative_count':2,
        'restoring_shifts':[4,3,2,1,0],'complete_module_instructions':10418,'unchanged_non_Divide_helpers':23,
        'compiler_commit':packet['compiler_commit'],'compiler_manifest':packet['compiler_manifest'],
        'original_failed_compiler_owner_passed':False,'original_failed_compiler_child_exit':1,
        'previous_distinct_UDiv_JIT_cost':'COST_PREMISE_SURVIVES_UNPROVEN, no arithmetic/benefit trial or performance failure fabricated.',
        'root_compiler_retirement':pin(work/'root-prerequisite-retirement.json'),
        'structural_restoring_graph_accepted':True,'native_or_driver_accepted':False,'arithmetic_or_runtime_or_cost_or_benefit_established':False,
        'production_adopted':False,'qualification_claimed':False,
        'scope':'One pure full Endpoint structural reader/control owner plus separate terminal acceptance. Actual five-step exact correction, guard/closedCFG/prefix/fullfallback/commonRound and23helper receipts bind immutable successful compiler5products. Source mutation negatives are not independently assembled/validated. No compiler/policy/tool/provider repeated; originalfalse1 stays separate. Native/cost and original numerical/integrated/benefit gates remain future.'}
    message = 'Preserve FP32 estimate five-step restoring full Endpoint graph and source controls without runtime admission (#44)\n'
elif kind == 'divide-fp32-estimate-restoring-compiler':
    import stat
    work = root / 'out/portable-divide-fp32-estimate-restoring-screen'
    products = root / 'bin/linux-gcc/portable-divide-fp32-estimate-restoring-screen'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp32-estimate-restoring-compiler'
    tag = 'evidence/issue-44-fp32-estimate-restoring-compiler'
    load = lambda path: json.loads(path.read_text())
    def restoring_record(path):
        return dict(pin(path), resolved_path=str(path.resolve(strict=True)),
                    mode=stat.S_IMODE(path.stat().st_mode),
                    symlink_target=os.readlink(path) if path.is_symlink() else None)
    def restoring_matches(path, expected):
        actual = restoring_record(path)
        assert {key:actual[key] for key in expected if key in actual} == expected, str(path)
    reviews = {}
    for name in ('independent-source-review.json','root-source-admission.json',
                 'independent-terminal-review.json','root-terminal-admission.json'):
        assert load(work/name)['passed'] is True, name
        reviews[name] = pin(work/name)
    assert load(work/'independent-terminal-review.json')['source_revision'] == parent
    assert reviews['independent-terminal-review.json'] == {'bytes':15873,'sha256':'691e29c9e9b86b45c14e5d071d539baf76b41537eeeded504bd61f90a3c5725c'}
    frozen = load(work/'source-freeze.json')
    assert (len(frozen['private_inputs']),len(frozen['external_inputs'])) == (44,58)
    for name, expected in frozen['private_inputs'].items(): restoring_matches(work/name,expected)
    for name, expected in frozen['external_inputs'].items(): restoring_matches(Path(name),expected)
    folder = work/'compile-endpoint-normal'
    owner = load(folder/'owner.json')
    actions = load(work/'actions.json')
    assert list(actions) == ['compile-endpoint-normal']
    action = actions['compile-endpoint-normal']
    assert action['revision'] == parent and len(action['frozen_inputs']) == 105
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert not owner['provider_mappings'] and owner['argv'] == action['argv']
    assert owner['limits'] == action['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    before = load(folder/'inputs-before.json')
    assert before == load(folder/'inputs-after.json') and len(before) == 636
    assert set(action['frozen_inputs']) <= {str(root/name) for name in before}
    for name, expected in before.items(): restoring_matches(root/name,expected)
    context = load(work/'operator-context.json')
    assert context['revision'] == parent and context['status'] == ''
    for key, name in [('actions','actions.json'),('preparation','preparation.json'),('source_freeze','source-freeze.json')]:
        assert context[key] == restoring_record(work/name)
    original_actions = load(work/'actions.unpinned.json')
    assert original_actions['compile-endpoint-normal']['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT'
    original_actions['compile-endpoint-normal']['revision'] = parent
    assert original_actions == actions
    result = load(work/'compiler-products.json')
    acceptance = load(work/'compiler-producer-acceptance.json')
    seal = load(folder/'output-seal.json')
    assert result['schema'] == 'sirius-fp32-estimate-compiler-products-v1' and result['passed'] is True
    assert acceptance['schema'] == 'sirius-fp32-estimate-compiler-producer-acceptance-v1' and acceptance['passed'] is True
    for value in (result,acceptance,seal): assert value['structural_graph_accepted'] is False
    assert result['runtime_or_numerical_admission'] is False and acceptance['runtime_or_numerical_admission'] is False
    assert result['benefit_established'] is False and result['original_compile_ssa_assemble_validate_complete']
    assert all(result[key] for key in ('stock_equals_current_canonical','complete_descriptor_and_workgroup_type_decorations_equal',
        'execution_modes_equal','ordered_barrier_operands_equal','exact_one_NoContraction_and_original_two_RTE_lines_inverse_equal_raw',
        'private_single_FDiv32_policy_and_original_normal_policy_passed'))
    for key, path in [('compiler_owner',folder/'owner.json'),('compiler_products',work/'compiler-products.json'),('output_seal',folder/'output-seal.json')]:
        assert acceptance[key] == restoring_record(path)
    assert seal['label'] == 'compile-endpoint-normal'
    for key, path in [('source_freeze',work/'source-freeze.json'),('operator_context',work/'operator-context.json'),
                      ('inputs_before',folder/'inputs-before.json'),('recipe',work/'compiler-recipes.json')]:
        assert seal[key] == restoring_record(path)
    expected_outputs = [products/'compiler/modules'/('candidate-endpoint-normal'+suffix)
                        for suffix in ('.compiler.spv','.raw.spvasm','.spvasm','.spv','.final.spvasm')]
    assert seal['outputs'] == {str(path):restoring_record(path) for path in [*expected_outputs,work/'compiler-products.json']}
    assert {path for path in products.rglob('*') if path.is_file()} == set(expected_outputs)
    for key in ('raw_post_SSA_module','raw_preannotation_assembly','actual_assembler_input'):
        restoring_matches(Path(result[key]['resolved_path']),result[key])
    for module in result['modules'].values():
        for key in ('module','assembly'): restoring_matches(Path(module[key]['resolved_path']),module[key])
    assert result['stock_reuse']['source_joins'] == 8 and result['stock_reuse']['tool_joins'] == 16
    for key in ('archive_manifest','original_owner','original_acceptance','guard_controls'):
        value = result['stock_reuse'][key]
        restoring_matches(Path(value['resolved_path']),value)
    previous = result['previous_actual_cost_premise']
    assert previous['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN'
    for key in ('manifest','cost_result','root_admission'):
        value = previous[key]
        restoring_matches(Path(value['resolved_path']),value)
    prior_cost = load(Path(previous['cost_result']['resolved_path']))
    prior_admission = load(Path(previous['root_admission']['resolved_path']))
    assert prior_cost['passed'] and prior_admission['passed']
    assert prior_cost['decision'] == prior_admission['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN'
    assert prior_cost['retained_arithmetic_executed_by_canary'] is False and prior_cost['production_adopted'] is False
    assert prior_admission['retained_arithmetic_or_performance_or_hitrate_or_adoption_or_qualification_admitted'] is False
    origins = load(work/'origins.json')
    failed_archive = Path(origins['previous_header_and_frozen_builder']['archive'])
    assert pin(failed_archive/'manifest.json')['sha256'] == origins['previous_header_and_frozen_builder']['manifest_sha256']
    failed_member = 'files/out/portable-divide-fp32-estimate-screen/compile-endpoint-normal/owner.json'
    assert pin(failed_archive/failed_member) == load(failed_archive/'manifest.json')['files'][failed_member]
    failed_owner = load(failed_archive/failed_member)
    assert failed_owner['passed'] is False and failed_owner['child_exit'] == 1
    operational = {}
    for name in ('root-prerequisite-retirement.json','root-prior-jit-issue-receipt.json'):
        receipt = load(work/name)
        assert receipt['passed'] is True and receipt['archive_commit'] == previous['commit']
        operational[name] = pin(work/name)
    mapping = [(work,'files/out/'+work.name),(products,'files/bin/linux-gcc/'+products.name),
               (root/'out/preserve-finite-result-fp32-restoring-compiler.patch','collector-source.patch')]
    summary = {'schema':'sirius-fp32-estimate-restoring-compiler-evidence-v1','passed':True,
        'source_revision':parent,'product_revision':context['product_revision'],'reviews':reviews,
        'compiler_owner':restoring_record(folder/'owner.json'),'producer_acceptance':restoring_record(work/'compiler-producer-acceptance.json'),
        'compiler_products':restoring_record(work/'compiler-products.json'),'output_seal':restoring_record(folder/'output-seal.json'),
        'static_inputs_verified':102,'private_inputs':44,'external_inputs':58,'declared_inputs':105,'actual_consumed_inputs':636,
        'five_shader_products':{str(path.relative_to(products)):restoring_record(path) for path in expected_outputs},
        'stock_scope':'Existing accepted complete stock module/program/effective recipe8source+16tool/24archive identities and12originalpolicy controls reused; no fresh stock compiler/policy rerun or all591-current metadata claim. No duplicate stock raw capture copied.',
        'prior_distinct_UDiv_JIT':{'commit':previous['commit'],'manifest':previous['manifest'],'cost_result':previous['cost_result'],
            'root_admission':previous['root_admission'],'decision':'COST_PREMISE_SURVIVES_UNPROVEN','numerical_or_benefit_trial_run':False},
        'original_failed_compiler_owner_passed':False,'original_failed_compiler_child_exit':1,
        'operational_root_receipts':operational,
        'structural_restoring_graph_accepted':False,'native_or_driver_accepted':False,
        'arithmetic_or_runtime_or_cost_or_benefit_established':False,'production_adopted':False,'qualification_claimed':False,
        'scope':'One successful original full normal Endpoint compiler/SSA/private exact annotation/original scoped policy/assembly/validation/readback producer. All5actualshader products+result and original input/context/recipe seal preserved. Source preparation/draft flags remain historical; actual terminal owner and producer acceptance are authoritative. Exact emitted restoring graph, native path/cost, original arithmetic/integrated controls/20+2benefit and production/fullqualification remain future.'}
    message = 'Preserve successful FP32 estimate five-step restoring compiler producer without graph or numerical admission (#44)\n'
elif kind == 'divide-fp32-estimate-ordinary-jit-screen':
    import stat
    work = root / 'out/portable-divide-fp32-estimate-jit'
    products = root / 'bin/linux-gcc/portable-divide-fp32-estimate-jit'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp32-estimate-ordinary-jit-screen'
    tag = 'evidence/issue-44-fp32-estimate-ordinary-jit-screen'
    load = lambda path: json.loads(path.read_text())
    def jit_record(path):
        return dict(pin(path), resolved_path=str(path.resolve(strict=True)),
                    mode=stat.S_IMODE(path.stat().st_mode),
                    symlink_target=os.readlink(path) if path.is_symlink() else None)
    def jit_matches(path, expected):
        actual = jit_record(path)
        assert {key:actual[key] for key in expected if key in actual} == expected, str(path)
    inputs = load(work/'inputs.json')
    acceptance = load(work/'capture-acceptance.json')
    independent = load(work/'independent-terminal-review.json')
    assert inputs['operator_revision'] == independent['source_revision'] == parent
    assert acceptance['passed'] is True and independent['passed'] is True
    reviews = {}
    for name in ('independent-source-review.json','root-source-admission.json',
                 'independent-terminal-review.json','root-terminal-admission.json',
                 'root-cost-premise-admission.json'):
        assert load(work/name)['passed'] is True, name
        reviews[name] = pin(work/name)
    assert reviews['independent-terminal-review.json'] == {'bytes':20099,'sha256':'b20dfa1e1ff2c7f6234f5bb368e6c258995ca4d0a8986ad650228bafc39f23b9'}
    assert reviews['root-terminal-admission.json'] == {'bytes':1933,'sha256':'bdef09dfd62c0d768c57da9eeff77aa3792dc963b2dcc5d1a264e812124a13b3'}
    assert reviews['root-cost-premise-admission.json'] == {'bytes':2351,'sha256':'641a304f3f3d0b8e7c40d00cbf24b4a967eb0eebb061a8ed3f49740bde429c40'}
    folder = work/'capture-candidate-endpoint-normal'
    owner = load(folder/'owner.json')
    actions = load(work/'capture-actions.json')
    assert list(actions) == ['capture-candidate-endpoint-normal']
    action = actions['capture-candidate-endpoint-normal']
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''} and action['revision'] == parent
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert len(owner['observed_births']) == 4 and all(row['absent'] for row in owner['observed_births'])
    assert owner['argv'] == [value.replace('{folder}',str(folder)) for value in action['argv']]
    assert owner['limits'] == action['limits'] == inputs['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert owner['wall_seconds'] < 300 and owner['peak_sampled_rss_kib'] <= 4194304 and owner['minimum_sampled_available_kib'] >= 2097152
    before = load(folder/'inputs-before.json')
    assert before == load(folder/'inputs-after.json') and len(before) == 939
    assert {str(root/name) for name in before} == set(action['frozen_inputs'])
    for name, expected in before.items(): jit_matches(root/name,expected)
    frozen = load(work/'source-freeze.json')['pins']
    assert len(frozen) == 933
    for name, expected in frozen.items(): jit_matches(Path(name),expected)
    context = load(work/'runtime-context.json')
    assert context['revision'] == parent
    for key, name in [('source_freeze','source-freeze.json'),('preparation','preparation.json'),
                      ('capture_actions','capture-actions.json'),('reuse_acceptance','reuse-acceptance.json')]:
        assert context[key] == jit_record(work/name)
    assert load(work/'capture-actions.unpinned.json')['capture-candidate-endpoint-normal']['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT'
    assert json.dumps(load(work/'capture-actions.unpinned.json')).replace('PENDING_OPERATOR_CLEAN_CONTEXT',parent) == json.dumps(actions)
    assert list(products.iterdir()) == [products/'short_input_jit_probe']
    assert context['host'] == jit_record(products/'short_input_jit_probe')
    assert pin(products/'short_input_jit_probe') == inputs['canary_identity'] == {'bytes':253224,'sha256':'c2b9093df67b85608ec3ca530ec910fe206bede457ac4224059ebc36c66364af'}
    assert context['modules']['candidate'] == jit_record(work/'module-inputs/candidate-endpoint-normal.spv')
    assert pin(work/'module-inputs/candidate-endpoint-normal.spv') == inputs['candidate_module'] == {'bytes':152844,'sha256':'0535c0f6a689d2d59896543882e1f5300fdd67fe70509b7dda0996f59ffd03d9'}
    for key, name in [('candidate_owner','owner.json'),('candidate_code_review','actual-code-review.json'),('candidate_feature_review','feature-review.json')]:
        assert acceptance[key] == jit_record(folder/name)
    code = load(folder/'actual-code-review.json')
    assert independent['actual_module_to_called_code']['one_completed_dispatch'] == code['single_completed_dispatch']
    assert code['single_completed_dispatch']['captured_ns'] == 1792029865229152350
    assert code['actual_loaded_module'] == load(folder/'actual-loaded-module.json')
    assert {key:code['actual_loaded_module'][key] for key in ('bytes','sha256')} == inputs['candidate_module']
    assert code['ordinary_pipelines'] == inputs['ordinary_pipelines']
    assert pin(folder/code['ordinary_ir_export']['path']) == {key:code['ordinary_ir_export'][key] for key in ('bytes','sha256')}
    assert len(code['linked_section_review']['live_sections']) == 6
    assert all(row['all_defined_relocations_and_uncovered_bytes_exact'] for row in code['linked_section_review']['live_sections'])
    bindings = load(folder/'bindings.json')
    events = [row for row in bindings['jit_events'] if row.get('file') == code['actual_entry']['file']]
    assert len(events) == 1
    event = events[0]
    assert pin(folder/event['file']) == {'bytes':event['object_bytes'],'sha256':event['sha256']}
    for row in event['live_sections_after_pointer_lookup']:
        assert pin(folder/row['path']) == {key:row[key] for key in ('bytes','sha256')}
    features = load(folder/'feature-review.json')
    assert features == load(folder/'feature-verification.json') and len(features['features']) == 55
    assert features['features']['shaderFloat64'] == 1 and features['features']['shaderInt64'] == 0
    assert all(value == 0 for name,value in features['features'].items() if name != 'shaderFloat64')
    assert features['successful_logical_device_used_for_actual_complete_shader']
    assert features['queue_family'] == 0 and features['enabled_extensions'] == [] and features['feature_chain_values'] == []
    assert acceptance['fresh_complete55_feature_queue_extension_order_and_chain_recorded']
    assert acceptance['fresh_inferior_and_called_JIT_maps_join_all_three_provider_files']
    assert acceptance['canonical_stock_raw_feature_arguments_available'] is False
    reuse = load(work/'reuse-acceptance.json')
    assert reuse['passed'] and reuse['stock'] == acceptance['stock_receipt'] and reuse['graph'] == acceptance['graph_receipt']
    chain = acceptance['graph_receipt']
    assert chain['original_failed_passed'] is False and chain['original_failed_child_exit'] == 1
    assert chain['graph_passed'] and chain['policy_passed'] and chain['arithmetic_execution'] is False
    for key in ('original_failed_owner','policy_owner','policy_result','policy_acceptance',
                'graph_owner','graph_result','graph_acceptance','graph_output_seal','graph_terminal_review'):
        jit_matches(Path(chain[key]['resolved_path']),chain[key])
    assert load(Path(chain['original_failed_owner']['resolved_path']))['passed'] is False
    assert load(Path(chain['original_failed_owner']['resolved_path']))['child_exit'] == 1
    cost_path = work/'cost-review/cost-premise-review.json'
    assert pin(cost_path) == {'bytes':76283,'sha256':'d1ada19a47e5f8ac4859f6a90718af23733e3e70d84de39aa36bbcfcc1b0d6ca'}
    cost = load(cost_path)
    assert cost['schema'] == 'sirius-fp32-estimate-ordinary-actual-cost-premise-review-v1'
    assert cost['passed'] and cost['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN'
    assert cost['no_arithmetic_or_benefit_admission'] and not cost['production_adopted']
    assert not cost['performance_or_coldness_or_hit_rate_claim'] and not cost['retained_arithmetic_executed_by_canary']
    assert len(cost['input_pins']) == 59 and len(cost['all_six_Divide_clones']) == 6
    for name, expected in cost['input_pins'].items(): jit_matches(Path(name),expected)
    assert cost['actual_called_object_join']['actual_entry'] == code['actual_entry']
    assert cost['actual_called_object_join']['called_object'] == jit_record(folder/event['file'])
    assert cost['actual_called_object_join']['ordinary_pipelines'] == code['ordinary_pipelines']
    stock = load(work/'stock-reuse-inputs.json')
    archive = Path(stock['archive'])
    manifest_path = archive/stock['manifest_name']
    assert pin(manifest_path)['sha256'] == stock['manifest_sha256']
    manifest = load(manifest_path)
    assert manifest['source_revision'] == inputs['product_revision']
    assert manifest['evidence_commit'] == stock['evidence_commit'] and manifest['evidence_tag'] == stock['evidence_tag']
    prefix = str(archive/'files') + '/'
    selected = {name[len(prefix):] for name in frozen if name.startswith(prefix) and name[len(prefix):] in manifest['files']}
    history = Path(stock['historical_work'])
    raw_prefix = str((history/'capture-stock-endpoint-normal').relative_to(archive/'files')) + '/'
    assert len(selected) == acceptance['stock_receipt']['selected_manifest_members'] == 56
    assert len([name for name in selected if name.startswith(raw_prefix)]) == 33
    assert {str(path.relative_to(archive/'files')) for path in (history/'capture-stock-endpoint-normal').iterdir() if path.is_file()} == {name for name in selected if name.startswith(raw_prefix)}
    lookup = {(row['original_key'],row['pin']['sha256']):row for row in stock['joins']}
    assert len(lookup) == len(stock['joins']) == acceptance['stock_receipt']['matched_historical_input_versions'] == 638
    for phase,count in [('compile-probe',545),('link-probe',547),('capture-stock-endpoint-normal',633)]:
        old = load(history/phase/'owner.json')
        assert old['passed'] and old['child_exit'] == 0 and old['whole_inputs_unchanged'] and old['source_unchanged']
        assert not old['cleanup_errors'] and not old['remaining_owned_processes']
        old_before = load(history/phase/'inputs-before.json')
        assert old_before == load(history/phase/'inputs-after.json') and len(old_before) == count
        for name, expected in old_before.items():
            row = lookup[(name,expected['sha256'])]
            assert row['pin'] == expected and phase in row['phases']
            assert pin(Path(row['actual'])) == {key:expected[key] for key in ('bytes','sha256')}
    required = {str((history/name).relative_to(archive/'files')) for name in ('input-recovery.json','preparation.json','short_input_jit_probe.cpp')}
    for phase in stock['phases']:
        required.update(str((history/phase/name).relative_to(archive/'files')) for name in ('owner.json','inputs-before.json','inputs-after.json'))
    required.update(row['manifest_member'] for row in stock['joins'] if row['method'] in ('whole_byte_manifest_archived_owned_file','exact_original_input_recovery'))
    assert required <= selected
    mapping = [(work,'files/out/'+work.name),(products/'short_input_jit_probe','files/bin/linux-gcc/'+products.name+'/short_input_jit_probe')]
    selected_records = []
    for name in sorted(selected):
        path = archive/'files'/name
        assert not path.is_symlink() and pin(path) == manifest['files'][name]
        relative = 'stock-comparison/files/' + name
        mapping.append((path,relative))
        selected_records.append({'archive_member':name,'source':str(path),'preserved_relative':relative,**jit_record(path)})
    mapping.append((manifest_path,'stock-comparison/'+stock['manifest_name']))
    canonical_records = []
    old_link = load(history/'link-probe/inputs-before.json')
    for name in ('bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o','bin/linux-gcc/src/sirius/libsirius_base.a'):
        row = lookup[(name,old_link[name]['sha256'])]
        path = Path(row['actual'])
        assert pin(path) == {key:old_link[name][key] for key in ('bytes','sha256')}
        relative = 'stock-comparison/canonical-producer/' + Path(name).name
        if not any(source == path for source,_ in mapping): mapping.append((path,relative))
        else: relative = next(destination for source,destination in mapping if source == path)
        canonical_records.append({'original_input':name,'source':str(path),'preserved_relative':relative,**jit_record(path)})
    assert canonical_records[0]['sha256'] == acceptance['stock_receipt']['canonical_device_object']['sha256']
    mapping.append((root/'out/preserve-finite-result-fp32-jit-screen.patch','collector-source.patch'))
    summary = {'schema':'sirius-fp32-estimate-ordinary-called-jit-screen-evidence-v1','passed':True,
        'source_revision':parent,'product_revision':inputs['product_revision'],
        'capture_acceptance':pin(work/'capture-acceptance.json'),'owner':pin(folder/'owner.json'),'reviews':reviews,
        'static_inputs_verified':933,'whole_owned_inputs_verified':939,'stock_distinct_input_versions':638,
        'candidate_host':context['host'],'candidate_module':context['modules']['candidate'],
        'actual_code_review':pin(folder/'actual-code-review.json'),'feature_review':pin(folder/'feature-review.json'),
        'ordinary_cost_review':{'receipt':pin(cost_path),'schema':cost['schema'],'decision':cost['decision'],
            'six_source_IR_native_loci':6,'input_pins':59,'premise_survives':True,'benefit_hypothesis_proven':False},
        'selected_stock_closure':{'original_archive':str(archive),'original_commit':stock['evidence_commit'],
            'original_tag':stock['evidence_tag'],'original_manifest':jit_record(manifest_path),
            'manifested_files':selected_records,'canonical_producer_files':canonical_records,
            'manifested_member_count':56,'stock_raw_member_count':33,
            'scope':'Only these selected stock raw comparison, compile/link/capture ownership, producer and recovery bytes enter this new Git evidence tree. Modes record the actual preserved raw files; generic Git100644 is not original executable-mode provenance. Remote availability requires root push/tag verification. No whole older multicapture archive or current toolchain-library backup claim.'},
        'original_failed_compiler_owner_passed':False,'original_failed_child_exit':1,
        'separate_policy_passed':True,'separate_graph_passed':True,'stock_recaptured':False,'host_recompiled':False,
        'timestamp_correction':independent['receipt_correction'],
        'retained_arithmetic_executed':False,'benefit_trial_executed':False,'benefit_trial_failed':False,
        'performance_or_hitrate_or_production_adoption_claimed':False,'qualification_claimed':False,
        'scope':'One candidate-only original short-input ordinary compiler/callable capture. Actual full32-only shader, fresh55features/Int64zero/unusedFloat64, original pipelines/export/lookup/called entry, all6allocatedsections/definedrelocs/uncoveredbytes and three-provider birth/session/dev/inode maps accepted. Undefined callbacks consistency only; old stock raw55 capture absent. Cost premise survives unproven, with countercosts preserved; next action is root source decision, not current arithmetic or benefit admission.'}
    message = 'Preserve FP32 estimate ordinary called-code capture, surviving unproven cost premise and selected recoverable stock closure (#44)\n'
elif kind == 'divide-fp32-estimate-graph-review':
    import stat
    work = root / 'out/portable-divide-fp32-estimate-graph-review'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp32-estimate-graph-review'
    tag = 'evidence/issue-44-fp32-estimate-graph-review'
    load = lambda path: json.loads(path.read_text())
    def graph_record(path):
        return dict(pin(path), mode=stat.S_IMODE(path.stat().st_mode),
                    resolved_path=str(path.resolve(strict=True)),
                    symlink_target=os.readlink(path) if path.is_symlink() else None)
    packet = load(work/'packet.json')
    folder = work/'read-graph'
    owner = load(folder/'owner.json')
    actions = load(work/'actions.json')
    assert list(actions) == ['read-graph'] and (work/'actions.json').read_bytes() == (work/'actions.unpinned.json').read_bytes()
    action = actions['read-graph']
    assert packet['operator_revision'] == action['revision'] == parent
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert not owner['provider_mappings'] and owner['argv'] == action['argv']
    assert owner['limits'] == action['limits'] == packet['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert (folder/'inputs-before.json').read_bytes() == (folder/'inputs-after.json').read_bytes()
    before = load(folder/'inputs-before.json')
    for name, expected in before.items():
        path = root/name
        assert str(path.resolve(strict=True)) == expected['resolved_path'], name
        assert pin(path) == {key:expected[key] for key in ('bytes','sha256')}, name
    frozen = load(work/'source-freeze.json')['pins']
    assert len(frozen) == 708
    for name, expected in frozen.items(): assert graph_record(Path(name)) == expected, name
    context = load(work/'operator-context.json')
    assert context['revision'] == parent
    for key, relative in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:
        assert context[key] == graph_record(work/relative)
    reviews = {}
    for name in ('independent-emitted-source-graph-review.json','independent-source-review.json',
                 'root-source-admission.json','independent-terminal-review.json','root-terminal-admission.json'):
        review = load(work/name)
        assert review['passed'] is True, name
        reviews[name] = pin(work/name)
    assert load(work/'independent-terminal-review.json')['source_revision'] == parent
    retirement = load(work/'root-prerequisite-retirement.json')
    assert retirement['passed'] and retirement['graph_inputs_must_use_immutable_archive']
    assert retirement['archive_commit'] == packet['recovery_commit']
    assert retirement['archive_manifest'] == packet['recovery_manifest']
    assert retirement['original_failed_compiler_owner_still_false_exit1'] is True
    assert retirement['separate_finisher_passed'] is True and retirement['finisher_child_exit'] == 0
    assert all(row['removed'] for row in retirement['roots'])
    for archive, expected in [(Path(packet['recovery_archive']),packet['recovery_manifest']),
                              (Path(packet['stop_archive']),{'sha256':packet['stop_manifest_sha256']})]:
        assert all(pin(archive/'manifest.json')[key] == value for key,value in expected.items())
        for name, expected_file in load(archive/'manifest.json')['files'].items():
            assert not (archive/name).is_symlink() and pin(archive/name) == expected_file, name
    result = load(work/'graph-result.json')
    controls = load(work/'parser-control-result.json')
    acceptance = load(work/'graph-acceptance.json')
    seal = load(folder/'output-seal.json')
    assert result['schema'] == 'sirius-fp32-estimate-full-endpoint-graph-result-v1' and result['passed'] is True
    assert controls['passed'] is True and controls['positive_count'] == len(controls['accepted']) == 3
    assert controls['negative_count'] == len(controls['rejected']) == 38
    assert len({row['name'] for row in controls['rejected']}) == 38
    assert controls['mutants_are_complete_actual_module'] and controls['mutated_modules_assembled_or_validated'] is False
    helper_controls = controls['helper_normalization_controls']
    assert helper_controls['passed'] and helper_controls['negative_count'] == len(helper_controls['rejected']) == 2
    helpers = result['unchanged_normalized_PB32_helper_instances']
    assert helpers['PB32_function_bodies'] == 24 and helpers['unique_PB32_names'] == 23
    assert helpers['unchanged_non_Divide_helper_bodies'] == len(helpers['all_present_helpers']) == 23
    assert [(row['name'],row['occurrence']) for row in helpers['all_present_helpers'] if row['name']=='PB32WideShiftJam'] == [('PB32WideShiftJam',0),('PB32WideShiftJam',1)]
    assert result['named_PB32WideProduct_presence'] == {'stock':False,'candidate':False}
    assert result['assembler_and_readback_complete_instructions_equal'] == 10339
    assert result['completed_original_policy_receipt_reused_without_rerun'] is True
    assert result['complete_interface_modes_and_ordered_barriers_equal_stock'] is True
    assert result['single_NoContraction_FDiv_scope_retained'] and result['original_normal_sum_policy_retained']
    assert result['arithmetic_execution'] is False and result['cost_or_hit_rate_or_adoption_established'] is False
    assert result['native_accuracy_or_RTE_assumption'] is False
    graph = result['actual_Divide']
    assert graph['original_complete_prefix_normalized_equal'] and graph['original_complete_fallback_suffix_normalized_equal']
    assert graph['original_common_round_sign_exponent_join'] and graph['exact_new_region_instruction_and_block_coverage']
    assert acceptance['schema'] == 'sirius-fp32-estimate-full-endpoint-graph-acceptance-v1'
    assert acceptance['passed'] is True and acceptance['structural_graph_accepted'] is True
    assert acceptance['arithmetic_or_cost_or_adoption_admitted'] is False
    assert acceptance['original_failed_passed'] is False and acceptance['original_failed_child_exit'] == 1
    for key, relative in [('graph_owner','read-graph/owner.json'),('output_seal','read-graph/output-seal.json'),
                          ('graph_result','graph-result.json'),('parser_controls','parser-control-result.json')]:
        assert acceptance[key] == graph_record(work/relative)
    assert seal['label'] == 'read-graph'
    for key, relative in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before','read-graph/inputs-before.json')]:
        assert seal[key] == graph_record(work/relative)
    assert seal['outputs'] == {str(work/name):graph_record(work/name) for name in ('graph-result.json','parser-control-result.json')}
    assert seal['arithmetic_execution'] is False and seal['cost_or_runtime_or_adoption_admission'] is False
    chain = acceptance['immutable_policy_chain']
    assert result['products'] == chain['candidate_products']
    assert result['stock_module'] == chain['stock_module'] and result['stock_assembly'] == chain['stock_assembly']
    assert result['candidate_header_sha256'] == packet['candidate_header_sha256'] and result['candidate_macro'] == packet['candidate_macro']
    assert {key:value for key,value in chain.items() if key!='fresh_saved_birth_absence'} == {key:value for key,value in result['completed_policy_chain'].items() if key!='fresh_saved_birth_absence'}
    assert chain['original_failed_passed'] is False and chain['original_failed_child_exit'] == 1
    assert chain['old96_static_and630_consumed_inputs_verified'] and chain['finisher662_static_and666_consumed_inputs_verified']
    for value in chain.values():
        if isinstance(value,dict) and 'resolved_path' in value:
            assert value == graph_record(Path(value['resolved_path']))
    for value in chain['candidate_products'].values(): assert value == graph_record(Path(value['resolved_path']))
    original_owner = load(Path(chain['old_failed_owner']['resolved_path']))
    policy_owner = load(Path(chain['policy_owner']['resolved_path']))
    assert original_owner['status']=='terminal' and original_owner['passed'] is False and original_owner['child_exit']==1
    assert policy_owner['status']=='terminal' and policy_owner['passed'] is True and policy_owner['child_exit']==0
    assert all(row['saved_birth_absent'] for rows in chain['fresh_saved_birth_absence'].values() for row in rows)
    findings = load(work/'source-findings.json')
    assert findings['held_unexecuted_prior_freeze_sha256'] == '543ba6a906bacff3ce259ff83c16b14966deb00538a1ad867154faad61d2f8fd'
    assert findings['operator_or_reader_executed'] is False
    assert pin(work/'source-freeze.json')['sha256'] != findings['held_unexecuted_prior_freeze_sha256']
    mapping = [(work,'files/out/portable-divide-fp32-estimate-graph-review')]
    summary = {'schema':'sirius-fp32-estimate-full-endpoint-graph-evidence-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
        'owner':owner,'graph_result':result,'controls':controls,'acceptance':acceptance,'reviews':reviews,
        'static_inputs_verified':708,'whole_owned_inputs_verified':len(before),
        'root_prerequisite_retirement':pin(work/'root-prerequisite-retirement.json'),
        'original_failed_owner_passed':False,'original_failed_child_exit':1,'separate_policy_owner_passed':True,
        'structural_graph_accepted':True,'arithmetic_or_runtime_or_cost_or_adoption_admitted':False,
        'qualification_claimed':False,
        'scope':'One pure complete normal Endpoint graph reader accepted actual quotient/correction/guard/closed-CFG/original-prefix/fallback shape, all23nonDivide helper bodies and source mutation controls. Immutable compiler false/1 and separate successful policy finisher stay distinct. No compiler/policy/tools/provider execution repeated; original arithmetic, full retained profile/integrated controls and benefit remain future.'}
    message = 'Preserve complete FP32 estimate Endpoint graph and strict source controls without arithmetic admission (#44)\n'
elif kind == 'divide-fp32-estimate-policy-recovery':
    import stat
    import re
    work = root / 'out/portable-divide-fp32-estimate-policy-recovery'
    products = root / 'bin/linux-gcc/portable-divide-fp32-estimate-policy-recovery'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp32-estimate-policy-recovery'
    tag = 'evidence/issue-44-fp32-estimate-policy-recovery'
    load = lambda path: json.loads(path.read_text())
    def recovery_record(path):
        return dict(pin(path), mode=stat.S_IMODE(path.stat().st_mode),
                    resolved_path=str(path.resolve(strict=True)),
                    symlink_target=os.readlink(path) if path.is_symlink() else None)
    packet = load(work / 'packet.json')
    folder = work / 'finish-policy'
    owner = load(folder / 'owner.json')
    action = load(work / 'actions.json')['finish-policy']
    assert packet['operator_revision'] == action['revision'] == parent
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['source'] == {'revision':parent,'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert not owner['provider_mappings'] and owner['argv'] == action['argv']
    assert owner['limits'] == action['limits'] == packet['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert (folder/'inputs-before.json').read_bytes() == (folder/'inputs-after.json').read_bytes()
    before = load(folder/'inputs-before.json')
    for name, expected in before.items():
        path = root/name
        assert str(path.resolve(strict=True)) == expected['resolved_path'], name
        assert pin(path) == {key:expected[key] for key in ('bytes','sha256')}, name
    frozen = load(work/'source-freeze.json')['pins']
    assert len(frozen) == 662
    for name, expected in frozen.items():
        assert recovery_record(Path(name)) == expected, name
    context = load(work/'operator-context.json')
    assert context['revision'] == parent
    for key, relative in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:
        assert context[key] == recovery_record(work/relative)
    reviews = {}
    for name in ('independent-source-review.json','root-source-admission.json',
                 'independent-terminal-review.json','root-terminal-admission.json'):
        review = load(work/name)
        assert review['passed'] is True, name
        reviews[name] = pin(work/name)
    assert load(work/'independent-terminal-review.json')['source_revision'] == parent
    retirement = load(work/'root-prerequisite-retirement.json')
    assert retirement['passed'] and retirement['recovery_inputs_must_use_immutable_archive']
    assert retirement['archive_commit'] == packet['stop_commit']
    assert retirement['archive_manifest']['sha256'] == packet['stop_manifest_sha256']
    assert retirement['original_owner_passed'] is False and retirement['original_child_exit'] == 1
    assert all(row['removed'] for row in retirement['roots'])
    archive = Path(packet['stop_archive'])
    assert pin(archive/'manifest.json')['sha256'] == packet['stop_manifest_sha256']
    for name, expected in load(archive/'manifest.json')['files'].items():
        assert not (archive/name).is_symlink() and pin(archive/name) == expected, name
    old = archive/'files'/packet['original_out']
    def archived_input(name, expected):
        original = Path(name)
        if not original.is_absolute(): original = root/original
        owned = any(original.is_relative_to(root/packet[key]) for key in ('original_out','original_bin'))
        actual = archive/'files'/original.relative_to(root) if owned else original
        assert pin(actual) == {key:expected[key] for key in ('bytes','sha256')}, name
        assert expected['resolved_path'] == str(original if owned else actual.resolve(strict=True)), name
        for key in ('mode','symlink_target'):
            if key in expected: assert recovery_record(actual)[key] == expected[key], name
    old_frozen = load(old/'source-freeze.json')
    assert len(old_frozen['private_inputs']) + len(old_frozen['external_inputs']) == 96
    for name, expected in old_frozen['private_inputs'].items():
        archived_input(root/packet['original_out']/name, expected)
    for name, expected in old_frozen['external_inputs'].items(): archived_input(name, expected)
    old_folder = old/'compile-endpoint-normal'
    old_owner = load(old_folder/'owner.json')
    assert old_owner['status'] == 'terminal' and old_owner['passed'] is False and old_owner['child_exit'] == 1
    assert old_owner['source'] == {'revision':parent,'status':''}
    assert old_owner['source_unchanged'] and old_owner['whole_inputs_unchanged']
    assert old_owner['stop_reason'] is None and not old_owner['cleanup_errors'] and not old_owner['remaining_owned_processes']
    assert old_owner['observed_births'] and all(row['absent'] for row in old_owner['observed_births'])
    assert (old_folder/'inputs-before.json').read_bytes() == (old_folder/'inputs-after.json').read_bytes()
    old_before = load(old_folder/'inputs-before.json')
    assert len(old_before) == 630
    for name, expected in old_before.items(): archived_input(name, expected)
    result = load(work/'policy-result.json')
    acceptance = load(work/'policy-acceptance.json')
    seal = load(folder/'output-seal.json')
    assert result['schema'] == 'sirius-fp32-estimate-policy-recovery-result-v1' and result['passed'] is True
    assert acceptance['schema'] == 'sirius-fp32-estimate-policy-recovery-acceptance-v1' and acceptance['passed'] is True
    assert acceptance['remaining_assembly_validation_complete'] is True
    assert result['old_failed_passed'] is False and result['old_failed_child_exit'] == 1
    assert acceptance['original_failed_passed'] is False and acceptance['original_failed_child_exit'] == 1
    proof = result['completed_stop']
    assert proof['old_failed_passed'] is False and proof['old_failed_child_exit'] == 1
    assert proof['old_failed_owner'] == acceptance['original_failed_owner'] == recovery_record(old_folder/'owner.json')
    for value in proof.values():
        if isinstance(value,dict) and 'resolved_path' in value:
            assert value == recovery_record(Path(value['resolved_path']))
    for document in (result,acceptance,seal): assert document['structural_graph_accepted'] is False
    assert result['runtime_or_numerical_admission'] is False and acceptance['runtime_or_numerical_admission'] is False
    assert result['benefit_established'] is False and result['compiler_or_SSA_or_old_disassembly_repeated'] is False
    assert result['one_annotation_and_original_two_RTE_lines_reversal_equal_raw'] is True
    assert result['unchanged_frozen_private_and_original_normal_policies_passed'] is True
    assert result['complete_type_interface_execution_modes_and_ordered_barriers_equal_stock'] is True
    assert result['recipes'] == load(work/'recipes.json')['only_new_tool_calls']
    for key, relative in [('reader_owner','finish-policy/owner.json'),('output_seal','finish-policy/output-seal.json'),('policy_result','policy-result.json')]:
        assert acceptance[key] == recovery_record(work/relative)
    assert seal['label'] == 'finish-policy'
    for key, relative in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),
                          ('inputs_before','finish-policy/inputs-before.json'),('result','policy-result.json')]:
        assert seal[key] == recovery_record(work/relative)
    output_paths = [products/name for name in ('candidate-endpoint-normal.spvasm','candidate-endpoint-normal.spv','candidate-endpoint-normal.final.spvasm')]
    assert {path for path in products.rglob('*') if path.is_file()} == set(output_paths)
    assert result['outputs'] == seal['outputs'] == {str(path):recovery_record(path) for path in output_paths}
    raw = Path(proof['raw_assembly']['resolved_path']).read_text()
    assembly = output_paths[0].read_text()
    entry = re.search(r'OpEntryPoint GLCompute (%\S+)',raw)[1]
    cap = 'OpCapability Shader\n               OpCapability RoundingModeRTE'
    mode = '\n               OpExecutionMode '+entry+' RoundingModeRTE 32'
    annotation = '\n               OpDecorate %proposal NoContraction'
    assert assembly.count(cap) == assembly.count(mode) == assembly.count(annotation) == 1
    assert assembly.replace(cap,'OpCapability Shader',1).replace(mode,'',1).replace(annotation,'',1) == raw
    mapping = [(work,'files/out/portable-divide-fp32-estimate-policy-recovery'),
               (products,'files/bin/linux-gcc/portable-divide-fp32-estimate-policy-recovery')]
    summary = {'schema':'sirius-fp32-estimate-policy-recovery-evidence-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
        'owner':owner,'result':result,'acceptance':acceptance,'reviews':reviews,
        'static_inputs_verified':662,'whole_owned_inputs_verified':len(before),
        'original_static_inputs_verified':96,'original_owned_inputs_verified':630,
        'root_prerequisite_retirement':pin(work/'root-prerequisite-retirement.json'),
        'original_failed_owner_passed':False,'original_failed_child_exit':1,
        'remaining_policy_assembly_validation_readback_complete':True,
        'structural_graph_accepted':False,'runtime_or_numerical_admission':False,
        'benefit_established':False,'qualification_claimed':False,
        'scope':'Separate exact annotation and unchanged policy/assembly/validation/readback finisher reuses immutable completed Slang/SSA/raw-disassembly outputs. Original false/child1 owner remains failed and retired. Quotient/correction/dominance/fallback graph, provider/JIT, original arithmetic, integrated controls and benefit remain future; no proposal accuracy, hit-rate, adoption or qualification claim.'}
    message = 'Preserve exact FP32 estimate policy and assembly recovery without graph admission (#44)\n'
elif kind == 'divide-fp32-estimate-compiler-policy-stop':
    import stat
    import re
    work = root / 'out/portable-divide-fp32-estimate-screen'
    products = root / 'bin/linux-gcc/portable-divide-fp32-estimate-screen/compiler'
    bundle = root / 'attestations/software-vulkan/4bcd987/fp32-estimate-compiler-policy-stop'
    tag = 'evidence/issue-44-fp32-estimate-compiler-policy-stop'
    load = lambda path: json.loads(path.read_text())
    def compiler_record(path):
        return dict(pin(path), mode=stat.S_IMODE(path.stat().st_mode),
                    resolved_path=str(path.resolve(strict=True)),
                    symlink_target=os.readlink(path) if path.is_symlink() else None)
    folder = work / 'compile-endpoint-normal'
    owner = load(folder / 'owner.json')
    actions = load(work / 'actions.json')
    assert list(actions) == ['compile-endpoint-normal']
    action = actions['compile-endpoint-normal']
    context = load(work / 'operator-context.json')
    assert context['revision'] == action['revision'] == parent and context['status'] == ''
    assert context['product_revision'] == '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
    assert owner['status'] == 'terminal' and owner['passed'] is False and owner['child_exit'] == 1
    assert owner['source'] == {'revision': parent, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert not owner['provider_mappings']
    assert owner['argv'] == action['argv']
    assert owner['limits'] == action['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    before = load(folder / 'inputs-before.json')
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    for name, expected in before.items():
        path = root / name
        assert str(path.resolve(strict=True)) == expected['resolved_path'], name
        assert pin(path) == {key:expected[key] for key in ('bytes','sha256')}, name
    for key, relative in [('actions','actions.json'), ('source_freeze','source-freeze.json'),
                          ('preparation','preparation.json')]:
        assert context[key] == compiler_record(work / relative)
    frozen = load(work / 'source-freeze.json')
    counts = {}
    for prefix, field in [(work,'private_inputs'), (root,'external_inputs')]:
        for name, expected in frozen[field].items():
            assert compiler_record(prefix / name) == expected, name
        counts[field] = len(frozen[field])
    preparation = load(work / 'preparation.json')
    assert preparation['private_input_count'] == counts['private_inputs']
    assert preparation['external_input_count'] == counts['external_inputs']
    reviews = {}
    for name in ('independent-arithmetic-source-review.json',
                 'independent-compiler-producer-source-review.json', 'root-source-admission.json',
                 'independent-compiler-policy-stop-review.json', 'root-terminal-admission.json'):
        review = load(work / name)
        assert review['passed'] is True, name
        reviews[name] = pin(work / name)
    stop = load(work / 'root-terminal-stop.json')
    assert stop['schema'] == 'sirius-fp32-estimate-root-compiler-policy-stop-v1'
    assert stop['stop_audit_passed'] is True and stop['passed'] is False
    assert stop['original_false1_owner_retained'] is True
    assert stop['unchanged_compiler_retry_permitted'] is False
    assert stop['no_graph_numerical_provider_cost_adoption_qualification_claim'] is True
    assert stop['all96_static_pins_match'] is True and sum(counts.values()) == 96
    assert stop['all_consumed_input_rows_current'] == len(before)
    for key, relative in [('owner','compile-endpoint-normal/owner.json'),
                          ('before','compile-endpoint-normal/inputs-before.json'),
                          ('after','compile-endpoint-normal/inputs-after.json'),
                          ('source_freeze','source-freeze.json')]:
        assert stop[key] == compiler_record(work / relative)
    assert stop['completed_operations'] == ['Slang compile','original SSA rewrite','optimized module disassembly']
    assert set(stop['not_completed']) == {'candidate policy acceptance','assembly','Vulkan1.2 validation',
                                         'final disassembly','producer output seal'}
    raw_module = products / 'modules/candidate-endpoint-normal.compiler.spv'
    raw_disassembly = products / 'modules/candidate-endpoint-normal.spvasm'
    assert {path for path in products.rglob('*') if path.is_file()} == {raw_module,raw_disassembly}
    assert stop['partial_products'] == {str(path):compiler_record(path) for path in (raw_module,raw_disassembly)}
    assert raw_module.read_bytes()[:4] == bytes.fromhex('03022307') and raw_module.stat().st_size % 4 == 0
    for path in (work/'compiler-products.json', work/'compiler-producer-acceptance.json',
                 folder/'output-seal.json', products/'modules/candidate-endpoint-normal.spv',
                 products/'modules/candidate-endpoint-normal.final.spvasm'):
        assert not path.exists(), path
    stderr = (folder / 'stderr.log').read_text()
    assert 'diagnostic_divide_fp32_normal_sum_controls' in stderr
    assert stderr.rstrip().endswith('ValueError: diagnostic FP32 proposal lost its single scoped division/order')
    text = raw_disassembly.read_text()
    fdiv = stop['actual_fdiv']
    assert len(fdiv) == 4 and tuple(fdiv) in re.findall(r'^\s*(%\S+) = OpFDiv (%\S+) (%\S+) (%\S+)\s*$',text,re.M)
    assert not re.search(r'\bOpDecorate '+re.escape(fdiv[0])+r' NoContraction\b',text)
    unexecuted = root / 'out/preserve-finite-result-fp32-unexecuted-producer.patch'
    assert unexecuted.is_file()
    mapping = [(work, 'files/out/portable-divide-fp32-estimate-screen'),
               (products, 'files/bin/linux-gcc/portable-divide-fp32-estimate-screen/compiler'),
               (unexecuted, 'prelaunch-collector/unexecuted-producer.patch')]
    summary = {'schema':'sirius-fp32-estimate-compiler-policy-stop-evidence-v1',
        'source_revision':parent, 'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
        'owner':owner, 'root_terminal_stop':stop, 'reviews':reviews,
        'static_inputs_verified':sum(counts.values()), 'whole_owned_inputs_verified':len(before),
        'raw_optimized_module':pin(raw_module), 'raw_disassembly':pin(raw_disassembly),
        'completed_compiler_producer':False, 'assembly_or_validation_complete':False,
        'producer_output_seal_created':False, 'structural_graph_accepted':False,
        'runtime_or_numerical_admission':False, 'benefit_established':False, 'qualification_claimed':False,
        'unexecuted_success_collector':{'patch':pin(unexecuted),'executed':False,
            'base':'attestations/software-vulkan/4498cc4/fp64-estimate-final-retirement/preserve-finite-result.py',
            'draft_sha256':'59f3f31e8c322875ef4599db58722764097dd0de10aeca7eb469440cbf980fb4'},
        'scope':'Original false/child1 owner and exact stable source/inputs preserved with successful Slang/SSA/disassembly outputs. Private diagnostic policy refused missing NoContraction on the untrusted FDiv32. Assembly, validation, completed producer seal and all graph/provider/numerical/cost/adoption/qualification remain unaccepted. No unchanged compiler retry is admitted.'}
    message = 'Preserve FP32 estimate compiler policy stop and exact partial outputs (#44)\n'
elif kind == 'divide-fp64-estimate-compiler-reader-stop':
    import stat
    work = root / 'out/portable-divide-fp64-estimate-screen/compiler'
    products = root / 'bin/linux-gcc/portable-divide-fp64-estimate-screen/compiler'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide-fp64-estimate-compiler-reader-stop'
    tag = 'evidence/issue-44-divide-fp64-estimate-compiler-reader-stop'
    load = lambda path: json.loads(path.read_text())
    folder = work / 'compile-endpoint-normal'
    owner = load(folder / 'owner.json')
    action = load(work / 'actions.json')['compile-endpoint-normal']
    assert owner['status'] == 'terminal' and owner['passed'] is False and owner['child_exit'] == 1
    assert owner['source'] == {'revision': parent, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert owner['argv'] == action['argv'] and owner['limits'] == action['limits'] == {'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    before = load(folder / 'inputs-before.json')
    assert before == load(folder / 'inputs-after.json')
    for name, expected in before.items():
        path = root / name
        assert str(path.resolve(strict=True)) == expected['resolved_path'], name
        assert pin(path) == {key:expected[key] for key in ('bytes','sha256')}, name
    frozen = load(work / 'source-freeze.json')
    count = 0
    for prefix, entries in [(work, frozen['private_inputs']), (root, frozen['external_inputs'])]:
        for name, expected in entries.items():
            path = prefix / name
            actual = dict(pin(path), mode=stat.S_IMODE(path.stat().st_mode), resolved_path=str(path.resolve(strict=True)), symlink_target=os.readlink(path) if path.is_symlink() else None)
            assert actual == expected, name
            count += 1
    assert count == 83
    assert load(work / 'independent-source-review.json')['passed']
    independent = load(work / 'independent-terminal-stop-review.json')
    assert independent['passed'] and independent['source_revision'] == parent
    assert not (work / 'compiler-result.json').exists() and not (work / 'compiler-acceptance.json').exists()
    assert not (root / 'out/portable-divide-fp64-estimate-screen/driver/runtime-context.json').exists()
    stderr = (folder / 'stderr.log').read_text()
    assert 'shape = estimate_shape' in stderr and 'estimate_graph.py' in stderr and 'line 70' in stderr
    assert stderr.endswith('AssertionError: []\n')
    candidate = products / 'modules/candidate-endpoint-normal.spv'
    assert candidate.read_bytes()[:4] == bytes.fromhex('03022307')
    assert (products / 'modules/candidate-endpoint-normal.final.spvasm').is_file()
    mapping = [(work, 'files/out/portable-divide-fp64-estimate-screen/compiler'),
               (products, 'files/bin/linux-gcc/portable-divide-fp64-estimate-screen/compiler')]
    for name in ['contracts.json','origins.json','commands.json','root-compiler-admission.json']:
        mapping.append((work.parent / name, 'files/out/portable-divide-fp64-estimate-screen/' + name))
    summary = {'schema':'sirius-divide-fp64-estimate-compiler-reader-stop-v1','source_revision':parent,
        'product_revision':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
        'owner':owner,'compiler_screen_accepted':False,'driver_executed':False,'numerical_execution_started':False,
        'qualification_claimed':False,'candidate_module':pin(candidate),'static_inputs_verified':83,
        'failure':'Complete Slang/SSA/assembly/validation returned; strict SSA reader rejected its expected guess-expression shape (no matching ID).',
        'scope':'Source and owner stability passed; structural compiler acceptance did not. No native estimate accuracy, executed retained arithmetic, benefit, cold-frame or full qualification claim.',
        'next_action':'Preserve original producer/reader and full emitted outputs. If exact compiled structure supports it, prepare a separate strict representation reader and reuse these outputs; do not recompile an unchanged candidate or proceed to device execution before structural acceptance.'}
    message = 'Preserve exact FP64 estimate compiler structural-reader stop (#44)\n'
elif kind == 'minimum-step-floor-focused-controls':
    import xml.etree.ElementTree as ET
    work = root / 'out/minimum-step-floor-current'
    controls = work / 'controls'
    bundle = root / 'attestations/scientific-core/4498cc4/minimum-step-floor-focused-controls'
    tag = 'evidence/issue-53-minimum-step-floor-focused-controls'
    load = lambda path: json.loads(path.read_text())
    result = load(work / 'focused-result.json')
    assert result['passed'] and result['source_revision'] == parent
    assert result['completed_case_count'] == result['required_case_count'] == 6
    assert result['qualification_claimed'] is False
    commands = load(controls / 'commands.json')
    actions = load(controls / 'actions.json')
    assert commands['order'] == ['source-governance', 'module-capabilities', 'fp32', 'fp32comp', 'fp64', 'ellis-physical-ray']
    assert set(actions) == set(commands['order'])
    assert [row['case'] for row in result['cases']] == commands['order']
    frozen = load(controls / 'source-freeze.json')
    assert frozen['revision'] == parent and len(frozen['files']) == 1033
    for name, expected in frozen['files'].items():
        path = root / name
        assert str(path.resolve(strict=True)) == expected['resolved_path'], name
        assert pin(path) == {key: expected[key] for key in ('bytes', 'sha256')}, name
    assert load(controls / 'root-runtime-admission.json')['passed']
    assert load(work / 'root-source-admission.json')['passed']
    independent = load(work / 'independent-terminal-review.json')
    assert independent['passed'] and independent['source_revision'] == parent
    build_actions = load(work / 'actions.json')
    owners = {}
    for label, folder, action in [(name, work / name, build_actions[name]) for name in ['configure', 'build-backend']] + [(name, controls / name, actions[name]) for name in commands['order']]:
        owner = load(folder / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        expected_argv = [arg.replace('{folder}', str(folder)) for arg in action['argv']]
        assert owner['argv'] == expected_argv and owner['limits'] == action['limits']
        assert owner['wall_seconds'] < action['limits']['seconds']
        assert owner['peak_sampled_rss_kib'] <= action['limits']['rss_kib']
        assert owner['minimum_sampled_available_kib'] >= action['limits']['reserve_kib']
        before = load(folder / 'inputs-before.json')
        assert before == load(folder / 'inputs-after.json')
        for name, expected in before.items():
            path = root / name
            assert str(path.resolve(strict=True)) == expected['resolved_path'], name
            assert pin(path) == {key: expected[key] for key in ('bytes', 'sha256')}, name
        if label in actions:
            assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
            assert action['disable_software_shader_cache'] is False
            assert action['software_provider'] == (label in {'fp32', 'fp32comp', 'fp64', 'ellis-physical-ray'})
            if action['software_provider']:
                assert owner['provider_mappings']
            xml = ET.parse(folder / ('ctest.xml' if label == 'source-governance' else 'gtest.xml')).getroot()
            cases = xml.findall('.//testcase')
            assert len(cases) == 1 and int(xml.get('tests', '-1')) == 1
            assert not any(xml.findall('.//' + key) for key in ('failure', 'error', 'skipped'))
            assert not any(int(node.get(key, '0')) for node in [xml, *xml.findall('.//testsuite')] for key in ('failures', 'errors', 'disabled', 'skipped'))
            if label == 'source-governance':
                assert cases[0].get('name') == 'OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims'
            else:
                expected = action['root_required_xml']
                assert cases[0].get('classname') == expected['suite'] and cases[0].get('name') == expected['name']
                assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
                if label in {'fp32', 'fp32comp', 'fp64'}:
                    props = {row.get('name'): row.get('value') for row in cases[0].findall('./properties/property')}
                    assert props.get(label + '_evidence') == 'real_adaptive_minimum_retry_controls_executed'
            recorded = next(row for row in result['cases'] if row['case'] == label)
            assert recorded['accepted'] and recorded['operator_exit'] == 0
            assert recorded['owned_wall_seconds'] == owner['wall_seconds'] and recorded['peak_sampled_rss_kib'] == owner['peak_sampled_rss_kib']
        owners[label] = {key: owner[key] for key in ('wall_seconds', 'peak_sampled_rss_kib', 'limits')}
    receipt = load(root / 'bin/linux-gcc/generated/sirius/alignment_receipt.json')
    assert receipt['source_revision'] == parent and receipt['alignment_mode'] == 'qualification'
    assert receipt['admitted_domains'] == receipt['records'] == []
    assert len(receipt['pending_domains']) == len(receipt['required_domains']) == 8
    assert receipt['ultimate_ideal']['satisfied'] is False
    assert 'Release packaging: DISABLED' in (work / 'configure/stdout.log').read_text()
    mapping = [(work, 'files/out/' + work.name)]
    for name in ['minimum-step-floor-registration', 'minimum-step-ellis-floor', 'macos-failure-diagnostics']:
        mapping.append((root / 'out' / name, 'files/out/' + name))
    build = root / 'bin/linux-gcc'
    products = {build / 'tests/backend/sirius_backend_tests', build / 'CMakeCache.txt', build / 'compile_commands.json',
        build / 'build.ninja', build / 'CMakeFiles/rules.ninja', build / 'generated/sirius/alignment_receipt.json',
        build / 'tests/source-governance-ctest.json', build / 'tests/evidence_governance.cmake'}
    products.update(build.rglob('CTestTestfile.cmake'))
    products.update((build / 'tests').rglob('*_tests.cmake'))
    products.update((build / 'tests').rglob('*_include.cmake'))
    products.update((build / 'generated/sirius').rglob('*.h'))
    for name in ['parity_probe.spv', 'parity_probe_fp32comp.spv', 'parity_probe_fp64.spv']:
        canonical = build / 'kernels' / name
        staged = build / 'tests/backend/resources/kernels' / name
        assert pin(canonical) == pin(staged)
        products.update([canonical, staged])
    products.update(build / 'kernels' / name for name in ['trace.spv', 'trace_fp32comp.spv', 'trace_fp64.spv', 'portability/trace.cu', 'portability/trace.metal'])
    mapping.extend((path, 'producer/' + str(path.relative_to(build))) for path in sorted(products))
    summary = {'schema': 'sirius-minimum-step-floor-focused-controls-v1', 'source_revision': parent,
        'result': result, 'owners': owners, 'static_input_count': 1033,
        'previous_configure_refusals': ['43ed3034d92a7b94a2d9a6711d8b2cda7a2e46cf', 'b220fa539ea6e3efd9161c6f04914d01b301e031'],
        'qualification_claimed': False, 'admitted_external_domains': 0, 'pending_external_domains': 8,
        'claim_limits': ['Live assertions and exact XML/properties passed; no raw float-buffer artifacts saved or byte-identity claim.',
            'Three precision controls cover the real adaptive attempt and central state/controller; all40 continuation fields and full frame publication are not established.',
            'Original Ellis finite ray passed separately; no Kerr frame success, CI-failure causality, allMandatory or full scientific/platform/release qualification.',
            'Build includes original source governance, Slang/SPIR-V validation and portability emissions; existing Slang initialization warnings remain.'],
        'prior_evidence_is_history_only': True}
    message = 'Preserve strict build and six exact-minimum focused controls (#53)\n'
elif kind == 'minimum-step-floor-alignment-stop':
    work = root / 'out/minimum-step-floor-verification'
    bundle = root / 'attestations/qualification/4498cc4/minimum-step-floor-alignment-stop'
    tag = 'evidence/issue-53-minimum-step-floor-alignment-stop'
    load = lambda p: json.loads(p.read_text())
    owner = load(work / 'configure/owner.json')
    assert owner['status'] == 'terminal' and not owner['passed'] and owner['child_exit'] == 1
    assert owner['source'] == {'revision': parent, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (work / 'configure/inputs-before.json').read_bytes() == (work / 'configure/inputs-after.json').read_bytes()
    assert not (work / 'build-backend').exists()
    stderr = ' '.join((work / 'configure/stderr.log').read_text().split())
    assert 'attested operating model differs from the exact source model' in stderr
    assert 'operating-model configuration failed' not in stderr
    assert load(work / 'alignment-refusal.json')['configuration_accepted'] is False
    mapping = [(work, 'files/out/' + work.name)]
    summary = {'schema':'sirius-minimum-step-stale-selected-attestation-refusal-v1', 'source_revision': parent,
        'owner': owner, 'configure_accepted': False, 'compilation_started': False, 'numerical_execution_started': False,
        'qualification_claimed': False, 'alignment_refusal':load(work / 'alignment-refusal.json'),
        'scope':'Corrected source/model policy accepted, selected cached older-model evidence rejected. Original frozen inputs exclude ignored selected attestations; post-refusal raw records preserved candidly, no complete external-input stability claimed.',
        'next_action':'Select no stale attestation input for existing qualification mode, preserving all8missing domains and all admission/release requirements; no unchanged input retry.'}
    message = 'Preserve fail-closed stale attestation rejection before minimum-step build (#53)\n'
elif kind == 'minimum-step-floor-configure-stop':
    work = root / 'out/minimum-step-floor-repair'
    bundle = root / 'attestations/qualification/27f5d9c/minimum-step-floor-configure-stop'
    tag = 'evidence/issue-53-minimum-step-floor-configure-stop'
    load = lambda p: json.loads(p.read_text())
    owner = load(work / 'configure/owner.json')
    assert owner['status'] == 'terminal' and not owner['passed'] and owner['child_exit'] == 1
    assert owner['source'] == {'revision': parent, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (work / 'configure/inputs-before.json').read_bytes() == (work / 'configure/inputs-after.json').read_bytes()
    for name, expected in load(work / 'configure/inputs-before.json').items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    for name, expected in load(work / 'source-freeze.json')['pins'].items():
        assert pin(root / name) == expected, name
    assert load(work / 'independent-review.json')['passed']
    assert not (work / 'build-backend').exists() and not (work / 'controls').exists()
    stderr = ' '.join((work / 'configure/stderr.log').read_text().split())
    assert 'semantic policy digest differs from verifier policy' in stderr
    assert 'enabled source GoogleTests are not Mandatory' in stderr
    assert 'evidence policy differs' in stderr
    assert (work / 'configure/stdout.log').read_text().endswith('-- Configuring incomplete, errors occurred!\n')
    mapping = [(work, 'files/out/' + work.name)]
    summary = {'schema':'sirius-minimum-step-floor-configure-stop-v1', 'source_revision': parent,
        'owner': owner, 'configure_accepted': False, 'compilation_started': False,
        'numerical_execution_started': False, 'qualification_claimed': False,
        'failure': 'Three new gtests were added to source/P1 model but omitted from explicit Mandatory registrations and the verifier P1 authority; semantic policy digest also stale.',
        'review_limit': 'Author, independent source reviewer and root missed the policy/registration joins. Earlier bounded source PASS does not admit configuration or runtime.',
        'next_action': 'Correct the existing Mandatory labels, verifier policy and legitimate semantic digest together; unchanged retries and gate relaxation prohibited.'}
    message = 'Preserve minimum-step repair configure rejection before compilation (#53)\n'
elif kind == 'multiply64-original-controls':
    import importlib.util
    work = root / 'out/portable-multiply64-controls-o0-lines'
    products = root / 'bin/linux-gcc/portable-multiply64-controls-o0-lines'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-original-controls'
    tag = 'evidence/issue-44-multiply64-original-controls'
    sys.path.insert(0, str(work))
    spec = importlib.util.spec_from_file_location('preserved_multiply64_original_controls', work / 'gate.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    gate.verify_static()
    packet = json.loads((work / 'packet.json').read_text())
    result = json.loads((work / 'prerequisites-result.json').read_text())
    assert result['passed'] and result['original_eight_controls_passed']
    assert not result['production_adopted'] and not result['benefit_executed']
    assert result['arithmetic'] == gate.verify_arithmetic()
    labels = [*packet['outputs'], *['control-' + str(index) for index in range(1, 9)]]
    assert len(labels) == len(set(labels)) == 18
    owners = {label: gate.accepted(label) for label in labels}
    runtime = json.loads((work / 'runtime-context.json').read_text())
    assert runtime['producer_output_seals'] == {label: gate.verify_output_seal(label) for label in packet['producer_labels']}
    for key, filename in [('profile', 'profile.json'), ('producer', 'producer-acceptance.json'),
                          ('feature', 'producer-feature-acceptance.json'), ('products', 'products-acceptance.json')]:
        assert runtime[key] == gate.record(work / filename)
    assert runtime['executable'] == gate.record(products / 'controls/original_candidate_controls')
    for index in range(1, 9):
        label = 'control-' + str(index)
        gate.join_runtime_inputs(label, owners[label])
        assert owners[label]['limits'] == packet['control_limits'][label]
    for label in ('verify-controls', 'verify-endpoint-routes'):
        gate.verify_output_seal(label)
    control = json.loads((work / 'control-acceptance.json').read_text())
    routes = json.loads((work / 'endpoint-route-acceptance.json').read_text())
    profile = json.loads((work / 'products-acceptance.json').read_text())
    assert control['passed'] and len(control['original_controls']) == 8 and routes['passed']
    assert control['timestamp_schema']['complete_readbacks'] == 32 and control['timestamp_schema']['ordered_markers'] == 51
    assert control['timestamp_schema']['actual_module_command_joins'] == 54
    assert len(routes['routes']) == 4 and routes['complete_readbacks'] == 56
    assert profile['passed'] and profile['changed_modules'] == 24 and profile['unchanged_modules'] == 11
    assert profile['main_arrays_changed'] == 24 and profile['main_arrays_unchanged'] == 16
    assert profile['owner_arrays_unchanged'] == 2 and profile['all_seven_tables_and_nonarray_ABI_constants_equal']
    reader = json.loads((work / 'parser-controls-result.json').read_text())
    assert reader['passed'] and len(reader['controls']) == 15 and sum(row['actual_acceptance'] for row in reader['controls']) == 2
    terminal = json.loads((work / 'independent-terminal-review.json').read_text())
    assert terminal['passed']
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-optional64-multiply-original-integrated-prerequisites-v1', 'result': result,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc', 'owners': owners,
               'parser_controls': reader, 'independent_terminal_review': pin(work / 'independent-terminal-review.json'), 'generated_profile': profile, 'control_acceptance': control, 'endpoint_routes': routes,
               'source_feature_scope': 'Private explicitly admitted Int64 diagnostic candidate; exact original stock producer preserved. Actual new ELF/header/module joins and numerical execution are distinct from historical public feature-argument capture reuse.',
               'stock_reuse_scope': '78 exact completed stock/candidate products and3accepted compiler-policy owners reused from immutable c420 archive. Failed old profile and failed fixture-construction attempts remain preserved. Fresh15parsercontrols, full35profile, actualproducer and8numerical controls are separately owned.',
               'performance_benefit_or_adoption_claimed': False, 'qualification_claimed': False,
               'next_dependency': 'Sole original fixed A1/B1/B2/A2 with three repeats per family and unchanged20+2 benefit gates. Production optional-capability/fallback admission and full qualification remain separate.'}
    message = 'Preserve multiply64 complete retained profile and original eight integrated controls (#44)\n'
elif kind == 'p1-kerr-observer-rss-stop':
    work = root / 'out/p1-kerr-minimum-step-observer'
    products = root / 'bin/linux-gcc/p1-kerr-minimum-step-observer'
    bundle = root / 'attestations/scientific-core/4bcd987/p1-kerr-minimum-step-observer-rss-stop'
    tag = 'evidence/issue-53-p1-kerr-minimum-step-observer-rss-stop'
    doc = lambda name: json.loads((work / name).read_text())
    summary = doc('incomplete-result.json')
    assert not summary['diagnostic_complete'] and not summary['qualification_claimed']
    freeze = doc('source-freeze.json')['pins']
    assert len(freeze) == 641
    for name, expected in freeze.items():
        assert pin(root / name) == {key:expected[key] for key in ('bytes','sha256')}, name
    for label in ('stage-trace-resource','compile-observer','link-observer','observe-original-kerr'):
        folder = work / label
        owner = doc(label + '/owner.json')
        assert owner['status'] == 'terminal' and owner['source'] == {'revision':parent,'status':''}
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in doc(label + '/inputs-before.json').items():
            assert pin(root / name) == {key:expected[key] for key in ('bytes','sha256')}, name
        if label == 'observe-original-kerr':
            assert not owner['passed'] and owner['child_exit'] == -15
            assert owner['stop_reason'] == 'diagnostic_sampled_rss_guard' and owner['peak_sampled_rss_kib'] == 4200640
        else:
            assert owner['passed'] and owner['child_exit'] == 0 and owner['stop_reason'] is None
    assert doc('compile-review.json')['passed'] and doc('link-review.json')['passed']
    assert doc('independent-terminal-review.json')['review_passed'] and not doc('independent-terminal-review.json')['diagnostic_complete']
    assert (work / 'observe-original-kerr/original-request-parameters.bin').read_bytes() == (work / 'original-request-parameters.expected.bin').read_bytes()
    assert not (work / 'observe-original-kerr/observation.json').exists() and not (work / 'observe-original-kerr/current-records.bin').exists()
    mapping = [(work,'files/out/' + work.name),(products,'files/bin/linux-gcc/' + products.name)]
    message = 'Preserve resource-limited original Kerr refusal observer without numerical claim (#53)\n'
elif kind == 'macos-current-provider-recovery':
    work = root / 'out/macos-current-provider-recovery'
    bundle = root / 'attestations/native-vulkan/5b5fad9/recovered-current-backend-modules'
    tag = 'evidence/issue-48-macos-current-backend-modules-recovery'
    doc = lambda name: json.loads((work / name).read_text())
    backend = work / 'sirius_backend_tests'
    assert pin(backend) == {'bytes':11560584,'sha256':'682b2ef3c91fe2ea46ab97a18096e509cab5d5eb3805a5b886e63acbe00577ad'}
    assert pin(work / 'native-provider-modules.json') == {'bytes':16451,'sha256':'c020b4c7b2e74979a2a094a9eaa7e0bd3afbaa839ffbf0bb125bf52d36acfef0'}
    modules = doc('native-provider-modules.json')['modules']
    assert len(modules) == 35
    for name, expected in modules.items():
        raw = backend.read_bytes()[expected['file_offset']:expected['file_offset'] + expected['bytes']]
        assert {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()} == {key:expected[key] for key in ('bytes','sha256')}
        assert (work / 'recovered-modules' / (name.replace('::','--') + '.spv')).read_bytes() == raw
    assert doc('module-sidecar-recovery.json')['passed'] and not doc('module-sidecar-recovery.json')['qualification_claimed']
    assert doc('independent-review.json')['passed']
    summary = {'schema':'sirius-exact-current-macos-backend-module-recovery-v1','recovery':doc('recovery.json'),
               'module_sidecar_recovery':doc('module-sidecar-recovery.json'),'independent_review':pin(work / 'independent-review.json'),
               'original_sidecar_digest_matched':True,'recovered_payloads':35,'native_provider_binding_sidecar_recovered':False,
               'only_pure_metadata_parser_executed_no_driver_shader_API_or_numerical_work':True,'qualification_claimed':False}
    mapping = [(work,'files/out/' + work.name)]
    message = 'Recover exact current macOS backend and original 35-module sidecar identity (#48)\n'
elif kind == 'macos-full-ci-incomplete-38577826600':
    work = root / 'out/native-provider-integration/full-ci-38577826600/macos-build'
    bundle = root / 'attestations/qualification/5b5fad9/macos-full-ci-incomplete-38577826600'
    tag = 'evidence/issue-48-macos-full-ci-incomplete-38577826600'
    doc = lambda name: json.loads((work / name).read_text())
    job, review, metadata = doc('job-api.json'), doc('independent-review.json'), doc('artifact-members.json')
    assert job['head_sha'] == '5b5fad91c08c0533e29c4f5c7ae249725162e405' and job['id'] == 115899499765
    assert job['status'] == 'completed' and job['conclusion'] == 'cancelled'
    assert review['passed'] and review['mandatory']['passed'] == 24 and review['mandatory']['completed'] == 25
    assert review['mandatory']['failed'] == 0 and review['mandatory']['timeout'] == 1
    assert review['mandatory']['at_cancellation']['id'] == 26 and not review['mandatory']['junit_present']
    assert pin(work / 'independent-review.json')['sha256'] == 'e766db58d19e41930545202bf9c8a277eedbfe40698b0672c77eacd63b174a8d'
    assert metadata['metadata']['id'] == 11924968602 and len(metadata['members']) == 12496
    assert pin(work / 'gate-diagnostics.zip') == {'bytes':93920417,'sha256':'0b67da9f14faed909b6051029e79ed6bf7b6ae6395c30523813a3bdf0ea4d6ad'}
    assert metadata['metadata']['digest'] == 'sha256:' + metadata['zip']['sha256']
    assert doc('artifact-transport.json')['accepted'] and doc('artifact-transport.json')['timeout_seconds'] == 300
    assert any('maximum execution time of 6h' in x['message'] for x in doc('check-run-annotations-api.json'))
    mapping = [(path, 'files/out/native-provider-integration/full-ci-38577826600/macos-build/' + path.name)
               for path in sorted(work.iterdir()) if path.is_file() and path.name != 'gate-diagnostics.zip']
    selected = review['selected_raw_member_bindings']
    assert len(selected) == 48
    for name, expected in selected.items():
        assert metadata['members'][name] == expected and pin(work / 'artifact' / name) == expected, name
        mapping.append((work / 'artifact' / name, 'selected-artifact/' + name))
    summary = {'schema':'sirius-macos-full-ci-incomplete-v1','source_revision':job['head_sha'],'job_id':job['id'],
               'run_id':job['run_id'],'completed_utc':job['completed_at'],'conclusion':'cancelled',
               'independent_audit':pin(work / 'independent-review.json'),'mandatory':review['mandatory'],
               'source_provider':review['source_provider'],'actual_provider_observations':review['actual_provider_observations'],
               'module_source_boundary':review['module_source_boundary'],'archive':metadata['zip'],
               'all_original_zip_members_recoverable':12496,'selected_raw_members_preserved':48,
               'archive_member_manifest':pin(work / 'artifact-members.json'),
               'qualification_claimed':False,'full_scientific_or_release_claimed':False,'unchanged_CI_retry_performed':False}
    parts = work / 'archive-parts'
    parts.mkdir()
    chunks = []
    with (work / 'gate-diagnostics.zip').open('rb') as stream:
        while block := stream.read(32 * 1024**2):
            part = parts / f'part-{len(chunks):03d}'
            part.write_bytes(block)
            name = 'artifact.zip.parts/' + part.name
            chunks.append({'path':name,**pin(part)})
            mapping.append((part,name))
    digest = hashlib.sha256()
    total = 0
    for row in chunks:
        raw = (parts / Path(row['path']).name).read_bytes()
        digest.update(raw); total += len(raw)
    assert {'bytes':total,'sha256':digest.hexdigest()} == metadata['zip']
    summary['archive_parts'] = chunks
    message = 'Preserve partial six-hour macOS runtime and exact diagnostic archive (#48)\n'
elif kind == 'linux-gate-incomplete-38577826600':
    work = root / 'out/native-provider-integration/full-ci-38577826600/linux-gate'
    bundle = root / 'attestations/qualification/5b5fad9/linux-gate-incomplete-38577826600'
    tag = 'evidence/issue-48-linux-gate-incomplete-38577826600'
    doc = lambda name: json.loads((work / name).read_text())
    job = doc('job-api.json')
    assert job['head_sha'] == '5b5fad91c08c0533e29c4f5c7ae249725162e405' and job['id'] == 115899499679
    assert job['status'] == 'completed' and job['conclusion'] == 'cancelled'
    assert pin(work / 'job-log.stdout') == {'bytes':124815,'sha256':'47d28bd65cc5a38e015a01c6dea10c1cea900202476156c32ddd8acfc1053a48'}
    assert pin(work / 'gate-diagnostics.zip') == {'bytes':42282,'sha256':'5545fa4a05902a266eb2d60e485772a011d44bcdf2d552a88de3ddd8589dd734'}
    artifact = doc('artifact-members.json')
    assert artifact['artifact_id'] == 11925267107 and len(artifact['members']) == 3
    assert artifact['metadata']['digest'] == 'sha256:' + artifact['zip']['sha256']
    for name, expected in artifact['members'].items():
        assert pin(work / 'artifact' / name) == expected, name
    review = doc('independent-review.json')
    assert review['passed'] and review['completed_results']['passed'] == 820
    assert review['completed_results']['completed'] == 824 and review['completed_results']['failed'] == 3
    assert review['completed_results']['timeout'] == 1 and not review['completed_results']['junit_present']
    assert review['elapsed_work']['largest_completed_case']['seconds'] == 15949.58
    assert any('maximum execution time of 6h' in x['message'] for x in doc('check-run-annotations-api.json'))
    summary = {'schema':'sirius-linux-full-ci-incomplete-v1','source_revision':job['head_sha'],'job_id':job['id'],
               'run_id':job['run_id'],'completed_utc':job['completed_at'],'conclusion':'cancelled',
               'independent_audit':pin(work / 'independent-review.json'),'observed_completed_results':review['completed_results'],
               'nonpasses':review['nonpasses'],'elapsed_work':review['elapsed_work'],
               'cause_boundary':review['cancellation']['cause_boundary'],'identity_limits':review['identity']['identity_limit'],
               'qualification_claimed':False,'full_scientific_or_release_claimed':False,'unchanged_retry_performed':False}
    mapping = [(work, 'files/out/native-provider-integration/full-ci-38577826600/linux-gate')]
    message = 'Preserve incomplete six-hour Linux Mandatory run and four observed nonpasses (#48)\n'
elif kind == 'native-windows-system-5b5fad9':
    work = root / 'out/native-windows-system-5b5fad9'
    bundle = root / 'attestations/native-vulkan/5b5fad9/fresh-windows-system-inventory'
    tag = 'evidence/issue-48-native-windows-system-5b5fad9'
    summary = json.loads((work / 'finite-result.json').read_text())
    doc = lambda name: json.loads((work / name).read_text(encoding='utf-8-sig'))
    assert summary['passed'] and not summary['qualification_claimed'] and not summary['runtime_domain_admitted']
    binding = doc('execution-bindings.json')
    assert binding['runner_source_revision'] == parent and binding['source_revision'] == '5b5fad91c08c0533e29c4f5c7ae249725162e405'
    assert len(binding['files']) == 115
    for name, expected in binding['files'].items():
        assert pin(root / name) == expected, name
    assert (work / 'inputs-before.json').read_bytes() == (work / 'inputs-after.json').read_bytes()
    before = doc('inputs-before.json')
    assert len(before) == 117
    for name, expected in before.items():
        assert pin(root / name) == expected, name
    owner, bridge = doc('native/inventory/owner.json'), doc('inventory-bridge.json')
    assert owner['status'] == 'completed' and owner['returncode'] == 0 and not owner['outer_stop'] and not owner['cleanup_errors']
    assert owner['system_inventory_passed'] and not owner['numerical_control_passed'] and owner['owned_process_absent']
    assert owner['source_build_gate_unchanged'] and owner['test_executable_unchanged']
    assert owner['outer_guard_seconds'] == 60 and owner['rss_guard_bytes'] == 4294967296
    assert owner['source_revision'] == binding['source_revision'] and owner['live_source_revision'] == parent
    assert bridge['accepted_owner_execution'] and bridge['whole_inputs_unchanged'] and bridge['linux_child_birth_absent'] and not bridge['cleanup_errors']
    terminal, audit = doc('terminal-inventory.json'), doc('terminal-inventory-bridge.json')
    assert audit['accepted'] and audit['returncode'] == 0 and audit['linux_child_birth_absent'] and not audit['errors']
    assert terminal['owned_births_absent'] and not terminal['actual_loaded_modules_observed'] and not terminal['loaded_module_post_seals']
    assert not terminal['consumer_scan']['matches'] and not terminal['consumer_scan']['accessible_compilers']
    assert doc('terminal-inventory-accepted.json')['pass']
    assert doc('independent-terminal-review.json')['passed']
    prior = summary['archive_reuse']
    assert subprocess.check_output(['git','rev-parse','refs/tags/' + prior['tag'] + '^{commit}'],text=True).strip() == prior['evidence_commit']
    previous = root / prior['bundle']
    assert pin(previous / 'manifest.json')['sha256'] == prior['manifest_sha256']
    for name, expected in json.loads((previous / 'manifest.json').read_text())['files'].items():
        assert pin(previous / name) == expected, name
        assert subprocess.check_output(['git','show',prior['evidence_commit'] + ':' + str((previous / name).relative_to(root))]) == (previous / name).read_bytes()
    mapping = [(path, 'files/out/' + work.name + '/' + str(path.relative_to(work)))
               for path in sorted(work.rglob('*')) if path.is_file() and path.name != 'export.zip']
    summary['producer_recovery'] = {'recovery_receipt': pin(work / 'export-recovery.json'), 'reused_bundle_all_Git_bytes_checked': True}
    message = 'Preserve one fresh native Windows Radeon inventory on exact archived producer (#48)\n'
elif kind == 'multiply64-o0-parser-fixture-stop':
    work = root / 'out/portable-multiply64-controls-o0'
    products = root / 'bin/linux-gcc/portable-multiply64-controls-o0'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-o0-parser-fixture-stop'
    tag = 'evidence/issue-44-multiply64-o0-parser-fixture-stop'
    summary = json.loads((work / 'parser-control-stop.json').read_text())
    assert not summary['complete_profile_admitted'] and not summary['producer_or_numerical_controls_executed']
    assert summary['literal_substring_count'] == 2 and summary['exact_instruction_line_count'] == 1
    assert not (work / 'parser-controls-result.json').exists() and not (work / 'verify-products').exists()
    assert not (work / 'prepare-producer').exists() and not any((work / f'control-{n}').exists() for n in range(1, 9))
    for name, expected in json.loads((work / 'source-freeze.json').read_text())['pins'].items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    for label in ('reuse-completed-products', 'check-o0-reader'):
        folder = work / label
        owner = json.loads((folder / 'owner.json').read_text())
        assert owner['status'] == 'terminal' and owner['source'] == {'revision': parent, 'status': ''}
        assert owner['passed'] == (label == 'reuse-completed-products')
        assert owner['child_exit'] == (0 if label == 'reuse-completed-products' else 1)
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in json.loads((folder / 'inputs-before.json').read_text()).items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    seal = json.loads((work / 'reuse-completed-products/output-seal.json').read_text())
    for name, expected in seal['outputs'].items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    assert not (work / 'check-o0-reader/output-seal.json').exists()
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    message = 'Preserve stopped multiply64 parser fixture before source-reader execution (#44)\n'
elif kind == 'multiply64-profile-graph-stop':
    work = root / 'out/portable-multiply64-controls'
    products = root / 'bin/linux-gcc/portable-multiply64-controls'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-profile-graph-stop'
    tag = 'evidence/issue-44-multiply64-profile-graph-stop'
    summary = json.loads((work / 'compiler-screen-stop.json').read_text())
    assert not summary['complete_profile_admitted'] and not summary['producer_or_numerical_controls_executed']
    assert summary['complete_candidate_35module_builder_passed']
    assert not (work / 'products-acceptance.json').exists() and not (work / 'prepare-producer').exists()
    assert not any((work / f'control-{n}').exists() for n in range(1, 9))
    for name, expected in json.loads((work / 'source-freeze.json').read_text())['pins'].items():
        assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    for label in ('reuse-products-stock', 'compile-products-candidate', 'native-policy-selfcheck', 'verify-products'):
        folder = work / label
        owner = json.loads((folder / 'owner.json').read_text())
        assert owner['status'] == 'terminal' and owner['source'] == {'revision': parent, 'status': ''}
        assert owner['child_exit'] == (1 if label == 'verify-products' else 0)
        assert owner['passed'] == (label != 'verify-products')
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        for name, expected in json.loads((folder / 'inputs-before.json').read_text()).items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
        if label != 'verify-products':
            seal = json.loads((folder / 'output-seal.json').read_text())
            for name, expected in seal['outputs'].items():
                assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    assert pin(products / 'candidate/retained_camera_fp64.spv')['sha256'] == summary['candidate_module_sha256']
    assert pin(products / 'analysis/candidate/retained_camera_fp64.spvasm')['sha256'] == summary['candidate_assembly_sha256']
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    message = 'Preserve complete multiply64 compiler products and stopped native O0 graph reader (#44)\n'
elif kind == 'multiply64-original-arithmetic':
    import importlib.util
    work = root / 'out/portable-multiply64-arithmetic'
    products = root / 'bin/linux-gcc/portable-multiply64-arithmetic'
    bundle = root / 'attestations/software-vulkan/4bcd987/multiply64-original-arithmetic'
    tag = 'evidence/issue-44-multiply64-original-arithmetic'
    sys.path.insert(0, str(work))
    spec = importlib.util.spec_from_file_location('preserved_multiply64_arithmetic', work / 'gate.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    gate.verify_static()
    load = lambda path: json.loads(path.read_text())
    packet = load(work / 'packet.json')
    result = load(work / 'arithmetic-result.json')
    runtime = load(work / 'runtime-context.json')
    assert result['passed'] and result['schema'] == 'sirius-multiply64-original-arithmetic-result-v1'
    assert load(work / 'independent-source-review.json')['passed']
    assert load(work / 'independent-terminal-review.json')['passed']
    assert result['original_records'] == 174268 and result['original_operation_count'] == 14
    assert result['complete_output_words_per_provider'] == 697072
    assert result['stock_result_reused'] and result['host_explicit64_and_both_full_gpu_outputs_equal_original_expected']
    assert result['candidate_macro'] == 'SIRIUS_PB32_MULTIPLY64=1'
    owners = {label: gate.accepted(label) for label in [*packet['outputs'], 'execute-host', 'execute-gpu-candidate']}
    assert len(owners) == 10
    for label, owner in owners.items():
        assert owner['source'] == {'revision': parent, 'status': ''}
        assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
        assert result['owners'][label] == gate.record(work / label / 'owner.json')
        for name, expected in load(work / label / 'inputs-before.json').items():
            assert pin(root / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    assert runtime['all_build_output_seals'] == {label: gate.verify_output_seal(label) for label in packet['outputs']}
    assert result['modules'] == runtime['modules'] == gate.module_checks()
    assert result['stock_arithmetic_receipt'] == gate.verify_stock_reuse()
    assert result['compiler_and_jit_receipts'] == gate.verify_reference_acceptances()
    assert result['original_corpus'] == gate.corpus_checks()
    for label in ('execute-host', 'execute-gpu-candidate'):
        gate.join_runtime_inputs(label, owners[label])
    gate.provider_mappings(owners['execute-gpu-candidate'])
    assert (work / 'execute-gpu-candidate/actual.bin').read_bytes() == (work / 'expected.bin').read_bytes()
    assert pin(work / 'expected.bin')['bytes'] == 697072 * 4
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-multiply64-original-arithmetic-corpus-v1', 'result': result,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc', 'owners': owners,
               'scope': 'Unchanged original 174268-record14-operation finite corpus, fresh candidate host and complete697072-word GPU output. Exact historical stock arithmetic/source/build/producer/provider receipt reused. No fresh stock dispatch or numeric-process feature-argument capture claimed.',
               'integrated_retained_programs_or_native_FP64_residual_profiles_admitted': False,
               'performance_benefit_or_adoption_claimed': False, 'qualification_claimed': False}
    message = 'Preserve unchanged full original arithmetic gate passing multiply64 candidate (#44)\n'
elif kind == 'divide64-original-arithmetic':
    import importlib.util
    work = root / 'out/portable-divide64-arithmetic'
    products = root / 'bin/linux-gcc/portable-divide64-arithmetic'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide64-original-arithmetic'
    tag = 'evidence/issue-44-divide64-original-arithmetic'
    spec = importlib.util.spec_from_file_location('preserved_divide64_arithmetic', work / 'gate.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    gate.verify_static()
    result = json.loads((work / 'arithmetic-result.json').read_text())
    assert result['passed'] and result['original_records'] == 174268
    assert result['complete_output_words_per_provider'] == 697072
    assert result['host_explicit64_and_both_full_gpu_outputs_equal_original_expected']
    packet = json.loads((work / 'packet.json').read_text())
    owners = {label: gate.accepted(label) for label in result['owners']}
    assert len(owners) == 17
    for label, expected in result['owners'].items():
        assert gate.record(work / label / 'owner.json') == expected
    runtime = json.loads((work / 'runtime-context.json').read_text())
    assert result['runtime_context'] == gate.record(work / 'runtime-context.json')
    assert runtime['build_review'] == gate.record(work / 'build-review.json')
    assert runtime['all_build_output_seals'] == {label: gate.verify_output_seal(label) for label in packet['outputs']}
    assert runtime['gpu'] == gate.record(products / 'gpu') and runtime['host'] == gate.record(products / 'host')
    for label in ('execute-host', 'execute-gpu-stock', 'execute-gpu-candidate'):
        gate.join_runtime_inputs(label, owners[label])
    for label in ('execute-gpu-stock', 'execute-gpu-candidate'):
        gate.provider_mappings(owners[label])
        assert (work / label / 'actual.bin').read_bytes() == (work / 'expected.bin').read_bytes()
    assert (products / 'stock.spv').read_bytes() == Path(packet['original_stock_module']).read_bytes()
    assert result['original_corpus'] == gate.corpus_checks()
    assert runtime['same_admitted_host_reuse'] == gate.verify_host_reuse()
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-optional64-original-arithmetic-corpus-v1', 'result': result,
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'owners': owners, 'complete_output_bytes_per_shader': pin(work / 'expected.bin'),
               'same_admitted_host_reuse': runtime['same_admitted_host_reuse'],
               'initial_unconsumed_proposal': 'Preserved before prospective repairs to direct disassembly input seals and final numerical-consumer input joins; no owner executed the proposal.',
               'arithmetic_scope': 'Original unchanged finite 174268-record14-operation corpus: candidate C++ host plus stock and optional64 selected-provider shaders, complete697072 output words per shader.13674 divide records,12263 finite candidate branch records. This is finite primitive arithmetic coverage, not exhaustive proof.',
               'integrated_retained_programs_tables_or_normal_sum_accepted': False,
               'performance_benefit_or_adoption_claimed': False, 'qualification_claimed': False,
               'next_dependency': 'Original eight integrated controls with exact optional64 module/producer/field/table/mode joins; original unchanged20+2 adoption gates remain required.'}
    message = 'Preserve unchanged original arithmetic corpus passing optional64 host and Vulkan candidate (#44)\n'
elif kind == 'divide64-compiler-jit-screen':
    import importlib.util
    work = root / 'out/portable-divide64-jit-preparation'
    screen_work = root / 'out/portable-divide64-screen-capability-order'
    pair_work = root / 'out/portable-divide64-jit-pair-review'
    feature_work = root / 'out/portable-divide64-feature-admission'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide64-compiler-jit-screen'
    tag = 'evidence/issue-44-divide64-compiler-jit-screen'
    def module(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        value = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(value)
        return value
    screen = module('preserved_divide64_compiler_screen', screen_work / 'screen.py')
    screen.verify_frozen()
    compiler_acceptance = json.loads((screen_work / 'compiler-acceptance.json').read_text())
    assert compiler_acceptance['passed']
    for name, filename in [('compiler_owner', 'compile-endpoint-normal/owner.json'), ('compiler_result', 'compiler-result.json')]:
        assert screen.record(screen_work / filename) == compiler_acceptance[name]
    compiler_owner = json.loads((screen_work / 'compile-endpoint-normal/owner.json').read_text())
    assert compiler_owner['passed'] and compiler_owner['child_exit'] == 0 and compiler_owner['stop_reason'] is None
    assert compiler_owner['source_unchanged'] and compiler_owner['whole_inputs_unchanged']
    assert not compiler_owner['cleanup_errors'] and not compiler_owner['remaining_owned_processes']
    assert compiler_owner['observed_births'] and all(row['absent'] for row in compiler_owner['observed_births'])
    host = module('preserved_divide64_host', work / 'host.py')
    host.verify_static()
    owners = {label: host.accepted(label) for label in
              ['compile-device', 'compile-probe', 'link', 'capture-stock-endpoint-normal', 'capture-candidate-endpoint-normal']}
    for label in ['compile-device', 'compile-probe']:
        host.verify_compile_outputs(label)
    pair = json.loads((pair_work / 'pair-review.json').read_text())
    assert pair['passed'] and pair['complete_55_enabled_fields_queue_extensions_order_and_feature_chain_equal']
    assert pair['reader'] == host.record(pair_work / 'review_pair.py') and pair['consumed_host'] == host.record(work / 'host.py')
    assert len(pair['negative_controls']) == 6 and all(row['rejected'] for row in pair['negative_controls'])
    for label, records in pair['captures'].items():
        folder = work / label
        actual = {str(path.relative_to(folder)): host.record(path) for path in sorted(folder.rglob('*')) if path.is_file()}
        assert actual == records
    controls = json.loads((root / 'out/portable-divide64-policy-controls.json').read_text())
    assert controls['passed'] and len(controls['controls']) == 12
    assert all(row['expected_acceptance'] == row['actual_acceptance'] for row in controls['controls'])
    assert pin(screen_work / 'candidate/scripts/build-retained-kernels.py')['sha256'] == controls['builder_sha256']
    feature = json.loads((feature_work / 'direct-query.json').read_text())
    assert feature['features']['shaderInt64'] == 1 and feature['before'] == feature['after']
    for path, expected in feature['after'].items():
        assert host.record(path) == expected
    independent = json.loads((pair_work / 'independent-reviews.json').read_text())
    assert independent['pair_reader']['accepted'] and independent['lowering_cost']['premise_survives_for_original_arithmetic_checks']
    assert independent['pair_reader']['reader'] == pin(pair_work / 'review_pair.py')
    assert independent['pair_reader']['receipt'] == pin(pair_work / 'pair-review.json')
    mapping = [(path, 'files/' + str(path.relative_to(root))) for path in
               [screen_work, root / 'bin/linux-gcc/portable-divide64-screen-capability-order',
                work, root / 'bin/linux-gcc/portable-divide64-jit', pair_work, feature_work,
                root / 'out/portable-divide64-policy-controls.py', root / 'out/portable-divide64-policy-controls.json']]
    summary = {'schema': 'sirius-optional64-compiler-native-jit-screen-v1',
               'product_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'original_stopped_compiler_attempt': {'commit': 'da2338b6f2a803c9aa06801a3699163b25f80d6f',
                   'tag': 'evidence/issue-44-divide64-compiler-capability-order'},
               'compiler_acceptance': compiler_acceptance, 'compiler_owner': compiler_owner,
               'five_terminal_host_and_capture_owners': owners, 'corrected_pair_review': pair,
               'candidate_builder_policy_controls': controls, 'independent_reviews': independent,
               'pair_reader_correction': 'Original consumed host compares integer-keyed Python relocation counts to string-keyed durable JSON. Separate reviewed pure reader compares complete durable JSON, retaining all fields, value types, sequences and six integrity negatives; neither capture was rerun.',
               'device_execution_scope': 'Both same-feature-host original one-word short-input refusal captures; complete Endpoint normal modules loaded and actual callable native code joined, but retained arithmetic and generated tables did not execute.',
               'stock_builder_and_original_32bit_fallback_preserved': True,
               'next_dependency': 'Original unchanged 174268-record arithmetic provider gate, then eight original integrated controls and original unchanged 20+2 benefit/adoption gates.',
               'retained_arithmetic_or_performance_benefit_or_adoption_claimed': False,
               'qualification_claimed': False}
    message = 'Preserve optional64 compiler and called native JIT screen without arithmetic or benefit claims (#44)\n'
elif kind == 'divide64-compiler-capability-order':
    import re
    work = root / 'out/portable-divide64-screen'
    products = root / 'bin/linux-gcc/portable-divide64-screen'
    bundle = root / 'attestations/software-vulkan/4bcd987/divide64-compiler-capability-order'
    tag = 'evidence/issue-44-divide64-compiler-capability-order'
    owner = json.loads((work / 'compile-endpoint-normal/owner.json').read_text())
    assert owner['status'] == 'terminal' and not owner['passed'] and owner['child_exit'] == 1 and owner['stop_reason'] is None
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert 'portable retained stage introduced an optional capability or non-32-bit integer' in (work / 'compile-endpoint-normal/stderr.log').read_text()
    assembly = (products / 'modules/candidate-endpoint-normal.spvasm').read_text()
    capabilities = re.findall(r'OpCapability (\S+)', assembly)
    assert capabilities == ['Int64', 'Shader']
    assert set(re.findall(r'OpTypeInt (\d+) [01]', assembly)) == {'32', '64'}
    assert re.findall(r'OpTypeFloat (\d+)', assembly) == ['32']
    stock = products / 'modules/stock-endpoint-normal.spv'
    canonical = root / 'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv'
    assert stock.read_bytes() == canonical.read_bytes()
    assert pin(stock) == {'bytes': 151400, 'sha256': '981bc85d13e168c6adf72c99adfe169157593dda1bee46e78ff867fc6b199ede'}
    assert not (work / 'compiler-result.json').exists() and not (work / 'compiler-acceptance.json').exists()
    mapping = [(work, 'files/out/' + work.name), (products, 'files/bin/linux-gcc/' + products.name)]
    summary = {'schema': 'sirius-divide64-compiler-capability-order-diagnostic-v1', 'owner': owner,
               'product_source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'stock_complete_module_equals_canonical': pin(stock),
               'observed_candidate_capability_order': capabilities,
               'observed_candidate_integer_widths': [32, 64],
               'failure': 'Private candidate-only guard expected declarations in Shader,Int64 order; actual Slang emits Int64,Shader. SPIR-V permits capability declaration order. Guard correction must retain the exact multiset and all original stock/type/float/layout/barrier/numerical gates.',
               'reference': 'https://registry.khronos.org/SPIR-V/specs/unified1/SPIRV.html#OpCapability',
               'candidate_final_validation_completed': False,
               'device_or_retained_arithmetic_executed': False,
               'performance_benefit_or_adoption_claimed': False,
               'qualification_claimed': False}
    message = 'Preserve stopped optional Int64 compiler screen and capability-order guard diagnosis (#44)\n'
elif kind == 'p4-doppler-scene2':
    work = root / 'out/science-p4-doppler-scene2-outer1200'
    rendered = root / 'renders/p4-doppler-scene2-outer1200'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-doppler-scene2'
    tag = 'evidence/issue-53-p4-doppler-scene2'
    result = json.loads((work / 'scene-2-result.json').read_text())
    owner = json.loads((work / 'scene-2/owner.json').read_text())
    prepared = json.loads((work / 'preparation.json').read_text())
    assert len(prepared['pins']) == 503
    assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert owner['limits'] == {'seconds': 1200, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    assert not result['original_assertions_changed'] and not result['qualification_claimed'] and not result['p4_complete']
    for name, frozen in prepared['pins'].items():
        path = root / name
        assert pin(path) == {k: frozen[k] for k in ('bytes', 'sha256')}, name
        assert str(path.resolve(strict=True)) == frozen['resolved_path'], name
    prefix = '[Vulkan] Render evidence: '
    frames = [json.loads(line.split(prefix, 1)[1]) for line in
              (work / 'scene-2/stdout.log').read_text().splitlines() if prefix in line]
    if result['passed']:
        assert owner['passed'] and owner['child_exit'] == 0 and owner['stop_reason'] is None
        assert len(frames) == 4 and result['complete_original_frames'] == frames
    else:
        assert not owner['passed'] and 'error' in result
    assert rendered.is_dir() and not rendered.is_symlink() and rendered.resolve(strict=True) == rendered
    outputs = list(rendered.iterdir())
    assert all(p.is_file() and not p.is_symlink() and p.name in prepared['cases'][0]['outputs'] for p in outputs)
    mapping = [(work, 'files/out/' + work.name), (rendered, 'files/renders/' + rendered.name),
               (root / 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'producer/sirius_render_tests'),
               (root / 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    summary = {'schema': 'sirius-original-p4-doppler-scene2-v1', 'result': result, 'owner': owner,
               'product_source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'whole_complete_frame_observations': frames, 'observed_complete_frames': len(frames),
               'original_full_four_frame_test_passed': result['passed'],
               'prior_scene1_incomplete_record_unchanged': True,
               'outer_diagnostic_policy': 'First originally unlaunched scene2, finite1200s private containment for the original four40x24/1SPP thin/volume Doppler frames; original submission feedback1000ms/ray attempts20000/precision/cache/source/numerical assertions unchanged. Required cold128x128/1SPP Kerr300s remains unchanged and unresolved.',
               'performance_or_optimization_adoption_claimed': False,
               'full_p4_science_mandatory_domain_or_release_qualification_claimed': False}
    message = 'Preserve complete raw scope and terminal disposition of original four-frame P4 Doppler scene (#53)\n'
elif kind in ('p4-ordinary-cases5-9', 'p4-ordinary-cases3-4-completion'):
    import math
    import xml.etree.ElementTree as ET
    completion = kind == 'p4-ordinary-cases3-4-completion'
    suffix = 'cases3-4-completion' if completion else 'cases5-9'
    work = root / ('out/science-p4-ordinary-fp32-' + suffix + '-outer600')
    bundle = root / ('attestations/scientific-core/4bcd987/p4-ordinary-fp32-' + suffix)
    tag = 'evidence/issue-53-p4-ordinary-fp32-' + suffix
    prepared = json.loads((work / 'preparation.json').read_text())
    assert len(prepared['pins']) == (585 if completion else 510)
    for name, frozen in prepared['pins'].items():
        path = root / name
        assert pin(path) == {k: frozen[k] for k in ('bytes', 'sha256')}, name
        assert str(path.resolve(strict=True)) == frozen['resolved_path'], name
    results, owners = {}, {}
    for case in prepared['cases']:
        label = case['label']
        result = json.loads((work / (label + '-result.json')).read_text())
        owner = json.loads((work / label / 'owner.json').read_text())
        assert result['passed'] and result['ordinary_fp32'] and result['name'] == case['name']
        assert not result['original_assertions_changed'] and not result['qualification_claimed'] and not result['p4_complete']
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0 and owner['stop_reason'] is None
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['limits'] == {'seconds': 600, 'rss_kib': 4194304, 'reserve_kib': 2097152}
        assert owner['wall_seconds'] <= 600 and owner['peak_sampled_rss_kib'] <= 4194304 and owner['minimum_sampled_available_kib'] >= 2097152
        for filename, key in (('owner.json', 'owner'), ('gtest.xml', 'gtest'), ('ctest.xml', 'ctest')):
            path = work / label / filename
            assert result[key] == dict(pin(path), resolved_path=str(path.resolve(strict=True)))
        for filename in ('gtest.xml', 'ctest.xml'):
            tree = ET.parse(work / label / filename).getroot()
            assert len(tree.findall('.//testcase')) == 1
            for suite in [tree, *tree.findall('.//testsuite')]:
                assert all(int(suite.get(key, '0')) == 0 for key in ('failures', 'errors', 'skipped', 'disabled'))
            assert not any(tree.findall('.//testcase/' + key) for key in ('failure', 'error', 'skipped'))
        props = result['original_properties']; prefix = 'parity_opcode_' + str(case['opcode'])
        assert props[prefix + '_completed_records'] == str(case['records']) and props[prefix + '_actual_submissions'] == '1'
        assert case['minimum_allocated_bytes'] <= int(props[prefix + '_allocated_bytes']) <= case['maximum_allocated_bytes']
        actual_ms = float(props[prefix + '_maximum_submit_ms'])
        assert math.isfinite(actual_ms) and 0 < actual_ms <= 1000
        assert len(result['original_dispatch_phases']) == 2 and result['no_performance_or_cache_warmness_claim']
        results[label], owners[label] = result, owner
    assert len(results) == (2 if completion else 5)
    mapping = [(work, 'files/out/' + work.name),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'producer/sirius_backend_tests'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    summary = {'schema': 'sirius-original-p4-ordinary-fp32-cases5-9-v1', 'results': results, 'owners': owners,
               'product_source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'original_case_ordinals': [5, 6, 7, 8, 9], 'original_cases_passed': 5,
               'previous_cases_1_and_2_passed': 'evidence/issue-53-p4-ordinary-fp32-isolated-partial',
               'original_cases_3_and_4_remain_incomplete_and_not_retried': True,
               'source_assertions_precision_provider_cache_and_dispatch_gates_unchanged': True,
               'outer_diagnostic_policy': 'First previously unlaunched cases5..9: finite600s whole-case diagnostic owner, distinct from original<=1000ms actual sample submissions; whole cold setup/teardown costs reported. Ordinary parity registration has no whole-case300s assertion. Required cold128x128/1SPP Kerr300s remains unchanged and unresolved.',
               'performance_or_optimization_adoption_claimed': False,
               'full_p4_science_mandatory_domain_or_release_qualification_claimed': False}
    message = 'Preserve five original ordinary FP32 P4 parity witnesses and full cold setup costs (#53)\n'
    if completion:
        summary = {'schema': 'sirius-original-p4-ordinary-fp32-cases3-4-completion-v1',
                   'results': results, 'owners': owners,
                   'product_source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
                   'original_case_ordinals': [3, 4], 'original_cases_passed': 2,
                   'second_attempts_of_original_scientific_inputs': True,
                   'prior_incomplete_attempts_preserved_unchanged': {
                       'case3': 'evidence/issue-53-p4-ordinary-fp32-isolated-partial at0d8beec5fdb1406d35add385cfa6f2d9f62848d2',
                       'case4': 'evidence/issue-53-p4-opacity-setup-incomplete atc901f1bcdf71533efbe61143b9437588ab37efdd'},
                   'previous_other_original_passes': {
                       'cases1_2': 'evidence/issue-53-p4-ordinary-fp32-isolated-partial',
                       'cases5_9': 'evidence/issue-53-p4-ordinary-fp32-cases5-9 at76a7d0281c1a47844d86170914dc18f72e27769a'},
                   'prospective_plan_revision': 'https://github.com/sunyaisblank/Sirius/issues/53#issuecomment-6172375641',
                   'source_assertions_precision_provider_cache_and_dispatch_gates_unchanged': True,
                   'outer_diagnostic_policy': 'Exactly one new completion opportunity per originalcase3/4 with finite600s private whole-case owner, after source audit established ordinary parity has no whole-case300s timeout and five matching cold controls completed. Actualsample<=1000ms, allocations/numerical/ownership/refusal unchanged. Full old and new wall costs retained; no causal speedup or coldness inference. Required cold128x128/1SPP Kerr300s remains unchanged and unresolved.',
                   'performance_or_optimization_adoption_claimed': False,
                   'full_p4_science_mandatory_domain_or_release_qualification_claimed': False}
        message = 'Preserve original Shakura-Sunyaev and opacity parity completion with prior timeout history (#53)\n'
elif kind == 'p4-opacity-setup-incomplete':
    work = root / 'out/science-p4-ordinary-fp32-parity-case4'
    held = root / 'out/science-p4-ordinary-fp32-parity-cases5-9'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-opacity-setup-incomplete'
    tag = 'evidence/issue-53-p4-opacity-setup-incomplete'
    result = json.loads((work / 'parity-4-result.json').read_text())
    owner = json.loads((work / 'parity-4/owner.json').read_text())
    assert not result['passed'] and not result['qualification_claimed'] and not result['original_assertions_changed']
    assert owner['status'] == 'terminal' and not owner['passed'] and owner['child_exit'] == -15
    assert owner['stop_reason'] == 'diagnostic_timeout' and 300 <= owner['wall_seconds'] < 310
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert all(row['absent'] for row in owner['observed_births'])
    assert '[ParityProbe] phase=' not in (work / 'parity-4/stdout.log').read_text()
    assert not list((work / 'parity-4').glob('*.xml'))
    assert all(not (held / ('parity-' + str(i))).exists() and not (held / ('parity-' + str(i) + '-result.json')).exists() for i in range(5, 10))
    preparation = json.loads((held / 'preparation.json').read_text())
    assert [c['records'] for c in preparation['cases']] == [5, 4, 2, 5, 2]
    assert len(preparation['pins']) == 470
    mapping = [(work, 'files/out/' + work.name), (held, 'files/out/' + held.name),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'producer/sirius_backend_tests'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    summary = {'schema': 'sirius-original-p4-opacity-setup-incomplete-v1', 'result': result, 'owner': owner,
               'product_source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
               'original_case_4_previously_unlaunched': True, 'original_case_3_unchanged_and_incomplete': True,
               'previous_two_passes_and_case3': 'evidence/issue-53-p4-ordinary-fp32-isolated-partial at0d8beec5fdb1406d35add385cfa6f2d9f62848d2',
               'held_cases_5_to_9': {'prepared_only': True, 'runtime_actions': 0, 'current_pins_independently_checked': 470},
               'numerical_verdict': None, 'setup_completed': False, 'original_assertions_and_limits_unchanged': True,
               'loop_policy_candidate_selected': False, 'performance_or_adoption_claimed': False,
               'full_p4_science_mandatory_domain_or_release_qualification_claimed': False}
    message = 'Preserve original opacity setup timeout and five unlaunched P4 witness preparations (#53)\n'
elif kind == 'p4-isolated-partial':
    work = root / 'out/science-p4-ordinary-fp32-parity-isolated'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-ordinary-fp32-isolated-partial'
    tag = 'evidence/issue-53-p4-ordinary-fp32-isolated-partial'
    results, owners = {}, {}
    for ordinal in (1, 2, 3):
        label = 'parity-' + str(ordinal)
        result = json.loads((work / (label + '-result.json')).read_text())
        owner = json.loads((work / label / 'owner.json').read_text())
        assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes'] and all(row['absent'] for row in owner['observed_births'])
        if ordinal < 3:
            assert result['passed'] and owner['passed'] and owner['child_exit'] == 0 and owner['stop_reason'] is None
            assert not result['original_assertions_changed'] and result['ordinary_fp32']
        else:
            assert not result['passed'] and not owner['passed'] and owner['stop_reason'] == 'diagnostic_timeout' and owner['child_exit'] == -15
            assert '[ParityProbe] phase=' not in (work / label / 'stdout.log').read_text()
            assert not (work / label / 'gtest.xml').exists() and not (work / label / 'ctest.xml').exists()
        results[label] = result; owners[label] = owner
    assert all(not (work / ('parity-' + str(i))).exists() for i in range(4, 10))
    mapping = [(work, 'files/out/science-p4-ordinary-fp32-parity-isolated'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'producer/sirius_backend_tests'),
               (root / 'bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    summary = {'schema': 'sirius-original-p4-isolated-partial-prefix-v1', 'results': results, 'owners': owners,
               'original_cases_passed': 2, 'original_cases_incomplete': 1, 'following_cases_unlaunched': 6,
               'source_assertions_precision_provider_cache_and_gates_unchanged': True,
               'allocation_reader_uses_original_inequality': True,
               'causal_speedup_or_optimization_adoption_claimed': False,
               'full_p4_or_scientific_or_runtime_domain_qualification_claimed': False}
    message = 'Preserve two original P4 parity passes and one isolated setup timeout (#53)\n'
elif kind == 'native-provider-integration':
    work = root / 'out/native-provider-integration'
    bundle = root / 'attestations/native-vulkan/provider-integration-5b5fad9'
    tag = 'evidence/issue-48-native-provider-integration-source'
    summary = json.loads((work / 'source-commit.json').read_text())
    runner_source = summary['commit']
    assert runner_source == '5b5fad91c08c0533e29c4f5c7ae249725162e405' and summary['working_checkout_unchanged'] and not summary['nativeexecution_or_qualification_claimed']
    assert not summary['independent_component_review']['material_blockers'] and not summary['independent_workflow_review']['material_blockers']
    checks = summary['checks']
    assert checks['source_only'] and checks['original18_preflight_order_continuity_controls'] and checks['new12_source_boundary_controls'] and checks['actual_candidate_workflow_policy']
    for name, row in checks['candidate_inputs'].items():
        assert pin(work / name) == row['candidate']
        assert subprocess.check_output(['git', 'show', runner_source + ':' + name]) == (work / name).read_bytes()
    mapping = [(work, 'files/out/native-provider-integration')]
    message = 'Preserve source and portable controls for current-source native provider integration (#48)\n'
elif kind == 'native-preflight-order':
    work = root / 'out/native-preflight-order'
    bundle = root / 'attestations/native-vulkan/preflight-order-54c0f3b'
    tag = 'evidence/issue-48-native-preflight-order'
    summary = json.loads((work / 'source-commit.json').read_text())
    runner_source = summary['commit']
    assert runner_source == '54c0f3bf45c5a42e1fb87933a5a82823f45b2ccc' and summary['canonical_checkout_unchanged'] and not summary['native_execution_or_qualification_claimed']
    checks = summary['source_only_checks']
    assert checks['source_only'] and checks['candidate_policy_valid_and_18_invalid_ordering_continuity_controls'] and checks['candidate_producer_original_self_tests'] and checks['integrated_preflight_reuse_viewer_policy'] and checks['existing_integration_boundary_valid_invalid_controls']
    for name, item in checks['candidate_inputs'].items():
        assert pin(work / name) == item['candidate']
        assert subprocess.check_output(['git', 'show', runner_source + ':' + name]) == (work / name).read_bytes()
    mapping = [(work, 'files/out/native-preflight-order')]
    message = 'Preserve source-only native preflight ordering repair controls (#48)\n'
elif kind == 'p4-setup-complete':
    work = root / 'out/p4-parity-setup-diagnostic'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-original-setup-diagnostic'
    tag = 'evidence/issue-53-p4-original-setup-diagnostic'
    owners, results = {}, {}
    for phase in ('compile-device', 'compile-driver', 'link', 'run'):
        result = json.loads((work / (phase + '-result.json')).read_text())
        owner = json.loads((work / phase / 'owner.json').read_text())
        assert result['passed'] and not result['numerical_acceptance'] and not result['qualification_claimed']
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0 and owner['stop_reason'] is None
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes'] and all(row['absent'] for row in owner['observed_births'])
        owners[phase] = owner; results[phase] = result
    assert json.loads((work / 'link-review.json').read_text())['original_vulkan_archive_member_extracted'] is False
    observation = results['run']['setup_observation']
    assert observation['summary']['setup_completed'] and observation['summary']['active_records'] == 0 and not observation['phase_observation']['unclosed_phases']
    mapping = [(work, 'files/out/p4-parity-setup-diagnostic'),
               (root / 'bin/linux-gcc/p4-parity-setup-diagnostic', 'files/bin/linux-gcc/p4-parity-setup-diagnostic')]
    summary = {'schema': 'sirius-original-parity-setup-only-diagnostic-v1', 'results': results,
               'owners': owners, 'source_instrumentation_only': True, 'full_original_spirv_unchanged': True,
               'active_record_dispatches': 0, 'numerical_or_scientific_or_runtime_domain_qualification_claimed': False,
               'original_allocation_contract': 'actual>=requested and<=1MiB; historical extra equality corrected only in new reader'}
    message = 'Preserve complete original full-module parity setup diagnosis (#53)\n'
elif kind == 'p4-scene-incomplete':
    work = root / 'out/science-p4-vulkan-session-scenes'
    bundle = root / 'attestations/scientific-core/4bcd987/p4-vulkan-scene-incomplete'
    tag = 'evidence/issue-53-p4-vulkan-scene-incomplete'
    result = json.loads((work / 'scene-1-result.json').read_text())
    owner = json.loads((work / 'scene-1/owner.json').read_text())
    assert not result['passed'] and not result['p4_complete']
    assert owner['status'] == 'terminal' and owner['stop_reason'] == 'diagnostic_timeout' and owner['child_exit'] == -15
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert all(row['absent'] for row in owner['observed_births'])
    assert not (work / 'scene-2').exists() and not (work / 'scene-1/gtest.xml').exists() and not (work / 'scene-1/ctest.xml').exists()
    text = (work / 'scene-1/stdout.log').read_text()
    assert text.count('0 active rays') == 6 and '[Vulkan] Render evidence:' not in text
    assert not list(work.rglob('*.exr'))
    progress = [line for line in text.splitlines() if '[Vulkan] Retained progress:' in line]
    assert progress
    mapping = [(work, 'files/out/science-p4-vulkan-session-scenes'),
               (root / 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'producer/sirius_render_tests'),
               (root / 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    summary = {'schema': 'sirius-incomplete-original-vulkan-scene-v1', 'result': result, 'owner': owner,
               'original_case': 'VulkanRenderSession.ProceduralVolumetricTurbulenceReachesLiveKernel',
               'software_preparations_completed': 6, 'software_preparation_wall_seconds_from_stdout': 104.158,
               'last_progress': progress[-1], 'complete_frames': 0, 'two_frame_numerical_verdict': False,
               'second_case_launched': False, 'raw_partial_log': 'files/out/science-p4-vulkan-session-scenes/scene-1/ctest.log.tmp715bf',
               'original_assertions_precision_scene_and_provider_unchanged': True,
               'full_p4_or_scientific_or_runtime_domain_qualification_claimed': False}
    message = 'Preserve incomplete original P4 Vulkan full-scene attempt (#53)\n'
elif kind == 'forensic-incomplete':
    work = root / 'out/forensic-alcubierre-png'
    case = work / 'legacy-alcubierre-original-size'
    bundle = root / 'attestations/forensic-images/b739f92/alcubierre-initial-incomplete'
    tag = 'evidence/issue-100-alcubierre-initial-incomplete'
    result = json.loads((case / 'result.json').read_text())
    owner = json.loads((case / 'owner.json').read_text())
    assert result['classification'] == 'INCOMPLETE_DIAGNOSTIC_LIMIT' and not result['complete_original_size_render']
    assert owner['status'] == 'terminal' and owner['stop_reason'] == 'diagnostic_timeout'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
    assert not (work / 'ported-alcubierre-original-size').exists()
    mapping = [(work, 'files/out/forensic-alcubierre-png'),
               (root / 'bin/linux-gcc/forensic-reconstruction-b739f92/scaffold-retry/renderers/bin/Sirius', 'producer/Sirius'),
               (root / 'bin/linux-gcc/forensic-reconstruction-b739f92/scaffold-retry/renderers/bin/sirius', 'producer/sirius')]
    summary = {'schema': 'sirius-incomplete-full-size-reconstruction-v1', 'evidence_kind': 'NEW_RECONSTRUCTION',
               'runner_revision': parent, 'source_export_revision': 'b739f920c2ae67847ed7961e011d2350a839e410',
               'source_and_build_evidence_reused': 'evidence/issue-100-kerr-exr-reconstruction-pair at38c6a264f7ac18d9b75f798337eaa2027e4264ac',
               'result': result, 'owner': owner, 'paired_comparison_performed': False,
               'old_ported1800_action_launched': False, 'historical_tape_identity_claimed': False,
               'qualification_claimed': False,
               'unused_pair_helper': 'Old unused helper has stale EXR512x512 argv; no comparator was executed. The fresh completion packet must use exact PNG1920x1080 before sealing.'}
    message = 'Preserve incomplete original-size Alcubierre reconstruction without a pixel verdict (#100)\n'
else:
    raise SystemExit(kind)

assert not bundle.exists()
bundle.mkdir(parents=True)
if kind == 'native-build-domains':
    summary['archive_parts'] = {}
    for role in ('macos', 'windows'):
        archive = work / (role + '.zip')
        parts = bundle / (role + '.zip.parts')
        parts.mkdir()
        hashes = []
        with archive.open('rb') as stream:
            number = 0
            while block := stream.read(32 * 1024**2):
                path = parts / f'part-{number:03d}'
                path.write_bytes(block)
                hashes.append({'path': str(path.relative_to(bundle)), **pin(path)})
                number += 1
        assert b''.join((bundle / row['path']).read_bytes() for row in hashes) == archive.read_bytes()
        summary['archive_parts'][role] = hashes
if kind in ('native-escape', 'native-side-effects', 'native-metadata-compile', 'native-header-boundary', 'native-four-controls'):
    parts = bundle / 'artifact.zip.parts'
    parts.mkdir()
    hashes = []
    with archive.open('rb') as stream:
        number = 0
        while block := stream.read(32 * 1024**2):
            path = parts / f'part-{number:03d}'
            path.write_bytes(block)
            hashes.append({'path': str(path.relative_to(bundle)), **pin(path)})
            number += 1
    assert b''.join((bundle / row['path']).read_bytes() for row in hashes) == archive.read_bytes()
    summary['archive_parts'] = hashes
for source, relative in mapping:
    destination = bundle / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    if kind == 'sanitizer-debug-o1-focused' and source.is_file() and source.stat().st_size > 32 * 1024**2:
        parts = destination.with_name(destination.name + '.parts')
        parts.mkdir()
        original = pin(source)
        chunks = []
        with source.open('rb') as stream:
            while block := stream.read(32 * 1024**2):
                path = parts / f'part-{len(chunks):03d}'
                path.write_bytes(block)
                chunks.append({'path': str(path.relative_to(bundle)), **pin(path)})
        digest = hashlib.sha256()
        total = 0
        for chunk in chunks:
            block = (bundle / chunk['path']).read_bytes()
            total += len(block)
            digest.update(block)
        assert {'bytes': total, 'sha256': digest.hexdigest()} == original == pin(source)
        summary.setdefault('chunked_inputs', {})[relative] = {'original': original, 'parts': chunks}
        continue
    if source.is_dir():
        shutil.copytree(source, destination, symlinks=True)
    else:
        shutil.copy2(source, destination)
for name in ('scripts/macos_moltenvk_integration.py', 'scripts/macos_vulkan_identity.py',
             'scripts/diagnose-macos-retained.py', 'scripts/macos_metal_probe.py',
             'scripts/runtime_identity.py', 'scripts/governance/build_policy.py', '.github/workflows/ci.yml'):
    destination = bundle / 'runner-source' / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(subprocess.check_output(['git', 'show', runner_source + ':' + name]))
shutil.copy2(Path(__file__), bundle / Path(__file__).name)
write(bundle / 'result.json', summary)
write(bundle / 'manifest.json', {'qualification_claimed': False,
      'files': {str(path.relative_to(bundle)): pin(path) for path in sorted(bundle.rglob('*')) if path.is_file()}})
fd, index = tempfile.mkstemp(dir=root / 'out', prefix='preserve-finite-')
os.close(fd)
os.unlink(index)
env = dict(os.environ, GIT_INDEX_FILE=index)
try:
    subprocess.run(['git', 'read-tree', '--empty'], env=env, check=True)
    rows = []
    paths = [path for path in sorted(bundle.rglob('*')) if path.is_file()]
    for path in paths:
        assert not path.is_symlink()
        blob = subprocess.check_output(['git', 'hash-object', '-w', '--', str(path)], text=True).strip()
        rows.append(f'100644 {blob}\t{path.relative_to(root)}\0')
    subprocess.run(['git', 'update-index', '-z', '--index-info'], input=''.join(rows).encode(), env=env, check=True)
    tree = subprocess.check_output(['git', 'write-tree'], env=env, text=True).strip()
    commit = subprocess.check_output(['git', 'commit-tree', tree, '-p', runner_source], input=message.encode(), env=env).decode().strip()
    subprocess.run(['git', 'tag', tag, commit], check=True)
    for path in paths:
        raw = subprocess.check_output(['git', 'show', commit + ':' + str(path.relative_to(root))])
        assert len(raw) == path.stat().st_size and hashlib.sha256(raw).hexdigest() == pin(path)['sha256']
    print(json.dumps({'bundle': str(bundle.relative_to(root)), 'tag': tag, 'commit': commit,
                      'files': len(paths), 'manifest': pin(bundle / 'manifest.json'), 'all_git_bytes_verified': True}))
finally:
    if Path(index).exists():
        Path(index).unlink()
