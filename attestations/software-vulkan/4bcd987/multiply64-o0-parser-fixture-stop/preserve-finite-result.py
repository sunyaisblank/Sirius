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
assert parent == ('46a9599a0cb30b0667e09838eac7a87690e982e4' if kind in ('sanitizer-debug-o1-focused', 'multiply64-compiler-screen', 'multiply64-ordinary-jit-screen', 'multiply64-original-arithmetic', 'multiply64-profile-graph-stop', 'multiply64-o0-parser-fixture-stop', 'windows-runner-loss-38577826600', 'p2-cpu-radial-congruence', 'p2-source-sky-differential', 'e3-dispatch-continuations') else 'f5d649bf5067d1d68dc9215462f37830122f0540')
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
