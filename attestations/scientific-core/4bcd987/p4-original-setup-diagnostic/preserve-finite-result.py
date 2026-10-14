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
assert parent == 'f5d649bf5067d1d68dc9215462f37830122f0540'
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
elif kind in ('e1', 'p1-e4'):
    import xml.etree.ElementTree as ET
    if kind == 'e1':
        label, work, case = ('e1-field-equations', root / 'out/science-e1-field-equations', 'original-two-field-equation-controls')
        expected = {'KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell',
                    'KerrSchildFieldEquations.KerrNewmanPotentialIsSourceFreeOutsideTheRing'}
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
