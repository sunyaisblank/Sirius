"""Seal completed current records for existing readers; never launch anything.

ROOT-only invocation after terminal acceptance:
  python3 -B materialize_reader_context.py timestamp|abba|all8 ARRAY_AUTHORITY
ARRAY_AUTHORITY is the current bound reader-form JSON with lowercase
fma_available keys transport/endpoint. All producer identities come from the
current phase binding/manifest; no future producer revision is guessed here.
"""
from pathlib import Path
import copy
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'out/native-constant-runtime-source'


def identity(path):
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    return {'bytes': size, 'sha256': digest.hexdigest()}


def document(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def pin(row):
    return {k: row[k] for k in ('bytes', 'sha256')}


def local(name):
    path = ROOT / name
    assert not Path(name).is_absolute() and path.resolve(strict=True) == path
    assert path.is_relative_to(ROOT / 'out') and path.is_file()
    return path


def main():
    phase, authority_name = sys.argv[1:]
    assert phase in ('timestamp', 'abba', 'all8')
    profile = 'remaining' if phase == 'all8' else phase
    work = ROOT / ('out/native-constant-runtime-' + profile)
    output = ROOT / ('out/native-constant-' +
                     {'timestamp': 'timestamp-evidence', 'abba': 'ABBA-evidence',
                      'all8': 'all8-evidence'}[phase])
    assert not output.exists()
    bindings = document(work / 'execution-bindings.json')
    plan = document(work / 'plan.json')
    assert bindings['schema'] == 'sirius-native-packed-runtime-bindings-v1'
    assert plan['schema'] == 'sirius-native-packed-runtime-plan-v1'
    assert bindings['admitted'] is True and plan['admitted'] is True
    producer, runner = bindings['source_revision'], bindings['runner_source_revision']
    assert re.fullmatch('[0-9a-f]{40}', producer) and re.fullmatch('[0-9a-f]{40}', runner)
    assert plan['source_revision'] == producer and plan['runner_source_revision'] == runner
    assert plan['control_profile'] == bindings['control_profile'] == profile
    stage = ROOT / bindings['stage']
    manifest, produced = (document(stage / name) for name in ('manifest.json', 'producer-result.json'))
    assert manifest['source_revision'] == produced['source_revision'] == producer
    assert identity(stage / 'manifest.json') == bindings['manifest_identity']
    assert identity(stage / 'producer-result.json') == bindings['producer_identity']
    assert manifest['passed'] is True and produced['passed'] is True
    labels = {'timestamp': ['control-8'], 'abba': ['a1', 'b1', 'b2', 'a2'],
              'all8': ['control-' + str(i) for i in range(1, 8)]}[phase]
    assert plan['order'] == labels
    files, cases = {}, {}

    def add(path, expected=None):
        path = Path(path)
        assert path.is_file()
        key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        current = identity(path)
        if expected is not None:
            assert current == pin(expected), key
        assert key not in files or files[key] == current
        files[key] = current

    for key, expected in bindings['files'].items():
        add(Path(key) if key.startswith('/') else ROOT / key, expected)
    for directory in (work, ROOT / 'renders' / work.name):
        for path in sorted(directory.rglob('*')):
            if path.is_file():
                assert not path.is_symlink()
                add(path)
    for name in ('native_evidence.py', 'numeric_gates.py', 'compare_native.py',
                 'read_controls_native.py', 'expected-readback-layout.json'):
        add(SOURCE / name)
    add(Path(__file__).resolve())
    authority = local(authority_name)
    add(authority, bindings['files'][authority_name])
    av = document(authority)
    assert av['schema'] == 'sirius-native-packed-matched-array-authority-v1'
    assert av['passed'] is True and av['producer_source_revision'] == producer
    assert av['manifest_identity'] == bindings['manifest_identity']
    assert av['producer_identity'] == bindings['producer_identity']
    for field in ('actual_PE_payload_and_selector_joins_passed',
                  'original27_test_bodies_and1492_assertion_lines_equal',
                  'seven_program_tables_and_public_constants_equal', 'owner2_payloads_equal'):
        assert av[field] is True
    assert set(av['arms']) == {'baseline', 'candidate'}
    left, right = (av['arms'][arm]['arrays'] for arm in ('baseline', 'candidate'))
    assert len(left) == len(right) == 40 and set(left) == set(right)
    assert {name for name in left if pin(left[name]) != pin(right[name])} == {'kTransportShader', 'kEndpointShader', 'kDenseShader', 'kInitializeShader', 'kDopriPhaseShader'}
    for arm in ('baseline', 'candidate'):
        relative = 'payload/' + arm + '/original_' + arm + '_controls.exe'
        assert av['arms'][arm]['PE'] == pin(manifest['files'][relative]) == identity(stage / relative)
    for arm in av['arms'].values():
        assert set(arm['fma_available']) == {'transport', 'endpoint'}
        assert all(type(value) is bool for value in arm['fma_available'].values())
    assert av['arms']['baseline']['fma_available'] == av['arms']['candidate']['fma_available']
    contexts = {}
    required = {'timestamp': [], 'abba': ['timestamp_acceptance'],
                'all8': ['timestamp_acceptance', 'benefit_acceptance']}[phase]
    assert set(plan['phase_prerequisites']) == set(required)
    for name in required:
        entry = plan['phase_prerequisites'][name]
        proof_path = local(entry['path'])
        add(proof_path, entry['identity'])
        proof = document(proof_path)
        assert proof['passed'] is True and proof['source_revision'] == producer
        if name == 'timestamp_acceptance':
            assert proof['schema'] == 'sirius-native-packed-original-controls-acceptance-v1'
            assert proof['control8_passed'] is True and proof['original_controls_passed'] == 1
        else:
            assert proof['schema'] == 'sirius-native-packed-original22-benefit-v1'
            assert proof['audit_passed'] is True and proof['benefit_passed'] is True
            assert proof['checks_passed'] == proof['checks_total'] == 22
        context_path, admission_path = (local(proof[field]['path'])
                                        for field in ('reader_bindings', 'root_reader_admission'))
        add(context_path, proof['reader_bindings'])
        add(admission_path, proof['root_reader_admission'])
        context, admission = document(context_path), document(admission_path)
        assert context['schema'] == 'sirius-native-packed-evidence-bindings-v1'
        assert context['admitted'] is True and context['producer_source_revision'] == producer
        assert context['runner_source_revision'] == runner
        assert context['phase'] == ('timestamp' if name == 'timestamp_acceptance' else 'abba')
        expected_cases = {'control-8'} if name == 'timestamp_acceptance' else {'control-8', 'a1', 'b1', 'b2', 'a2'}
        assert set(context['cases']) == expected_cases
        assert admission['schema'] == 'sirius-native-packed-evidence-root-admission-v1'
        assert admission['passed'] is True and admission['reader_bindings'] == proof['reader_bindings']
        assert admission['scope'] == 'native matched retained controls and sole original22 ABBA only'
        for key, expected in context['files'].items():
            assert pin(bindings['files'][key]) == pin(expected), ('prior closure not frozen', key)
            add(Path(key) if key.startswith('/') else ROOT / key, expected)
        contexts[name] = context
    if required:
        cases['control-8'] = copy.deepcopy(contexts['timestamp_acceptance']['cases']['control-8'])
    for label in labels:
        accepted = work / ('terminal-' + label + '-accepted.json')
        item = document(accepted)
        assert item['pass'] is True and item['case'] == item['label'] == label
        cases[label] = {'work': str(work.relative_to(ROOT)),
                        'arm': 'baseline' if label in ('a1', 'a2') else 'candidate',
                        'profile': profile, 'terminal_accepted': identity(accepted)}
    assert set(cases) == ({'control-8'} if phase == 'timestamp' else
                         {'control-8', 'a1', 'b1', 'b2', 'a2'} if phase == 'abba' else
                         {'control-' + str(i) for i in range(1, 9)})
    loaders, hardware = {}, None
    for label, row in cases.items():
        current_work = local(row['work'] + '/execution-bindings.json').parent
        current_binding = document(current_work / 'execution-bindings.json')
        assert current_binding['source_revision'] == producer and current_binding['runner_source_revision'] == runner
        accepted_path = current_work / ('terminal-' + label + '-accepted.json')
        add(accepted_path, row['terminal_accepted'])
        accepted = document(accepted_path)
        assert accepted['pass'] is True and accepted['case'] == accepted['label'] == label
        for member, expected in accepted['record']['whole_inputs'].items():
            assert not Path(member).is_absolute() and '..' not in Path(member).parts
            add(current_work / member, expected)
        modules = document(current_work / 'native' / label / 'driver-modules.json')
        loader, = [m for m in modules if m['name'].casefold() == 'vulkan-1.dll']
        amd, = [m for m in modules if 'amdvlk' in m['name'].casefold()]
        assert amd == current_binding['expected_amd_module']
        arm = row['arm']
        assert arm not in loaders or loaders[arm] == loader
        expected_sdk = 'payload/' + arm + '/vulkan-1.dll'
        assert pin(loader) == pin(manifest['files'][expected_sdk])
        expected_path = bindings['windows_checkout_unc'] + '\\' + (bindings['stage'] + '/' + expected_sdk).replace('/', '\\')
        assert loader['path'].replace('/', '\\').casefold() == expected_path.casefold()
        loaders[arm] = loader
        if label == 'control-8' or label in ('a1', 'b1', 'b2', 'a2'):
            props = {}
            for node in ET.parse(current_work / 'native' / label / 'gtest.xml').findall('.//property'):
                key = node.get('name')
                assert key not in props
                props[key] = node.get('value')
            observed = {k: props[k] for k in ('selected_device', 'selected_driver',
                        'retained_scalar_route', 'fma_fp32_enabled',
                        'preserves_fp32_signed_zero_inf_nan', 'timestamp_valid_bits', 'timestamp_period_ns')}
            assert hardware is None or hardware == observed
            hardware = observed
    assert hardware is not None
    if phase == 'timestamp':
        baseline = copy.deepcopy(loaders['candidate'])
        baseline['path'] = bindings['windows_checkout_unc'] + '\\' + (bindings['stage'] + '/payload/baseline/vulkan-1.dll').replace('/', '\\')
        assert identity(stage / 'payload/baseline/vulkan-1.dll') == identity(stage / 'payload/candidate/vulkan-1.dll')
        loaders['baseline'] = baseline
    elif phase == 'all8':
        prior_loaders = contexts['benefit_acceptance']['loader_modules']
        assert prior_loaders['candidate'] == loaders['candidate']
        loaders['baseline'] = copy.deepcopy(prior_loaders['baseline'])
    assert set(loaders) == {'baseline', 'candidate'}
    hardware.update(supports_fp64=True, rounds_fp64_to_nearest=True)
    value = {'schema': 'sirius-native-packed-evidence-bindings-v1', 'admitted': True,
             'phase': phase, 'producer_source_revision': producer, 'runner_source_revision': runner,
             'observed_instance_api_minimum_improvement_percent': 0, 'files': files,
             'array_authority': {'path': authority_name, **identity(authority)},
             'native_identity': hardware,
             'native_identity_scope': 'Current accepted control8/ABBA XML observations. FP64/RTE64 true are expected same-provider capability requirements,not observed wide execution;the unchanged all8 Factored reader must prove its actual wide route.',
             'loader_modules': loaders,
             'loader_module_scope': ('Candidate observed;baseline path/version expected for identical stagedSDKbytes and still needs actualABBA observation.' if phase == 'timestamp' else 'Both arms observed in current acceptedABBA children,or its still-current sealed predecessor for all8.'),
             'cases': cases}
    if phase == 'all8':
        entry = plan['phase_prerequisites']['benefit_acceptance']
        value['benefit_acceptance'] = {'path': entry['path'], **pin(entry['identity'])}
    # Complete all verification before creating the exclusive output directory.
    for key, expected in files.items():
        assert identity(Path(key) if key.startswith('/') else ROOT / key) == expected
    output.mkdir()
    context = output / 'reader-bindings.json'
    context.write_text(json.dumps(value, indent=2) + '\n')
    admission = {'schema': 'sirius-native-packed-evidence-root-admission-v1', 'passed': True,
                 'reader_bindings': {'path': str(context.relative_to(ROOT)), **identity(context)},
                 'scope': 'native matched retained controls and sole original22 ABBA only'}
    (output / 'root-admission.json').write_text(json.dumps(admission, indent=2) + '\n')
    print(json.dumps({'phase': phase, 'files': len(files), 'output': str(output.relative_to(ROOT)),
                      'workload_executed': False}))


if __name__ == '__main__':
    assert not sys.flags.optimize
    main()
