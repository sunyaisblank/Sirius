"""ROOT-only data materialization; never launch the frame or a provider.

CLI: python3 -B materialize.py CI_ADMISSION_JSON ARTIFACTS_JSON WINDOWS_BUILD_JSON
All three are checkout-relative, current out/native-constant-current-ci files.
Admission requires passed:true, source_revision, stage, files:{path:{bytes,sha256}},
and true canonical_source_verified/native_gate_verified/all42_PE_payloads_verified.
Its files must include the two actual downloaded ZIPs and selected CI provenance.
Both output admission flags remain false until independent review and ROOT release.
"""
from pathlib import Path, PurePosixPath
import copy
import hashlib
import json
import re
import subprocess
import sys

sys.dont_write_bytecode = True
if sys.flags.optimize:
    raise RuntimeError('Assertions must remain enabled')
WORK = Path(__file__).resolve().parent
ROOT = Path.cwd()
assert ROOT == WORK.parents[1]
SOURCE = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
CI = ROOT / 'out/native-constant-current-ci'
ORIGIN = ROOT / ('attestations/qualification/5b2df58b/native-original-frame-current/'
                 'files/out/native-current-retained-frame-corrected/execution-bindings.json')
HELPERS = ('run_system.py', 'run_terminal.py', 'run_native_control.ps1',
           'terminal-audit.ps1', 'emergency_native_cleanup.ps1',
           'run_owned.py', 'owned_guard.py')


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def identity(path):
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            size += len(block)
            digest.update(block)
    return {'bytes': size, 'sha256': digest.hexdigest()}


def local(name):
    relative = PurePosixPath(name)
    assert not relative.is_absolute() and '..' not in relative.parts
    path = ROOT / relative
    assert path.resolve(strict=True) == path and path.is_file()
    return path


def emit(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def main():
    admission_path, artifacts_path, native_path = map(local, sys.argv[1:])
    assert all(path.is_relative_to(CI) for path in (admission_path, artifacts_path, native_path))
    for name in ('plan.json', 'execution-bindings.json', 'preparation.json',
                 'canonical-source-freeze.json'):
        assert not (WORK / name).exists()
    assert not (WORK / 'native').exists() and not (WORK / 'temp').exists()
    assert not (ROOT / 'renders/constant-original-frame-2c8ef086').exists()
    plan = read(WORK / 'plan.pending.json')
    assert plan['admitted'] is False
    assert plan['source_revision'] == plan['runner_source_revision'] == SOURCE
    assert plan['stage'] == 'bin/windows-msvc/constant-original-frame-2c8ef086'
    assert plan['order'] == ['native-frame'] and plan['actions'] == {
        'native-frame': {'kind': 'Control',
                         'test': 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent',
                         'expected_exit': 0, 'owner_seconds': 300, 'bridge_outer_seconds': 1380}}
    assert (plan['owner_rss_bytes'], plan['minimum_wsl_available_kib'],
            plan['terminal_observer_seconds']) == (4294967296, 2097152, 60)
    assert plan['child_environment'] == {'MESA_SHADER_CACHE_DISABLE': 'true',
        'TEMP': 'renders/constant-original-frame-2c8ef086/native-frame',
        'TMP': 'renders/constant-original-frame-2c8ef086/native-frame'}
    admission = read(admission_path)
    assert admission['passed'] is True and admission['source_revision'] == SOURCE
    assert admission['stage'] == plan['stage']
    for name in ('canonical_source_verified', 'native_gate_verified', 'all42_PE_payloads_verified'):
        assert admission[name] is True
    files = {}

    def add(path, expected=None):
        assert path.is_file()
        key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        current = identity(path)
        if expected is not None:
            assert set(expected) == {'bytes', 'sha256'}
            assert type(expected['bytes']) is int and expected['bytes'] >= 0
            assert re.fullmatch('[0-9a-f]{64}', expected['sha256'])
            assert current == expected, key
        assert key not in files or files[key] == current
        files[key] = current
        return current

    assert type(admission['files']) is dict and admission['files']
    for name, expected in admission['files'].items():
        path = local(name)
        assert path.is_relative_to(CI)
        add(path, expected)
    for path in (admission_path, artifacts_path, native_path, ORIGIN):
        add(path)
    stage = ROOT / plan['stage']
    assert stage.resolve(strict=True) == stage and stage.is_dir()
    manifest_path = stage / 'generated/sirius/windows-diagnostic-build-products.json'
    gate_path = stage / 'generated/sirius/native_build_gate.json'
    manifest, gate, native = read(manifest_path), read(gate_path), read(native_path)
    assert manifest['schema'] == 'sirius-windows-diagnostic-build-products-v1'
    assert manifest['source_revision'] == SOURCE and manifest['repository'] == 'sunyaisblank/Sirius'
    assert manifest['run_id'] == '39116392750' and manifest['run_attempt'] == '1'
    assert manifest['diagnostic_build_only'] is True and manifest['qualification_claimed'] is False
    assert len(manifest['files']) == 138
    expected_stage = {manifest_path, gate_path}
    for name, expected in manifest['files'].items():
        logical = PurePosixPath(name)
        assert name.startswith('bin/windows-msvc/') and '..' not in logical.parts
        assert '\\' not in name and ':' not in name
        path = stage / logical.relative_to('bin/windows-msvc')
        assert path.resolve(strict=True) == path and path.is_file()
        add(path, expected)
        expected_stage.add(path)
    assert gate['kind'] == 'sirius-native-build-gate' and gate['schema_version'] == 2
    assert gate['status'] == 'passed' and gate['source'] == {'revision': SOURCE, 'clean': True}
    assert type(gate['ctest']['executed']) is int and gate['ctest']['executed'] == 9
    assert all(type(gate['ctest'][name]) is int and gate['ctest'][name] == 0
               for name in ('failures', 'errors', 'skipped'))
    for role, relative in (('sirius', 'src/sirius/app/Release/sirius.exe'),
                           ('sirius_backend_tests', 'tests/backend/Release/sirius_backend_tests.exe')):
        item = gate['tested_artifacts'][role]
        assert item['root'] == 'build' and item['path'] == relative
        assert {k: item[k] for k in ('bytes', 'sha256')} == manifest['files']['bin/windows-msvc/' + relative]
    extra_relative = 'src/sirius/app/Release/sirius_render_tests.exe'
    extra = gate['tested_artifacts']['sirius_render_tests']
    assert set(extra) == {'root', 'path', 'bytes', 'sha256'}
    assert extra['root'] == 'build' and extra['path'] == extra_relative
    extra_identity = {k: extra[k] for k in ('bytes', 'sha256')}
    extra_path = stage / extra_relative
    add(extra_path, extra_identity)
    expected_stage.add(extra_path)
    assert len(expected_stage) == 141
    assert {p for p in stage.rglob('*') if p.is_file()} == expected_stage
    assert all(not p.is_symlink() for p in stage.rglob('*'))
    gate_identity, manifest_identity = add(gate_path), add(manifest_path)
    assert native['schema_version'] == 1 and native['status'] == 'pass'
    assert native['domains'] == ['windows-native-build'] and native['source_revision'] == SOURCE
    assert native['platform']['os'] == 'windows' and native['platform']['wsl2'] is False
    assert native['claims']['runtime_ready'] is False
    assert native['claims']['test_estate_kind'] == 'native-build-non-render-authority'
    assert native['claims']['test_estate_passed'] is True
    assert native['claims']['test_report'] == {'cases': 9, 'failures': 0, 'errors': 0, 'skipped': 0}
    for name, expected in (('native_build_gate.json', gate_identity),
                           ('native-build-tested-sirius_render_tests', extra_identity)):
        assert native['artifacts'][name] == dict(path=name, **expected)
        add(native_path.parent / name, expected)
    artifacts = read(artifacts_path)['artifacts']

    def artifact(name):
        item, = [value for value in artifacts if value['name'] == name]
        assert type(item['id']) is int and item['id'] > 0 and item['expired'] is False
        assert str(item['workflow_run']['id']) == manifest['run_id']
        assert item['workflow_run']['head_sha'] == SOURCE
        assert item['digest'].startswith('sha256:')
        assert type(item['size_in_bytes']) is int and item['size_in_bytes'] > 0
        digest = item['digest'].removeprefix('sha256:')
        assert re.fullmatch('[0-9a-f]{64}', digest)
        archive_identity = {'bytes': item['size_in_bytes'], 'sha256': digest}
        assert any(value == archive_identity and key.endswith('.zip') for key, value in files.items())
        return {'repository': manifest['repository'], 'run_id': manifest['run_id'],
                'run_attempt': manifest['run_attempt'], 'name': name, 'id': item['id'],
                'zip_identity': archive_identity}

    diagnostic_artifact = artifact('sirius-windows-diagnostic-build-products')
    native_artifact = artifact('sirius-windows-native-attestation')
    old = read(ORIGIN)
    host = {name: value for name, value in old['files'].items() if name.startswith('/')}
    assert len(host) == 148
    for name, expected in host.items():
        add(Path(name), expected)  # Refuse unexplained host changes; never silently repin.
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
    assert not subprocess.check_output(['git', 'status', '--porcelain'])
    tracked = [name for name in subprocess.check_output(['git', 'ls-files', '-z'], text=True).split('\0') if name]
    assert len(tracked) == 538
    source_rows = {name: add(local(name)) for name in tracked}
    adaptation = read(WORK / 'source-adaptation.json')
    for name in HELPERS:
        add(WORK / name, adaptation['helpers'][name]['adapted'])
    for name in ('plan.pending.json', 'source-adaptation.json', 'materialize.py'):
        add(WORK / name)
    replacements = ((old['source_revision'], SOURCE),
                    ('current-original-frame-5b2df58b', 'constant-original-frame-2c8ef086'),
                    ('native-current-retained-frame-corrected', WORK.name),
                    (old['gate_sha256'], gate_identity['sha256']))

    def retarget(value):
        if isinstance(value, str):
            for previous, current in replacements:
                value = value.replace(previous, current)
            return value
        if isinstance(value, list):
            return [retarget(v) for v in value]
        return {key: retarget(v) for key, v in value.items()}

    binding = {key: copy.deepcopy(plan[key]) for key in
               ('source_revision', 'runner_source_revision', 'control_profile', 'stage', 'windows_checkout_unc')}
    binding.update(schema='sirius-native-windows-current-retained-frame-bindings-v1', admitted=False,
        artifact=diagnostic_artifact, native_artifact=native_artifact,
        manifest_identity=manifest_identity, native_gate_identity=gate_identity,
        gate_sha256=gate_identity['sha256'], native_attestation_path=str(native_path.relative_to(ROOT)),
        native_attestation_identity=identity(native_path),
        additional_test_executable={'role': 'sirius_render_tests', 'stage_relative_path': extra_relative,
                                    'identity': extra_identity},
        expected_amd_module=old['expected_amd_module'], files=files,
        actions=retarget(old['actions']), terminal_actions=retarget(old['terminal_actions']),
        pending_root_fields=['Both admitted flags remain false until independent concrete review and ROOT release.'],
        scope='One original native128x128/1SPP Kerr frame; current canonical source/CI stage141. No private shader override/fullqualification.')
    assert set(binding['actions']) == set(binding['terminal_actions']) == {'native-frame'}
    for name, expected in files.items():
        assert identity(Path(name) if name.startswith('/') else ROOT / name) == expected
    source_freeze = WORK / 'canonical-source-freeze.json'
    emit(source_freeze, {'source_revision': SOURCE, 'files': source_rows})
    add(source_freeze)
    emit(WORK / 'plan.json', plan)
    emit(WORK / 'execution-bindings.json', binding)
    emit(WORK / 'preparation.json', {'schema': 'sirius-native-constant-original-frame-materialization-v1',
        'source_revision': SOURCE, 'materialized': True, 'runtime_admitted': False,
        'CI_admission': {'path': str(admission_path.relative_to(ROOT)), **identity(admission_path)},
        'historical_binding': {'path': str(ORIGIN.relative_to(ROOT)), **identity(ORIGIN)},
        'source_files': len(source_rows), 'matched_historical_external_host_inputs': len(host),
        'stage_files': len(expected_stage), 'bound_files': len(files),
        'materializer': identity(Path(__file__)), 'binding': identity(WORK / 'execution-bindings.json'),
        'plan': identity(WORK / 'plan.json'),
        'failure_policy': 'Firstnonpass stops;positive-only Python terminal reader not run. ROOT uses unchanged bounded PS terminal audit directly for failure teardown. No retry.',
        'activation': 'Held plan.admitted=false and binding.admitted=false;ROOT+independent concrete review before activation. No owner/terminal/temp/native/render directory created.'})
    print(json.dumps({'materialized': True, 'runtime_admitted': False, 'bound_files': len(files)}))


if __name__ == '__main__':
    main()
