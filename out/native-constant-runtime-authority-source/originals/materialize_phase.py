"""Materialize one immutable native diagnostic phase; never launch a workload."""
from pathlib import Path
import ast
import hashlib
import json
import shutil
import subprocess
import sys

ROOT = Path.cwd()
AUTH = ROOT / 'out/native-transport-validation-runtime-authority-source'
READERS = ROOT / 'out/native-transport-validation-runtime-source'


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write(path, value):
    assert not path.exists(), path
    path.write_text(json.dumps(value, indent=2) + '\n')


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def quote(word):
    return "'" + word.replace("'", "''") + "'"


phase, stage_name, artifacts_name, pe_authority_name, *prior_names = sys.argv[1:]
assert phase in ('timestamp', 'abba', 'remaining')
stage = ROOT / stage_name
assert stage.is_relative_to(ROOT / 'bin/windows-msvc') and stage.resolve(strict=True) == stage
artifacts_path, pe_authority_path = ROOT / artifacts_name, ROOT / pe_authority_name
artifacts = read(artifacts_path)['artifacts']
assert len(artifacts) == 1 and artifacts[0]['name'] == 'sirius-native-packed-matched-producer'
artifact = artifacts[0]
producer = read(stage / 'producer-result.json')
manifest = read(stage / 'manifest.json')
assert producer['passed'] is True and manifest['passed'] is True
source = producer['source_revision']
assert source == manifest['source_revision'] == producer['github_sha']
assert artifact['workflow_run']['head_sha'] == source
assert str(artifact['workflow_run']['id']) == producer['run_id']
assert artifact['expired'] is False
pe_authority = read(pe_authority_path)
assert pe_authority['passed'] is True
assert pe_authority['producer_source_revision'] == source
assert pe_authority['actual_PE_payload_and_selector_joins_passed'] is True
assert pe_authority['original27_test_bodies_and1492_assertion_lines_equal'] is True
assert pe_authority['seven_program_tables_and_public_constants_equal'] is True
assert pe_authority['owner2_payloads_equal'] is True
assert pe_authority['manifest_identity'] == identity(stage / 'manifest.json')
assert pe_authority['producer_identity'] == identity(stage / 'producer-result.json')
for arm in ('baseline', 'candidate'):
    relative = 'payload/' + arm + '/original_' + arm + '_controls.exe'
    selected = identity(stage / relative)
    assert pe_authority['arms'][arm]['PE'] == selected
    assert selected == {k: manifest['files'][relative][k] for k in ('bytes', 'sha256')}
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert head == '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'
assert not subprocess.check_output(['git', 'status', '--porcelain'])

# Parse the held source's static tables without importing its execution body.
constants = {}
for node in ast.parse((AUTH / 'run_system.py').read_text()).body:
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
        if node.targets[0].id in ('ACTION_NAMES', 'ACTION_ARMS', 'PROFILES'):
            constants[node.targets[0].id] = ast.literal_eval(node.value)
order = constants['PROFILES'][phase]
work = ROOT / ('out/native-transport-validation-runtime-' + phase)
assert not work.exists()

prior_keys = {'timestamp': [], 'abba': ['timestamp_acceptance'],
              'remaining': ['timestamp_acceptance', 'benefit_acceptance']}[phase]
assert len(prior_names) == len(prior_keys)
prerequisites = {}
prior_closure = {}
for name, relative in zip(prior_keys, prior_names):
    path = ROOT / relative
    proof = read(path)
    assert proof['passed'] is True and proof['source_revision'] == source
    if name == 'timestamp_acceptance':
        assert proof['schema'] == 'sirius-native-packed-original-controls-acceptance-v1'
        assert proof['control8_passed'] is True
    else:
        assert proof['schema'] == 'sirius-native-packed-original22-benefit-v1'
        assert proof['audit_passed'] is True
        assert proof['benefit_passed'] is True and proof['checks_passed'] == proof['checks_total'] == 22
    proof_inputs = {}
    for field in ('reader_bindings', 'root_reader_admission'):
        row = proof[field]
        item = ROOT / row['path']
        assert item.resolve(strict=True) == item and item.is_relative_to(ROOT / 'out')
        expected = {k: row[k] for k in ('bytes', 'sha256')}
        assert identity(item) == expected
        proof_inputs[field] = read(item)
        prior_closure[key(item)] = expected
    context = proof_inputs['reader_bindings']
    admission = proof_inputs['root_reader_admission']
    assert context['schema'] == 'sirius-native-packed-evidence-bindings-v1'
    assert context['admitted'] is True
    assert context['producer_source_revision'] == source and context['runner_source_revision'] == head
    assert context['phase'] == ('timestamp' if name == 'timestamp_acceptance' else 'abba')
    assert admission['schema'] == 'sirius-native-packed-evidence-root-admission-v1'
    assert admission['passed'] is True and admission['reader_bindings'] == proof['reader_bindings']
    assert admission['scope'] == 'native matched retained controls and sole original22 ABBA only'
    for member, row in context['files'].items():
        item = Path(member) if member.startswith('/') else ROOT / member
        assert identity(item) == row, ('Prior raw reader closure changed', member)
        assert member not in prior_closure or prior_closure[member] == row
        prior_closure[member] = row
    prerequisites[name] = {'path': relative, 'identity': identity(path)}

# Rebind current host/tool bytes from the mature accepted input closure. This
# retains its bounded scope; it is not a complete Windows/SDK tool backup.
old_binding = ROOT / 'attestations/qualification/5b2df58b/native-original-frame-current/files/out/native-current-retained-frame-corrected/execution-bindings.json'
old = read(old_binding)
files = dict(prior_closure)
for name, row in old['files'].items():
    if name.startswith('/'):
        path = Path(name)
        assert identity(path) == row, ('Host input changed; assess before materializing', name)
        files[name] = row
    elif not name.startswith(('out/', 'bin/', 'attestations/', 'renders/')):
        path = ROOT / name
        assert identity(path) == row, ('Mature source input changed', name)
        files[name] = row

paths = [ROOT / p for p in subprocess.check_output(['git', 'ls-files', '-z']).decode().split('\0') if p]
for directory in (AUTH, READERS, stage):
    paths.extend(p for p in directory.rglob('*') if p.is_file())
paths.extend([Path(__file__).resolve(), artifacts_path, pe_authority_path, old_binding])
paths.extend(ROOT / item['path'] for item in prerequisites.values())
for path in paths:
    assert path.is_file() and not path.is_symlink(), path
    files[key(path)] = identity(path)

work.mkdir()
for name in ('run_system.py', 'run_terminal.py', 'run_native_control.ps1',
             'terminal-audit.ps1', 'emergency_native_cleanup.ps1'):
    destination = work / name
    shutil.copyfile(AUTH / name, destination)
    assert identity(destination) == identity(AUTH / name)
    files[key(destination)] = identity(destination)

unc = r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'


def win(path):
    assert path.is_relative_to(ROOT)
    return unc + '\\' + str(path.relative_to(ROOT)).replace('/', '\\')


ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
actions, terminal_actions, plan_actions, environments = {}, {}, {}, {}
for label in order:
    arm = constants['ACTION_ARMS'][label]
    test = constants['ACTION_NAMES'][label]
    exe = stage / 'payload' / arm / ('original_' + arm + '_controls.exe')
    output = work / 'native' / label
    temporary = ROOT / 'renders' / work.name / label
    owner_args = ['-Stage', win(exe.parent), '-Output', win(output), '-Filter', test, '-Seconds', '300',
                  '-RssBytes', '4294967296', '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage),
                  '-ExpectedManifestSha256', identity(stage / 'manifest.json')['sha256'],
                  '-ExpectedProducerSha256', identity(stage / 'producer-result.json')['sha256'],
                  '-Arm', arm, '-LiveRevision', head, '-ExecutableName', exe.name, '-InvocationKind', 'Control']
    command = "$ErrorActionPreference='Stop';foreach($v in Get-ChildItem Env:){if($v.Name -match '^(LVP_|SLANG_|MESA_)'){throw 'Unexpected inherited compiler/cache override'}};& " + quote(win(work / 'run_native_control.ps1')) + ' ' + ' '.join(word if i % 2 == 0 else quote(word) for i, word in enumerate(owner_args)) + ';exit $LASTEXITCODE'
    argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-Command', command]
    cleanup = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(work / 'emergency_native_cleanup.ps1'),
               '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(work / (label + '-emergency-cleanup.json'))]
    actions[label] = {'argv': argv, 'cleanup_argv': cleanup}
    prefix = 'terminal-' + label
    temp = work / 'temp' / prefix
    command = "$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath " + quote(win(temp / 'bootstrap.json')) + " -Encoding UTF8;& " + quote(win(work / 'terminal-audit.ps1')) + ' -Workspace ' + quote(win(work)) + ' -RunName ' + quote('native/' + label) + ' -Receipt ' + quote(prefix + '.json') + ' -BuildRoot ' + quote(win(stage))
    argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-Command', command]
    cleanup = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(work / 'emergency_native_cleanup.ps1'),
               '-TemporaryDirectory', win(temp), '-Output', win(temp), '-Receipt', win(work / (prefix + '-emergency-cleanup.json'))]
    terminal_actions[label] = {'argv': argv, 'cleanup_argv': cleanup}
    plan_actions[label] = {'kind': 'Control', 'test': test, 'arm': arm, 'expected_exit': 0, 'owner_seconds': 300, 'bridge_outer_seconds': 1380}
    environments[label] = {'TEMP': key(temporary), 'TMP': key(temporary)}

amd = old['expected_amd_module']
plan = {'schema': 'sirius-native-packed-runtime-plan-v1', 'admitted': True, 'source_revision': source,
        'runner_source_revision': head, 'stage': stage_name, 'windows_checkout_unc': unc,
        'control_profile': phase, 'order': order, 'actions': plan_actions, 'owner_rss_bytes': 4294967296,
        'minimum_wsl_available_kib': 2097152, 'terminal_observer_seconds': 60,
        'resource_semantics': '300s Windows child elapsed/4GiB sampled child WorkingSet; sampled2GiB WSL available reserve; no Windows aggregate/GPU-memory claim',
        'child_environment': environments, 'phase_prerequisites': prerequisites}
binding = {'schema': 'sirius-native-packed-runtime-bindings-v1', 'admitted': True, 'source_revision': source,
           'runner_source_revision': head, 'stage': stage_name, 'windows_checkout_unc': unc, 'control_profile': phase,
           'manifest_identity': identity(stage / 'manifest.json'), 'producer_identity': identity(stage / 'producer-result.json'),
           'artifact': {'repository': producer['repository'], 'name': artifact['name'], 'id': artifact['id'],
                        'run_id': producer['run_id'], 'run_attempt': producer['run_attempt']},
           'selected_executables': {arm: {'relative_path': 'payload/' + arm + '/original_' + arm + '_controls.exe',
                                        'identity': identity(stage / 'payload' / arm / ('original_' + arm + '_controls.exe'))}
                                    for arm in ('baseline', 'candidate')},
           'expected_amd_module': amd, 'private_producer_authority_kind': 'private-ci-producer-result-and-manifest',
           'native_build_gate_authority_used': False, 'files': dict(sorted(files.items())),
           'actions': actions, 'terminal_actions': terminal_actions}
write(work / 'plan.json', plan)
write(work / 'execution-bindings.json', binding)
write(work / 'preparation.json', {'schema': 'sirius-native-packed-runtime-materialization-v1', 'phase': phase,
      'source_revision': source, 'runner_source_revision': head, 'plan': identity(work / 'plan.json'),
      'binding': identity(work / 'execution-bindings.json'), 'materializer': identity(Path(__file__).resolve()),
      'runtime_executed': False, 'scope': 'Concrete immutable phase inputs; ROOT must review and release execution separately.'})
print(json.dumps({'work': key(work), 'phase': phase, 'files': len(files), 'runtime_executed': False}))
