"""Prepare reviewed seven-case bindings; no native invocation or provider query."""
from pathlib import Path, PurePosixPath
import hashlib, json, re, subprocess, sys

sys.dont_write_bytecode = True
assert not sys.flags.optimize and len(sys.argv) == 1
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
SOURCE = 'cd99fd06627cf1f2e1ece74fcf90abf129178940'
assert ROOT == WORK.parents[1] and ROOT.resolve(strict=True) == ROOT
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain'])
def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def record(value):
    assert set(value) == {'bytes', 'sha256'} and type(value['bytes']) is int and value['bytes'] >= 0
    assert type(value['sha256']) is str and re.fullmatch('[0-9a-f]{64}', value['sha256'])
def key(path): return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
def path(name):
    logical = PurePosixPath(name)
    assert type(name) is str and '\\' not in name and ':' not in name and str(logical) == name
    assert all(part not in ('', '.', '..') for part in name.split('/')[1:] if logical.is_absolute())
    assert '..' not in logical.parts and '.' not in logical.parts
    return Path(name) if logical.is_absolute() else ROOT / logical
def checked(location, expected):
    record(expected)
    assert location.is_file() and location.resolve(strict=True) == location and pin(location) == expected
def encoded(document): return (json.dumps(document, indent=2, sort_keys=True) + '\n').encode()
def create(location, data):
    with location.open('xb') as stream: stream.write(data)
for name in ['plan.json', 'execution-bindings.json', 'source-freeze.json', 'native', 'temp']:
    assert not (WORK / name).exists(), name
plan = read(WORK / 'plan.regressions.pending.json')
assert plan['admitted'] is False and plan['source_revision'] == plan['runner_source_revision'] == SOURCE
assert plan['control_profile'] == 'floor-regressions-then-original-kerr'
order = ['floor-fp32', 'floor-fp32comp', 'floor-fp64', 'captured-fp32', 'captured-fp32comp', 'kerr-image', 'kerr-precision']
assert plan['order'] == order and set(plan['actions']) == set(order)
tests = ['KernelParity.RejectedNeighbourRetriesExactMinimumFp32', 'KernelParity.RejectedNeighbourRetriesExactMinimumFp32Comp', 'KernelParity.RejectedNeighbourRetriesExactMinimumFp64', 'KernelParity.CapturedKerrFloorCompositeMatchesCanonicalEndpointFp32', 'KernelParity.CapturedKerrFloorCompositeMatchesCanonicalEndpointFp32Comp', 'KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow', 'KernelTrace.Fp64RungAgreesWithFp32OnKerrScene']
assert [plan['actions'][case]['test'] for case in order] == tests
assert plan['schema'] == 'sirius-native-windows-kerr-composite-plan-v1'
admission_path = WORK / 'root-ci-stage-admission.json'
admission_bytes = admission_path.read_bytes(); admission = json.loads(admission_bytes)
admission_pin = {'bytes': len(admission_bytes), 'sha256': hashlib.sha256(admission_bytes).hexdigest()}
assert admission['passed'] is True and admission['source_revision'] == SOURCE
assert admission['stage'] == plan['stage'] == 'bin/windows-msvc/diagnostic-cd99fd06'
files = dict(admission['files']); assert files and key(admission_path) not in files
for name, expected in files.items(): checked(path(name), expected)
for name in ['plan.json', 'execution-bindings.json', 'source-freeze.json']:
    assert key(WORK / name) not in files
assert not any(path(name).is_relative_to(WORK / leaf) for name in files for leaf in ['native', 'temp'])
stage = ROOT / plan['stage']; manifest_path = stage / 'generated/sirius/windows-diagnostic-build-products.json'
checked(manifest_path, admission['manifest_identity']); manifest = read(manifest_path)
for field, expected in {'schema': 'sirius-windows-diagnostic-build-products-v1', 'source_revision': SOURCE, 'repository': 'sunyaisblank/Sirius', 'run_id': '38735683288', 'run_attempt': '1'}.items():
    assert type(manifest[field]) is str and manifest[field] == expected
assert manifest['diagnostic_build_only'] is True and manifest['qualification_claimed'] is False
artifact = admission['artifact']; assert type(artifact['id']) is int and artifact['id'] > 0
assert artifact['name'] == 'sirius-windows-diagnostic-build-products'; record(artifact['zip_identity'])
assert all(artifact[field] == manifest[field] for field in ['repository', 'run_id', 'run_attempt'])
stage_files = {manifest_path, stage / 'generated/sirius/native_build_gate.json'}
for name, expected in manifest['files'].items():
    logical = PurePosixPath(name)
    assert str(logical) == name and '\\' not in name and ':' not in name and '..' not in logical.parts
    location = stage / logical.relative_to('bin/windows-msvc'); checked(location, expected)
    assert key(location) not in files or files[key(location)] == expected
    files[key(location)] = expected; stage_files.add(location)
assert {location for location in stage.rglob('*') if location.is_file()} == stage_files
assert not any(location.is_symlink() for location in stage.rglob('*'))
files[key(manifest_path)] = admission['manifest_identity']
gate_path = stage / 'generated/sirius/native_build_gate.json'; checked(gate_path, admission['native_gate_identity'])
gate = read(gate_path)
assert type(gate['schema_version']) is int and gate['schema_version'] == 2
assert gate['kind'] == 'sirius-native-build-gate' and gate['status'] == 'passed' and gate['alignment_mode'] == 'qualification'
assert gate['source']['revision'] == SOURCE and gate['source']['clean'] is True
counts = gate['ctest']; assert counts['selected'] == counts['executed'] == 9
assert all(type(counts[name]) is int for name in ['selected', 'executed', 'registered'])
assert all(type(counts[name]) is int and counts[name] == 0 for name in ['failures', 'errors', 'skipped'])
for role, relative in [('sirius', 'src/sirius/app/Release/sirius.exe'), ('sirius_backend_tests', 'tests/backend/Release/sirius_backend_tests.exe')]:
    product = gate['tested_artifacts'][role]
    assert product['root'] == 'build' and product['path'] == relative
    assert {name: product[name] for name in ['bytes', 'sha256']} == manifest['files']['bin/windows-msvc/' + relative]
for relative in ['src/sirius/app/Release/vulkan-1.dll', 'tests/backend/Release/vulkan-1.dll']:
    assert 'bin/windows-msvc/' + relative in manifest['files']
alignment = read(stage / 'src/sirius/app/Release/resources/model/alignment_receipt.json')
assert alignment['source_revision'] == SOURCE and alignment['source_tree_clean'] is True and alignment['alignment_mode'] == 'qualification'
files[key(gate_path)] = admission['native_gate_identity']; files[key(admission_path)] = admission_pin
for name in ['run_system.py', 'run_terminal.py', 'run_native_control.ps1', 'terminal-audit.ps1', 'emergency_native_cleanup.ps1', 'prepare_bindings.py', 'preparation.json', 'source-contracts.json', 'plan.pending.json', 'plan.regressions.pending.json', 'execution-bindings.pending.json']:
    assert key(WORK / name) in files
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'; assert ps in files
amd = admission['expected_amd_module']; amd_pin = {name: amd[name] for name in ['bytes', 'sha256']}; record(amd_pin)
assert 'amdvlk' in amd['name'].lower() and type(amd['path']) is str and amd['path']
assert any(not path(name).is_relative_to(ROOT) and path(name).name.casefold() == amd['name'].casefold() and expected == amd_pin for name, expected in files.items())
def win(location):
    location = location.resolve(); assert location.is_relative_to(ROOT)
    return plan['windows_checkout_unc'] + '\\' + str(location.relative_to(ROOT)).replace('/', '\\')
def quote(word): return "'" + word.replace("'", "''") + "'"
binding = read(WORK / 'execution-bindings.pending.json')
assert binding['admitted'] is False and binding['schema'] == 'sirius-native-windows-kerr-composite-bindings-v1'
assert binding['source_revision'] == binding['runner_source_revision'] == SOURCE and binding['stage'] == plan['stage']
assert binding['windows_checkout_unc'] == plan['windows_checkout_unc'] == r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'
binding.update(admitted=True, control_profile=plan['control_profile'], artifact=artifact, manifest_identity=admission['manifest_identity'], native_gate_identity=admission['native_gate_identity'], gate_sha256=admission['native_gate_identity']['sha256'], expected_amd_module=amd, actions={}, terminal_actions={}, pending_root_fields=[])
for case in order:
    action = plan['actions'][case]; seconds = 28800 if case.startswith('kerr-') else 900
    assert set(action) == {'kind', 'test', 'expected_exit', 'owner_seconds', 'bridge_outer_seconds'}
    assert action['kind'] == 'Control' and type(action['expected_exit']) is int and action['expected_exit'] == 0
    assert type(action['owner_seconds']) is int and action['owner_seconds'] == seconds and type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds'] == seconds + 480
    output = WORK / 'native' / case; temporary = WORK / 'temp' / case; exe = stage / 'tests/backend/Release/sirius_backend_tests.exe'
    argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'), '-Stage', win(exe.parent), '-Output', win(output), '-Filter', action['test'], '-Seconds', str(seconds), '-RssBytes', str(plan['owner_rss_bytes']), '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage), '-ExpectedGateSha256', binding['gate_sha256'], '-LiveRevision', SOURCE, '-ExecutableName', exe.name, '-InvocationKind', 'Control']
    cleanup = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'), '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(WORK / (case + '-emergency-cleanup.json'))]
    binding['actions'][case] = {'argv': argv, 'cleanup_argv': cleanup}
    prefix = 'terminal-' + case; temp = WORK / 'temp' / prefix; receipt = WORK / (prefix + '.json')
    command = "$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath " + quote(win(temp / 'bootstrap.json')) + " -Encoding UTF8;& " + quote(win(WORK / 'terminal-audit.ps1')) + " -Workspace " + quote(win(WORK)) + " -RunName " + quote('native/' + case) + " -Receipt " + quote(receipt.name) + ' -BuildRoot ' + quote(win(stage))
    argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-Command', command]
    cleanup = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'), '-TemporaryDirectory', win(temp), '-Output', win(temp), '-Receipt', win(WORK / (prefix + '-emergency-cleanup.json'))]
    binding['terminal_actions'][case] = {'argv': argv, 'cleanup_argv': cleanup}
assert plan['owner_rss_bytes'] == 4294967296 and plan['minimum_wsl_available_kib'] == 2097152 and plan['terminal_observer_seconds'] == 60
binding['scope'] = 'Seven explicitly selected original controls; current native producer authority only. No fullMandatory, horizon-domain, performance or release qualification.'
plan['admitted'] = True; create(WORK / 'plan.json', encoded(plan)); files[key(WORK / 'plan.json')] = pin(WORK / 'plan.json')
create(WORK / 'source-freeze.json', encoded({'schema': 'sirius-native-windows-kerr-composite-source-freeze-v1', 'source_revision': SOURCE, 'files': files}))
files[key(WORK / 'source-freeze.json')] = pin(WORK / 'source-freeze.json')
for name, expected in files.items(): checked(path(name), expected)
binding['files'] = files; create(WORK / 'execution-bindings.json', encoded(binding))
print(json.dumps({'source_revision': SOURCE, 'control_profile': plan['control_profile'], 'input_count': len(files), 'plan': pin(WORK / 'plan.json'), 'bindings': pin(WORK / 'execution-bindings.json'), 'freeze': pin(WORK / 'source-freeze.json')}))
