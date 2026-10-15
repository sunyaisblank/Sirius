"""Capture one expected original Kerr failure using unchanged native PowerShell ownership."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import math
import struct
import os
import re
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'): raise RuntimeError('Python assertions must be enabled')
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
if ROOT != WORK.parents[1]: raise RuntimeError('Run from the canonical checkout')
ACTION_NAMES = {'kerr-image':'KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow'}
ORDER = list(ACTION_NAMES)
OWNER_SECONDS = {'kerr-image':28800}
kind, = sys.argv[1:]
assert kind in ACTION_NAMES
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['runner_source_revision']
SOURCE = plan['source_revision']
assert plan['admitted'] is True and binding['admitted'] is True
assert SOURCE == HEAD == '02ee51f85e45e8342923bae1d9f26d50f1d39e6d'
assert plan['stage'] == binding['stage'] == 'bin/windows-msvc/diagnostic-02ee51f'
assert plan['windows_checkout_unc']==binding['windows_checkout_unc']==r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'
signals = []
predecessors = {}
for number in (signal.SIGINT, signal.SIGTERM):
    signal.signal(number, lambda number, frame: signals.append(number))

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, d):
    temp = p.with_suffix(p.suffix + '.partial')
    temp.write_text(json.dumps(d, indent=2) + '\n')
    temp.replace(p)

def clean():
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
    assert not subprocess.check_output(['git', 'status', '--porcelain'])

def read_failure_context(xml_path):
    xml = ET.parse(xml_path).getroot()
    suites = xml.findall('.//testsuite')
    cases = xml.findall('.//testcase')
    assert len(suites) == len(cases) == 1
    for node in [xml, *suites]:
        assert int(node.get('tests', '-1')) == 1 and int(node.get('failures', '-1')) == 1
        assert all(int(node.get(key, '0')) == 0 for key in ['errors', 'skipped', 'disabled'])
    case = cases[0]
    assert case.get('classname')+'.'+case.get('name') == ACTION_NAMES['kerr-image']
    assert case.get('status') == 'run' and case.get('result') == 'completed'
    failures = case.findall('failure')
    assert len(failures) == 1 and len(xml.findall('.//failure')) == 1
    assert not any(xml.findall('.//'+key) for key in ['error', 'skipped'])
    message = failures[0].get('message')
    assert isinstance(message, str) and 'output.has_value()' in message
    assert 'kernel: validate trace continuation: trace termination 5 cannot publish radiance' in message
    marker = '{"schema":"sirius-trace-continuation-failure-v1",'
    assert message.count(marker) == 1
    start = message.index(marker)
    context, end = json.JSONDecoder().raw_decode(message[start:])
    raw_context = message[start:start+end]
    assert set(context) == {'schema','scalar_bytes','packet_scalar_bytes','record_bytes','byte_order',
                            'action','row','active','active_rectangle','expected_identity',
                            'physical_submissions','maximum_attempts','packet_hex',
                            'previous_record_hex','current_record_hex','device_name_bytes_hex'}
    assert context['schema'] == 'sirius-trace-continuation-failure-v1'
    for key in ['scalar_bytes','packet_scalar_bytes','record_bytes','action','row','active',
                'physical_submissions','maximum_attempts']:
        assert type(context[key]) is int
    assert context['scalar_bytes'] == context['packet_scalar_bytes'] == 4
    assert context['record_bytes'] == 272 and context['byte_order'] == 'little'
    assert context['action'] == 1 and context['active'] == 64 and 0 <= context['row'] < 64
    assert context['physical_submissions'] > 0 and context['maximum_attempts'] == 2000
    def raw(key, size):
        text = context[key]
        assert isinstance(text, str) and re.fullmatch('[0-9a-f]{'+str(size*2)+'}', text)
        return bytes.fromhex(text)
    packet_bytes = raw('packet_hex', 70*4)
    previous = raw('previous_record_hex', 272)
    current = raw('current_record_hex', 272)
    packet = struct.unpack('<70f', packet_bytes)
    words = struct.unpack('<70I', packet_bytes)
    rectangle = context['active_rectangle']
    identity = context['expected_identity']
    assert type(rectangle) is list and len(rectangle) == 4 and all(type(x) is int for x in rectangle)
    assert type(identity) is list and len(identity) == 4 and all(type(x) is int for x in identity)
    x, y, width, height = rectangle
    assert x == 0 and 0 <= y < 64 and width == 64 and height == 1
    assert tuple(packet[31:35]) == tuple(rectangle) and packet[0] == packet[1] == 64.0
    assert packet[21] == 2000.0 and packet[46] == 0.0 and packet[68] == 1.0
    assert identity == [x+context['row'] % width, y+context['row'] // width, 0, words[69]] and words[69] > 0
    previous_control = struct.unpack_from('<4I', previous, 240)
    current_control = struct.unpack_from('<4I', current, 240)
    assert list(struct.unpack_from('<4I', previous, 256)) == identity
    assert list(struct.unpack_from('<4I', current, 256)) == identity
    assert current_control[0] == 2 and current_control[1] == 5
    assert current_control[3] <= current_control[2] <= context['maximum_attempts']
    assert all(math.isfinite(x) for x in struct.unpack_from('<60f', current))
    name_hex = context['device_name_bytes_hex']
    assert isinstance(name_hex, str) and re.fullmatch('(?:[0-9a-f]{2})+', name_hex)
    assert bytes.fromhex(name_hex).decode('utf-8') == 'AMD Radeon 780M Graphics'
    rollback = {
        'physical_vectors_0_through_11_equal': previous[:192] == current[:192],
        'accepted_affine_length_last_interval_beam_seed_equal': previous[196:208] == current[196:208],
        'source_and_volume_equal': previous[208:240] == current[208:240],
        'accepted_attempt_count_equal': previous_control[3] == current_control[3],
        'one_attempt_added': current_control[2] == previous_control[2]+1,
        'previous_running': previous_control[0] == 1 and previous_control[1] == 0,
    }
    return {'schema':'sirius-native-windows-kerr-failure-context-v1',
            'raw_context_json':raw_context, 'previous_control':list(previous_control),
            'current_control':list(current_control), 'rollback_comparisons':rollback,
            'rollback_comparison_is_observation_not_acceptance':True,
            'numerical_test_passed':False, 'radiance_publication_admitted':False}

def verify_manifest():
    stage = ROOT / plan['stage']
    manifest_path = stage / 'generated/sirius/windows-diagnostic-build-products.json'
    assert {'bytes': manifest_path.stat().st_size, 'sha256': digest(manifest_path)} == binding['manifest_identity']
    assert binding['files'][str(manifest_path.relative_to(ROOT))] == binding['manifest_identity']
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    for key, value in {'schema':'sirius-windows-diagnostic-build-products-v1', 'source_revision':SOURCE,
                       'repository':'sunyaisblank/Sirius', 'run_id':'38711720394', 'run_attempt':'1'}.items():
        assert type(manifest[key]) is str and manifest[key] == value
    assert manifest['diagnostic_build_only'] is True and manifest['qualification_claimed'] is False
    assert binding['artifact']['repository'] == manifest['repository'] and binding['artifact']['run_id'] == manifest['run_id']
    assert binding['artifact']['run_attempt'] == manifest['run_attempt'] and binding['artifact']['name'] == 'sirius-windows-diagnostic-build-products'
    assert type(binding['artifact']['id']) is int and binding['artifact']['id'] > 0
    assert isinstance(manifest['files'], dict) and manifest['files']
    for key, row in manifest['files'].items():
        logical = PurePosixPath(key)
        assert isinstance(key, str) and key.startswith('bin/windows-msvc/') and '\\' not in key and ':' not in key
        assert not logical.is_absolute() and str(logical) == key and all(part not in ('', '.', '..') for part in key.split('/'))
        assert set(row) == {'bytes', 'sha256'} and type(row['bytes']) is int and row['bytes'] >= 0 and re.fullmatch('[0-9a-f]{64}', row['sha256'])
        relative = logical.relative_to('bin/windows-msvc')
        path = stage / relative
        assert path.is_file() and path.resolve(strict=True) == path and path.is_relative_to(stage)
        assert binding['files'][str(path.relative_to(ROOT))] == row
    for relative in ['src/sirius/app/Release/sirius.exe', 'tests/backend/Release/sirius_backend_tests.exe',
                     'src/sirius/app/Release/vulkan-1.dll', 'tests/backend/Release/vulkan-1.dll']:
        assert 'bin/windows-msvc/' + relative in manifest['files']
    alignment_path = stage / 'src/sirius/app/Release/resources/model/alignment_receipt.json'
    alignment = json.loads(alignment_path.read_text(encoding='utf-8-sig'))
    assert alignment['source_revision'] == SOURCE and alignment['source_tree_clean'] is True and alignment['alignment_mode'] == 'qualification'
    assert binding['files'][str(alignment_path.relative_to(ROOT))] == manifest['files']['bin/windows-msvc/src/sirius/app/Release/resources/model/alignment_receipt.json']
    assert type(binding['manifest_identity']['bytes']) is int and binding['manifest_identity']['bytes'] > 0
    assert re.fullmatch('[0-9a-f]{64}', binding['manifest_identity']['sha256'])
    assert binding['manifest_identity'] == {'bytes':16928,'sha256':'03fe9af83a424b003965137d48cbfd0aa8276739836f4272298f87ceeff05864'}
    assert binding['artifact']['id'] == 11961403981
    assert len(manifest['files']) == 74
    assert binding['ci_native_gate_failed'] is True and binding['ci_native_gate_accepted'] is False
    assert binding['diagnostic_build_manifest_authority_used'] is True and binding['native_build_gate_authority_used'] is False
    return manifest

def verify():
    assert not signals
    clean()
    for rel, item in binding['files'].items():
        p = Path(rel) if rel.startswith('/') else ROOT / rel
        assert {'bytes': p.stat().st_size, 'sha256': digest(p)} == item, rel
    for rel, item in predecessors.items():
        p = ROOT / rel
        assert {'bytes': p.stat().st_size, 'sha256': digest(p)} == item, rel
    assert digest(binding_path) == binding_hash
    verify_manifest()

def win(p):
    p = p.resolve()
    assert p.is_relative_to(ROOT)
    return plan['windows_checkout_unc'] + '\\' + str(p.relative_to(ROOT)).replace('/', '\\')

def available_kib():
    rows = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(rows['MemAvailable'].split()[0])

def input_records():
    paths = [ROOT / rel for rel in binding['files']] + [binding_path, WORK / 'plan.json', WORK / 'preparation.json'] + [ROOT / rel for rel in predecessors]
    return {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in paths}

def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except FileNotFoundError:
        return None

binding_hash = digest(binding_path)
assert binding['source_revision'] == SOURCE and binding['runner_source_revision'] == HEAD
assert plan['order'] == ORDER and set(plan['actions']) == set(ORDER)
for label in ORDER:
    action = plan['actions'][label]
    assert set(action) == {'kind','test','expected_exit','owner_seconds','bridge_outer_seconds'}
    assert action['kind'] == 'Control' and action['test'] == ACTION_NAMES[label]
    assert type(action['expected_exit']) is int and action['expected_exit'] == 1
    assert type(action['owner_seconds']) is int and action['owner_seconds'] == OWNER_SECONDS[label]
    assert type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds'] == OWNER_SECONDS[label]+480
assert plan['owner_rss_bytes'] == 4294967296 and plan['terminal_observer_seconds'] == 60
bridge_seconds = plan['actions'][kind]['bridge_outer_seconds']
assert plan['minimum_wsl_available_kib'] == 2097152

verify()
stage = ROOT / plan['stage']
output = WORK / 'native' / kind
temporary = WORK / 'temp' / kind
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
assert not (WORK / (kind+'-inputs-before.json')).exists() and not (WORK / (kind+'-inputs-after.json')).exists()
exe = stage / 'tests/backend/Release/sirius_backend_tests.exe'
case = ACTION_NAMES[kind]
seconds = plan['actions'][kind]['owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'),
        '-Stage', win(exe.parent), '-Output', win(output), '-Filter', case, '-Seconds', str(seconds),
        '-RssBytes', str(plan['owner_rss_bytes']), '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage),
        '-ExpectedManifestSha256', binding['manifest_identity']['sha256'], '-LiveRevision', HEAD, '-ExecutableName', exe.name,
        '-InvocationKind', 'Control']
cleanup_argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'),
                '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(WORK / (kind + '-emergency-cleanup.json'))]
assert argv == binding['actions'][kind]['argv'] and cleanup_argv == binding['actions'][kind]['cleanup_argv']
before = input_records()
write(WORK / (kind+'-inputs-before.json'), before)
record = {'source_revision': SOURCE, 'runner_source_revision': HEAD, 'kind': kind, 'argv': argv, 'binding_sha256': binding_hash,
          'linux_controller_pid': os.getpid(), 'linux_controller_birth': birth(os.getpid()),
          'predecessor_inputs': predecessors,
          'bridge_outer_seconds': bridge_seconds, 'windows_owner_seconds': seconds, 'minimum_wsl_available_kib': plan['minimum_wsl_available_kib'],
          'resource_scope': plan['resource_semantics'],
          'accepted_owner_execution': False, 'numerical_test_passed': False, 'expected_child_exit': 1, 'completed': False, 'cleanup_errors': [], 'full_qualification_claimed': False}
write(record_path, record)
started = time.monotonic()
child = None
try:
    assert not signals
    initial_available = available_kib()
    record['minimum_sampled_wsl_available_kib'] = initial_available
    assert initial_available >= plan['minimum_wsl_available_kib'], 'WSL sampled available-memory reserve below 2GiB'
    with (WORK / (kind + '-bridge.stdout')).open('wb') as out, (WORK / (kind + '-bridge.stderr')).open('wb') as err:
        child = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
        record['linux_child_pid'] = child.pid
        record['linux_child_birth'] = birth(child.pid)
        assert record['linux_child_birth'] is not None
        write(record_path, record)
        while child.poll() is None:
            assert not signals
            available = available_kib()
            record['minimum_sampled_wsl_available_kib'] = min(record.get('minimum_sampled_wsl_available_kib', available), available)
            assert available >= plan['minimum_wsl_available_kib'], 'WSL sampled available-memory reserve below 2GiB'
            assert time.monotonic() - started < bridge_seconds, 'bridge outer diagnostic stop'
            try:
                child.wait(timeout=0.25)
            except subprocess.TimeoutExpired:
                pass
        record['returncode'] = child.wait(timeout=5)
except BaseException as error:
    record['error'] = repr(error)
finally:
    if child is not None and (record.get('error') or signals or child.poll() is None):
        try:
            with (WORK / (kind + '-cleanup.stdout')).open('wb') as out, (WORK / (kind + '-cleanup.stderr')).open('wb') as err:
                result = subprocess.run(cleanup_argv, stdout=out, stderr=err, timeout=60)
            assert result.returncode == 0
            receipt = json.loads((WORK / (kind + '-emergency-cleanup.json')).read_text(encoding='utf-8-sig'))
            assert receipt['recorded_births_absent'] and not receipt['errors']
        except BaseException as error:
            record['cleanup_errors'].append(repr(error))
    if child is not None:
        try:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except BaseException as error:
            record['cleanup_errors'].append('signal: ' + repr(error))
        try:
            record['returncode'] = child.wait(timeout=15)
        except BaseException as error:
            record['cleanup_errors'].append('reap: ' + repr(error))
        try:
            record['linux_child_birth_absent'] = birth(child.pid) != record.get('linux_child_birth')
        except BaseException as error:
            record['cleanup_errors'].append('birth: ' + repr(error))
try:
    assert record.get('returncode') == 1 and not record.get('error') and not record['cleanup_errors'] and not signals
    assert record['linux_child_birth_absent']
    owner = json.loads((output / 'owner.json').read_text(encoding='utf-8-sig'))
    assert owner['status'] == 'completed' and owner['returncode'] == 1 and not owner['outer_stop'] and not owner['cleanup_errors']
    assert owner['owned_process_absent'] and owner['diagnostic_manifest_unchanged'] and owner['test_executable_unchanged']
    assert owner['source_revision'] == SOURCE and owner['live_source_revision'] == HEAD and owner['source_diagnostic_manifest_sha256'] == binding['manifest_identity']['sha256']
    record['diagnostic_build_manifest_authority_used'] = True
    record['native_build_gate_authority_used'] = False
    record['ci_native_gate_failed'] = True
    record['ci_native_gate_accepted'] = False
    record['mandatory_or_full_runtime_admitted'] = False
    assert owner['diagnostic_build_manifest_authority_used'] is True and owner['native_build_gate_authority_used'] is False and owner['ci_native_gate_failed'] is True and owner['ci_native_gate_accepted'] is False
    assert owner['numerical_control_passed'] is False and owner['system_inventory_passed'] is False and owner['loaded_driver_modules_recorded'] is True
    record['failure_context'] = read_failure_context(output/'gtest.xml')
    record['numerical_test_passed'] = False
    record['expected_completed_failure_captured'] = True
    record['loaded_provider_evidence_complete'] = owner['loaded_driver_modules_recorded'] is True
    if owner['loaded_driver_modules_recorded']:
        modules=json.loads((output/'driver-modules.json').read_text(encoding='utf-8-sig'))
        assert isinstance(modules,list) and len(modules) == 2
        loader=[row for row in modules if row['name'].lower() == 'vulkan-1.dll']
        assert len(loader) == 1 and any('amdvlk' in row['name'].lower() for row in modules)
        assert [row for row in modules if 'amdvlk' in row['name'].lower()] == [binding['expected_amd_module']]
        expected=verify_manifest()['files']['bin/windows-msvc/'+str(exe.parent.relative_to(stage)).replace('\\','/')+'/vulkan-1.dll']
        assert loader[0]['path'].replace('/','\\').casefold() == win(exe.parent/'vulkan-1.dll').casefold()
        assert {key:loader[0][key] for key in ['bytes','sha256']} == expected
    verify()
    assert not signals
    after = input_records()
    write(WORK / (kind+'-inputs-after.json'), after)
    assert before == after
    record['whole_inputs_unchanged'] = True
    record['accepted_owner_execution'] = True
except BaseException as error:
    record['postcondition_error'] = repr(error)
try:
    after = input_records()
    write(WORK / (kind+'-inputs-after.json'), after)
    record['whole_inputs_unchanged'] = before == after
    if not record['whole_inputs_unchanged']:
        record['accepted_owner_execution'] = False
except BaseException as error:
    record['accepted_owner_execution'] = False
    record['input_postcondition_error'] = repr(error)
record.update(completed=True, elapsed_seconds=time.monotonic() - started, signals=list(signals))
if signals:
    record['accepted_owner_execution'] = False
write(record_path, record)
print(json.dumps({k: record.get(k) for k in ('kind', 'returncode', 'elapsed_seconds', 'accepted_owner_execution', 'error', 'postcondition_error')}), flush=True)
sys.exit(0 if record['accepted_owner_execution'] else 1)
