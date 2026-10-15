"""Own exact diagnostic-manifest Windows inventory or original smoke; no qualification."""
from pathlib import Path, PurePosixPath
import hashlib
import json
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
kind, = sys.argv[1:]
assert kind in ('inventory', 'smoke')
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['runner_source_revision']
SOURCE = plan['source_revision']
assert plan['admitted'] is True and binding['admitted'] is True
assert SOURCE == HEAD == '23d406b8ca99695f362b88e0c5f9a706f08a7a39'
assert plan['stage'] == binding['stage'] == 'bin/windows-msvc/diagnostic-23d406b'
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

def verify_manifest():
    stage = ROOT / plan['stage']
    manifest_path = stage / 'generated/sirius/windows-diagnostic-build-products.json'
    assert {'bytes': manifest_path.stat().st_size, 'sha256': digest(manifest_path)} == binding['manifest_identity']
    assert binding['files'][str(manifest_path.relative_to(ROOT))] == binding['manifest_identity']
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    for key, value in {'schema':'sirius-windows-diagnostic-build-products-v1', 'source_revision':SOURCE,
                       'repository':'sunyaisblank/Sirius', 'run_id':'38701994420', 'run_attempt':'1'}.items():
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
assert plan['inventory_owner_seconds'] == 60 and plan['smoke_owner_seconds'] == 320 and plan['owner_rss_bytes'] == 4294967296
assert plan['linux_bridge_outer_seconds'] == 480 and plan['terminal_observer_seconds'] == 60
assert plan['minimum_wsl_available_kib'] == 2097152

verify()
stage = ROOT / plan['stage']
output = WORK / 'native' / kind
temporary = WORK / 'temp' / kind
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
assert not (WORK / (kind+'-inputs-before.json')).exists() and not (WORK / (kind+'-inputs-after.json')).exists()
exe = stage / ('src/sirius/app/Release/sirius.exe' if kind == 'inventory' else 'tests/backend/Release/sirius_backend_tests.exe')
case = 'SystemInventory' if kind == 'inventory' else 'VulkanBackend.SlangKernelMatchesCpuReference'
seconds = plan[kind+'_owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'),
        '-Stage', win(exe.parent), '-Output', win(output), '-Filter', case, '-Seconds', str(seconds),
        '-RssBytes', str(plan['owner_rss_bytes']), '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage),
        '-ExpectedManifestSha256', binding['manifest_identity']['sha256'], '-LiveRevision', HEAD, '-ExecutableName', exe.name,
        '-InvocationKind', 'SystemInventory' if kind == 'inventory' else 'Control']
cleanup_argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'),
                '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(WORK / (kind + '-emergency-cleanup.json'))]
assert argv == binding['actions'][kind]['argv'] and cleanup_argv == binding['actions'][kind]['cleanup_argv']
if kind == 'smoke':
    prior = json.loads((WORK/'terminal-inventory-accepted.json').read_text())
    assert prior['pass'] is True and prior['case'] == 'inventory'
    assert prior['record']['bridge']['binding_sha256'] == binding_hash
    assert prior['terminal_receipt_sha256'] == digest(WORK/'terminal-inventory.json')
    assert prior['terminal_bridge_sha256'] == digest(WORK/'terminal-inventory-bridge.json')
    prior_owner = prior['record']['owner']
    assert prior_owner == json.loads((WORK/'native/inventory/owner.json').read_text(encoding='utf-8-sig'))
    assert prior_owner['source_revision'] == SOURCE and prior_owner['system_inventory_passed'] is True
    for relative in ['inventory-bridge.json', 'inventory-inputs-before.json', 'inventory-inputs-after.json',
                     'native/inventory/owner.json', 'native/inventory/stdout.log', 'native/inventory/stderr.log', 'native/inventory/memory.json',
                     'terminal-inventory.json', 'terminal-inventory-bridge.json', 'terminal-inventory-accepted.json']:
        path = WORK/relative
        predecessors[str(path.relative_to(ROOT))] = {'bytes':path.stat().st_size, 'sha256':digest(path)}
    if (WORK/'native/inventory/driver-modules.json').exists():
        path=WORK/'native/inventory/driver-modules.json';predecessors[str(path.relative_to(ROOT))]={'bytes':path.stat().st_size,'sha256':digest(path)}
    verify()
before = input_records()
write(WORK / (kind+'-inputs-before.json'), before)
record = {'source_revision': SOURCE, 'runner_source_revision': HEAD, 'kind': kind, 'argv': argv, 'binding_sha256': binding_hash,
          'linux_controller_pid': os.getpid(), 'linux_controller_birth': birth(os.getpid()),
          'predecessor_inputs': predecessors,
          'bridge_outer_seconds': plan['linux_bridge_outer_seconds'], 'windows_owner_seconds': seconds, 'minimum_wsl_available_kib': plan['minimum_wsl_available_kib'],
          'resource_scope': plan['resource_semantics'],
          'accepted_owner_execution': False, 'completed': False, 'cleanup_errors': [], 'full_qualification_claimed': False}
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
            assert time.monotonic() - started < plan['linux_bridge_outer_seconds'], 'bridge outer diagnostic stop'
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
    assert record.get('returncode') == 0 and not record.get('error') and not record['cleanup_errors'] and not signals
    assert record['linux_child_birth_absent']
    owner = json.loads((output / 'owner.json').read_text(encoding='utf-8-sig'))
    assert owner['status'] == 'completed' and owner['returncode'] == 0 and not owner['outer_stop'] and not owner['cleanup_errors']
    assert owner['owned_process_absent'] and owner['source_product_manifest_unchanged'] and owner['test_executable_unchanged']
    assert owner['source_revision'] == SOURCE and owner['live_source_revision'] == HEAD and owner['source_product_manifest_sha256'] == binding['manifest_identity']['sha256']
    assert owner['native_build_or_mandatory_gate_admitted'] is False
    if kind == 'inventory':
        assert owner['system_inventory_passed'] is True and owner['numerical_control_passed'] is False
        inventory = json.loads((output / 'stdout.log').read_text(encoding='utf-8-sig'))
        assert inventory['platform'] == 'Windows' and inventory['wsl2'] is False
        assert inventory['backends']['vulkan']['device_count'] == 1
        device = inventory['backends']['vulkan']['devices'][0]
        assert device['name'] == 'AMD Radeon 780M Graphics' and device['kind'] in ('integrated', 'discrete')
        assert inventory['backends']['vulkan']['selected_device_index'] == device['index'] == 0
        assert not any(word in (device['driver_name'] + ' ' + device['driver_info']).lower() for word in ('dozen', 'dzn', 'llvmpipe', 'lavapipe'))
    else:
        assert owner['numerical_control_passed'] is True and owner['system_inventory_passed'] is False and owner['loaded_driver_modules_recorded'] is True
        xml=ET.parse(output/'gtest.xml').getroot();cases=xml.findall('.//testcase')
        assert int(xml.get('tests','-1')) == 1 and len(cases) == 1
        assert all(int(node.get(key,'0')) == 0 for node in [xml,*xml.findall('.//testsuite')] for key in ['failures','errors','skipped','disabled'])
        assert not any(xml.findall('.//'+key) for key in ['failure','error','skipped'])
        assert cases[0].get('classname')+'.'+cases[0].get('name') == case and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    record['loaded_provider_evidence_complete'] = owner['loaded_driver_modules_recorded'] is True
    if owner['loaded_driver_modules_recorded']:
        modules=json.loads((output/'driver-modules.json').read_text(encoding='utf-8-sig'))
        assert isinstance(modules,list) and len(modules) == 2
        loader=[row for row in modules if row['name'].lower() == 'vulkan-1.dll']
        assert len(loader) == 1 and any('amdvlk' in row['name'].lower() for row in modules)
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
