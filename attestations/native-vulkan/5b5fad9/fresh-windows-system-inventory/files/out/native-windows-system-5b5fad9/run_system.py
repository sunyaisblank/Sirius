"""Own one archived-export native Windows inventory; no test/render/qualification."""
from pathlib import Path
import hashlib
import json
import os
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
assert kind == 'inventory'
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['runner_source_revision']
SOURCE = plan['source_revision']
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

def win(p):
    p = p.resolve()
    assert p.is_relative_to(ROOT)
    return plan['windows_checkout_unc'] + '\\' + str(p.relative_to(ROOT)).replace('/', '\\')

def available_kib():
    rows = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(rows['MemAvailable'].split()[0])

def input_records():
    paths = [ROOT / rel for rel in binding['files']] + [binding_path, WORK / 'preparation.json']
    return {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in paths}

def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except FileNotFoundError:
        return None

binding_hash = digest(binding_path)
assert binding['source_revision'] == SOURCE and binding['runner_source_revision'] == HEAD
assert plan['inventory_owner_seconds'] == 60 and plan['owner_rss_bytes'] == 4294967296
assert plan['minimum_wsl_available_kib'] == 2097152

verify()
stage = ROOT / binding['stage']
output = WORK / 'native' / kind
temporary = WORK / 'temp' / kind
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
assert not (WORK / 'inputs-before.json').exists() and not (WORK / 'inputs-after.json').exists()
executable = 'sirius'
exe = stage / json.loads((stage / 'generated/sirius/native_build_gate.json').read_text())['tested_artifacts'][executable]['path']
case = 'SystemInventory'
seconds = plan['inventory_owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'),
        '-Stage', win(exe.parent), '-Output', win(output), '-Filter', case, '-Seconds', str(seconds),
        '-RssBytes', str(plan['owner_rss_bytes']), '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage),
        '-ExpectedGateSha256', binding['gate_sha256'], '-LiveRevision', HEAD, '-ExecutableName', exe.name,
        '-InvocationKind', 'SystemInventory' if kind == 'inventory' else 'Control']
cleanup_argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'),
                '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(WORK / (kind + '-emergency-cleanup.json'))]
assert argv == binding['inventory_argv'] and cleanup_argv == binding['inventory_cleanup_argv']
before = input_records()
write(WORK / 'inputs-before.json', before)
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
    assert owner['owned_process_absent'] and owner['source_build_gate_unchanged'] and owner['test_executable_unchanged']
    assert owner['source_revision'] == SOURCE and owner['live_source_revision'] == HEAD and owner['source_build_gate_sha256'] == binding['gate_sha256']
    assert owner['system_inventory_passed'] and not owner['numerical_control_passed']
    inventory = json.loads((output / 'stdout.log').read_text(encoding='utf-8-sig'))
    assert inventory['platform'] == 'Windows' and inventory['wsl2'] is False
    assert inventory['backends']['vulkan']['device_count'] == 1
    device = inventory['backends']['vulkan']['devices'][0]
    assert device['name'] == 'AMD Radeon 780M Graphics' and device['kind'] in ('integrated', 'discrete')
    assert not any(word in (device['driver_name'] + ' ' + device['driver_info']).lower() for word in ('dozen', 'dzn', 'llvmpipe', 'lavapipe'))
    verify()
    assert not signals
    after = input_records()
    write(WORK / 'inputs-after.json', after)
    assert before == after
    record['whole_inputs_unchanged'] = True
    record['accepted_owner_execution'] = True
except BaseException as error:
    record['postcondition_error'] = repr(error)
try:
    after = input_records()
    write(WORK / 'inputs-after.json', after)
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
