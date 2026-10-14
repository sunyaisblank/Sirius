"""Own one original current native invocation; fail closed on uncertainty."""
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
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
kind, = sys.argv[1:]
assert kind in ('inventory', 'science', 'render')
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['source_revision']
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
    return subprocess.check_output(['wslpath', '-w', str(p.resolve())], text=True).strip()

def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except FileNotFoundError:
        return None

binding_hash = digest(binding_path)
assert binding['source_revision'] == HEAD and binding['whole_array_joins'] == 120
owner_review = json.loads((WORK / 'independent-owner-review.json').read_text())
assert owner_review['pass'] and owner_review['source_revision'] == HEAD
for name in ('run_current.py', 'run_native_control.ps1', 'emergency_native_cleanup.ps1', 'bind_current.py',
             'run_terminal_current.py', 'terminal-audit.ps1', 'plan.json'):
    p = WORK / name
    assert owner_review['reviewed_files'][name] == {'bytes': p.stat().st_size, 'sha256': digest(p)}
assert json.loads((WORK / 'native-preregister-readback.json').read_text())['body'] == (WORK / 'native-preregister-comment.md').read_text()
verify()
stage = ROOT / binding['stage']
output = WORK / 'native' / kind
temporary = WORK / 'temp' / kind
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
if kind != 'inventory':
    inventory = json.loads((WORK / 'native/inventory/stdout.log').read_text(encoding='utf-8-sig'))
    inventory_bridge = json.loads((WORK / 'inventory-bridge.json').read_text())
    assert inventory_bridge['accepted_owner_execution'] and inventory_bridge['completed'] and inventory_bridge['returncode'] == 0
    assert inventory_bridge['source_revision'] == HEAD and inventory_bridge['kind'] == 'inventory' and inventory_bridge['binding_sha256'] == binding_hash
    inventory_owner = json.loads((WORK / 'native/inventory/owner.json').read_text(encoding='utf-8-sig'))
    assert inventory_owner['source_revision'] == inventory_owner['live_source_revision'] == HEAD
    assert inventory_owner['source_build_gate_sha256'] == binding['gate_sha256'] and inventory_owner['system_inventory_passed']
    device = inventory['backends']['vulkan']['devices'][0]
    assert device['name'] == 'AMD Radeon 780M Graphics' and device['supports_fp64'] and device['rounds_fp64_to_nearest']
    assert device['preserves_fp32_denormals'] and device['rounds_fp32_to_nearest']
if kind == 'render':
    science_bridge = json.loads((WORK / 'science-bridge.json').read_text())
    assert science_bridge['accepted_owner_execution'] and science_bridge['completed'] and science_bridge['returncode'] == 0
    assert science_bridge['source_revision'] == HEAD and science_bridge['kind'] == 'science' and science_bridge['binding_sha256'] == binding_hash
    terminal = json.loads((WORK / 'terminal-science-accepted.json').read_text())
    assert terminal['pass'] and terminal['record']['bridge'] == science_bridge
    assert terminal['record']['owner']['source_revision'] == HEAD
    assert terminal['record']['gtest_sha256'] == digest(WORK / 'native/science/gtest.xml')
    review = json.loads((WORK / 'independent-science-review.json').read_text())
    assert review['pass'] and review['source_revision'] == HEAD and review['binding_sha256'] == binding_hash
    for rel in ('science-bridge.json', 'native/science/owner.json', 'native/science/gtest.xml',
                'native/science/driver-modules.json', 'terminal-science-accepted.json'):
        p = WORK / rel
        assert review['reviewed_files'][rel] == {'bytes': p.stat().st_size, 'sha256': digest(p)}
if kind != 'inventory':
    paths = [WORK / 'inventory-bridge.json', WORK / 'inventory-bridge.stdout', WORK / 'inventory-bridge.stderr',
             *[p for p in (WORK / 'native/inventory').rglob('*') if p.is_file()]]
    if kind == 'render':
        paths += [WORK / 'science-bridge.json', WORK / 'science-bridge.stdout', WORK / 'science-bridge.stderr',
                  WORK / 'independent-science-review.json', WORK / 'terminal-science-accepted.json',
                  WORK / 'terminal-science.json', WORK / 'terminal-science-bridge.json',
                  *[p for p in (WORK / 'native/science').rglob('*') if p.is_file()]]
    for p in paths:
        predecessors[str(p.relative_to(ROOT))] = {'bytes': p.stat().st_size, 'sha256': digest(p)}
    verify()
executable = 'sirius' if kind == 'inventory' else ('sirius_backend_tests' if kind == 'science' else 'sirius_render_tests')
exe = stage / json.loads((stage / 'generated/sirius/native_build_gate.json').read_text())['tested_artifacts'][executable]['path']
case = 'SystemInventory' if kind == 'inventory' else plan[kind + '_case']
seconds = plan['inventory_owner_seconds'] if kind == 'inventory' else plan[kind + '_owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'),
        '-Stage', win(exe.parent), '-Output', win(output), '-Filter', case, '-Seconds', str(seconds),
        '-RssBytes', str(plan['owner_rss_bytes']), '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage),
        '-ExpectedGateSha256', binding['gate_sha256'], '-LiveRevision', HEAD, '-ExecutableName', exe.name,
        '-InvocationKind', 'SystemInventory' if kind == 'inventory' else 'Control']
if kind == 'render':
    renders = ROOT / plan['render_output_root']
    assert not renders.exists()
    argv += ['-RenderDirectory', win(renders)]
cleanup_argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'),
                '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(WORK / (kind + '-emergency-cleanup.json'))]
record = {'source_revision': HEAD, 'kind': kind, 'argv': argv, 'binding_sha256': binding_hash,
          'linux_controller_pid': os.getpid(), 'linux_controller_birth': birth(os.getpid()),
          'predecessor_inputs': predecessors,
          'bridge_outer_seconds': plan['linux_bridge_outer_seconds'], 'windows_owner_seconds': seconds,
          'accepted_owner_execution': False, 'completed': False, 'cleanup_errors': [], 'full_qualification_claimed': False}
write(record_path, record)
started = time.monotonic()
child = None
try:
    assert not signals
    with (WORK / (kind + '-bridge.stdout')).open('wb') as out, (WORK / (kind + '-bridge.stderr')).open('wb') as err:
        child = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
        record['linux_child_pid'] = child.pid
        record['linux_child_birth'] = birth(child.pid)
        assert record['linux_child_birth'] is not None
        write(record_path, record)
        while child.poll() is None:
            assert not signals
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
    assert owner['source_revision'] == owner['live_source_revision'] == HEAD and owner['source_build_gate_sha256'] == binding['gate_sha256']
    if kind == 'inventory':
        assert owner['system_inventory_passed']
    else:
        assert owner['numerical_control_passed'] and owner['loaded_driver_modules_recorded']
        modules = json.loads((output / 'driver-modules.json').read_text(encoding='utf-8-sig'))
        assert {m['name'].lower() for m in modules} == {'amdvlk64.dll', 'vulkan-1.dll'}
        for m in modules:
            p = '/mnt/' + m['path'][0].lower() + m['path'][2:].replace('\\', '/')
            expected = next(value for path, value in binding['providers'].items() if path.lower() == p.lower())
            assert {'bytes': m['bytes'], 'sha256': m['sha256']} == expected
        xml = ET.parse(output / 'gtest.xml')
        assert len(xml.findall('.//testcase')) == 1
        for node in [xml.getroot(), *xml.findall('.//testsuite')]:
            for key in ('failures', 'errors', 'skipped', 'disabled'):
                assert int(node.get(key, '0')) == 0
        node, = xml.findall('.//testcase')
        assert node.get('classname') + '.' + node.get('name') == case and node.get('status') == 'run' and node.get('result') == 'completed'
        assert not any(xml.findall('.//' + tag) for tag in ('failure', 'error', 'skipped'))
    verify()
    assert not signals
    record['accepted_owner_execution'] = True
except BaseException as error:
    record['postcondition_error'] = repr(error)
record.update(completed=True, elapsed_seconds=time.monotonic() - started, signals=list(signals))
if signals:
    record['accepted_owner_execution'] = False
write(record_path, record)
print(json.dumps({k: record.get(k) for k in ('kind', 'returncode', 'elapsed_seconds', 'accepted_owner_execution', 'error', 'postcondition_error')}), flush=True)
sys.exit(0 if record['accepted_owner_execution'] else 1)
