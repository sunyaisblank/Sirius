"""Run the reviewed isolated stateful numerical gate once."""
from pathlib import Path
import hashlib
import json
import os
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
QUERY = TASK / 'numerical'


def binding(path):
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(items):
    for item in items:
        assert binding(Path(item['path'])) == item, item['path']


def win(path):
    return r'\\wsl.localhost\Ubuntu' + str(path.resolve()).replace('/', '\\')


preparation = json.loads((QUERY / 'preparation.json').read_text())
assert preparation['device_dispatches_authorized'] == 54
assert preparation['query_pipelines'] == 2
assert preparation['seconds_guard'] == 90
assert preparation['rss_bytes_guard'] == 4294967296
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == preparation['revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
assert not (QUERY / 'bridge-query.json').exists()
assert not (QUERY / 'run-query').exists()
reviews = [QUERY / name for name in ('root-prelaunch-review.json', 'independent-prelaunch-review.json')]
for path in reviews:
    review = json.loads(path.read_text())
    assert review['clear_for_single_numerical_gate'] is True
    assert review['module_sha256'] == '721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27'
    verify(review['apparatus_bindings'])
inputs = preparation['input_bindings'] + preparation['provider_bindings']
inputs += [binding(QUERY / 'preparation.json'), binding(Path(__file__).resolve())]
inputs += [binding(path) for path in reviews]
verify(inputs)
args = [str(Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe')),
        '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File',
        win(QUERY / 'query-owner.ps1'), '-Stage', win(QUERY), '-Output', win(QUERY / 'run-query'),
        '-Filter', 'isolated-native-paired-transport', '-Expected', 'isolated-native-paired-transport',
        '-Seconds', '90', '-RssBytes', '4294967296',
        '-TemporaryDirectory', win(QUERY / 'query-owner-temp')]
record = {'revision': preparation['revision'], 'arguments': args, 'input_bindings': inputs,
          'linux_owner_pid': os.getpid(),
          'linux_owner_start_ticks': Path('/proc/self/stat').read_text().split(') ', 1)[1].split()[19],
          'expected_device_dispatches': 54, 'completed': False,
          'bridge_outer_seconds': 180,
          'scope': preparation['scope']}
(QUERY / 'bridge-query.json').write_text(json.dumps(record, indent=2) + '\n')
clock = time.monotonic()
cleanup_args = [args[0], '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass',
                '-File', win(QUERY / 'emergency-cleanup.ps1'),
                '-TemporaryDirectory', win(QUERY / 'query-owner-temp'),
                '-Output', win(QUERY / 'run-query'), '-TaskWorkspace', win(QUERY),
                '-Receipt', win(QUERY / 'emergency-cleanup.json')]


def interrupted(number, _frame):
    raise InterruptedError(f'bridge_signal_{number}')


for number in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
    signal.signal(number, interrupted)
bridge_error = cleanup_error = None
outer_stop = False
code = None
child = None
with (QUERY / 'bridge-query.stdout').open('wb') as stdout, (QUERY / 'bridge-query.stderr').open('wb') as stderr:
    try:
        record['spawn_attempted'] = True
        (QUERY / 'bridge-query.json').write_text(json.dumps(record, indent=2) + '\n')
        child = subprocess.Popen(args, cwd=ROOT, stdout=stdout, stderr=stderr, start_new_session=True)
        record['linux_bridge_child_pid'] = child.pid
        record['linux_bridge_child_start_ticks'] = Path('/proc', str(child.pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
        (QUERY / 'bridge-query.json').write_text(json.dumps(record, indent=2) + '\n')
        code = child.wait(timeout=180)
    except BaseException as error:
        bridge_error = repr(error)
        outer_stop = True
    finally:
        for number in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(number, signal.SIG_IGN)
        if outer_stop or (child is not None and child.poll() is None):
            try:
                with (QUERY / 'emergency-cleanup.stdout').open('wb') as out, (QUERY / 'emergency-cleanup.stderr').open('wb') as err:
                    result = subprocess.run(cleanup_args, stdout=out, stderr=err, timeout=60)
                record['emergency_cleanup_returncode'] = result.returncode
                assert result.returncode == 0
                receipt = json.loads((QUERY / 'emergency-cleanup.json').read_text(encoding='utf-8-sig'))
                assert receipt['recorded_births_absent'] and not receipt['errors']
            except BaseException as error:
                cleanup_error = repr(error)
        if child is not None:
            try:
                code = child.wait(timeout=15)
            except BaseException as error:
                cleanup_error = cleanup_error or repr(error)
                try:
                    # Signal only the exact still-owned Linux child. Native cleanup
                    # above uses recorded Windows births independently of this PID.
                    current = Path('/proc', str(child.pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
                    assert current == record['linux_bridge_child_start_ticks']
                    os.killpg(child.pid, signal.SIGKILL)
                except (FileNotFoundError, ProcessLookupError):
                    pass
                except BaseException as error:
                    cleanup_error = repr(error)
                try:
                    code = child.wait(timeout=10)
                except BaseException as error:
                    cleanup_error = repr(error)
                    code = child.returncode
            try:
                current = Path('/proc', str(child.pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
            except (FileNotFoundError, ProcessLookupError):
                current = None
            saved_birth = record.get('linux_bridge_child_start_ticks')
            record['linux_bridge_child_birth_absent'] = saved_birth is not None and current != saved_birth
        else:
            record['linux_bridge_child_birth_absent'] = False
record.update(returncode=code, elapsed_seconds=time.monotonic() - clock, completed=True,
              bridge_error=bridge_error, cleanup_error=cleanup_error, outer_stop=outer_stop)
try:
    verify(inputs)
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == preparation['revision']
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
    record['all_selected_inputs_unchanged'] = True
except Exception as error:
    record['all_selected_inputs_unchanged'] = False
    record['post_seal_error'] = str(error)
(QUERY / 'bridge-query.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({key: record[key] for key in
                 ('returncode', 'elapsed_seconds', 'all_selected_inputs_unchanged')}), flush=True)
accepted = (code == 0 and not outer_stop and bridge_error is None and cleanup_error is None and
            record['linux_bridge_child_birth_absent'] and record['all_selected_inputs_unchanged'])
raise SystemExit(0 if accepted else 1)
