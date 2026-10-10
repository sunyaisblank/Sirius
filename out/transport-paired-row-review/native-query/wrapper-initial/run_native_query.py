"""Run the reviewed isolated compiler query once; never dispatch a shader."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
QUERY = TASK / 'native-query'


def binding(path):
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(items):
    for item in items:
        assert binding(Path(item['path'])) == item, item['path']


def win(path):
    return r'\\wsl.localhost\Ubuntu' + str(path.resolve()).replace('/', '\\')


preparation = json.loads((QUERY / 'preparation.json').read_text())
assert preparation['device_dispatches_authorized'] == 0
assert preparation['query_pipelines'] == 1
assert preparation['seconds_guard'] == 90
assert preparation['rss_bytes_guard'] == 4294967296
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == preparation['revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
assert not (QUERY / 'bridge-query.json').exists()
assert not (QUERY / 'run-query').exists()
reviews = [QUERY / name for name in ('root-prelaunch-review.json', 'independent-prelaunch-review.json')]
for path in reviews:
    review = json.loads(path.read_text())
    assert review['clear_for_single_zero_dispatch_query'] is True
    assert review['module_sha256'] == '721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27'
    verify(review['apparatus_bindings'])
inputs = preparation['input_bindings'] + preparation['provider_bindings']
inputs += [binding(QUERY / 'preparation.json'), binding(Path(__file__).resolve())]
inputs += [binding(path) for path in reviews]
verify(inputs)
args = [str(Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe')),
        '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File',
        win(QUERY / 'query-owner.ps1'), '-Stage', win(QUERY), '-Output', win(QUERY / 'run-query'),
        '-Filter', 'ordinary-native-shader-info', '-Expected', 'ordinary-native-shader-info',
        '-Seconds', '90', '-RssBytes', '4294967296',
        '-TemporaryDirectory', win(QUERY / 'query-owner-temp')]
record = {'revision': preparation['revision'], 'arguments': args, 'input_bindings': inputs,
          'linux_owner_pid': os.getpid(),
          'linux_owner_start_ticks': Path('/proc/self/stat').read_text().split(') ', 1)[1].split()[19],
          'device_dispatches': 0, 'completed': False,
          'scope': preparation['scope']}
(QUERY / 'bridge-query.json').write_text(json.dumps(record, indent=2) + '\n')
clock = time.monotonic()
with (QUERY / 'bridge-query.stdout').open('wb') as stdout, (QUERY / 'bridge-query.stderr').open('wb') as stderr:
    child = subprocess.run(args, cwd=ROOT, stdout=stdout, stderr=stderr)
record.update(returncode=child.returncode, elapsed_seconds=time.monotonic() - clock, completed=True)
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
raise SystemExit(child.returncode if record['all_selected_inputs_unchanged'] else 1)
