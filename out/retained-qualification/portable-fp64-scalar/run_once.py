#!/usr/bin/env python3
"""One authorized bounded default-Lavapipe dispatch; terminal evidence required."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
import struct
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REPORT = HERE/'gpu-report.json'
TIMEOUT_SECONDS = 180
RSS_GUARD_KIB = 4 * 1024 * 1024

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save(report):
    temporary = REPORT.with_suffix('.tmp')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(REPORT)

def process_tree_rss(pid):
    pending, seen, total = [pid], set(), 0
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            text = Path(f'/proc/{current}/status').read_text()
            for line in text.splitlines():
                if line.startswith('VmRSS:'):
                    total += int(line.split()[1])
            children = Path(f'/proc/{current}/task/{current}/children').read_text().split()
            pending.extend(map(int, children))
        except (FileNotFoundError, ProcessLookupError):
            pass
    return total

assert not REPORT.exists(), 'single-run record already exists; do not dispatch again'
assert not (HERE/'gpu-actual.bin').exists(), 'existing output would violate fresh-result binding'
preparation = json.loads((HERE/'preparation-report.json').read_text())
runtime_bound_sources = {name: value for name, value in preparation['bound_sources'].items()
                         if not name.endswith('.a')}
for relative, expected in runtime_bound_sources.items():
    assert sha(ROOT/relative) == expected, f'preparation input changed: {relative}'
# Driver compiler/shading/cache options must genuinely remain default. Device
# selection is the sole environment change, binding this execution to Lavapipe.
driver_options = {key: value for key, value in os.environ.items()
                  if key.startswith(('GALLIVM_', 'LP_', 'MESA_'))}
assert not driver_options, f'nondefault driver options: {driver_options}'
icd = Path('/usr/share/vulkan/icd.d/lvp_icd.json')
icd_data = json.loads(icd.read_text())
library = Path(icd_data['ICD']['library_path'])
if not library.is_absolute():
    # Bare ICD library names resolve through the loader's library search. Bind
    # the exact cache entry now, then require dladdr's actual resident identity.
    soname = str(library)
    entries = subprocess.run(['/sbin/ldconfig', '-p'], capture_output=True, text=True, check=True).stdout.splitlines()
    candidates = [line.split('=>', 1)[1].strip() for line in entries
                  if line.strip().startswith(soname + ' ') and 'x86-64' in line and '=>' in line]
    assert len(candidates) == 1
    library = Path(candidates[0]).resolve()
assert library.is_absolute() and library.is_file()
environment = os.environ.copy()
environment['VK_ICD_FILENAMES'] = str(icd)
environment['SIRIUS_VULKAN_DEVICE'] = '0'
assert not environment.get('VK_DRIVER_FILES') and not environment.get('VK_ADD_DRIVER_FILES')
command = [str(HERE/'gpu_probe'), str(HERE/'probe.spv'), str(HERE/'shader-input.bin'),
           str(HERE/'shader-expected.bin'), str(HERE/'gpu-actual.bin'), '--dispatch']
report = {
    'scope': 'one finite isolated raw-word prototype; not retained-stage or universal binary64 qualification',
    'status': 'prepared', 'started_utc': timestamp(), 'command': command,
    'preparation_sha256': sha(HERE/'preparation-report.json'),
    'runner_source_sha256': sha(Path(__file__)),
    'source_head': preparation['source_head'], 'bounds': {'timeout_seconds': TIMEOUT_SECONDS,
                 'process_tree_rss_guard_kib': RSS_GUARD_KIB, 'explicit_device_bytes_limit': 8*1024*1024},
    'driver_options': {'GALLIVM_': 'unset', 'LP_': 'unset', 'MESA_': 'unset'},
    'selection_environment': {'VK_ICD_FILENAMES': str(icd), 'SIRIUS_VULKAN_DEVICE': '0'},
    'icd': {'path': str(icd), 'sha256': sha(icd), 'contents': icd_data},
    'driver_library': {'path': str(library), 'sha256': sha(library)},
    'executable_sha256': sha(HERE/'gpu_probe'), 'module_sha256': sha(HERE/'probe.spv'),
    'raw_input_sha256': sha(HERE/'shader-input.bin'), 'expected_sha256': sha(HERE/'shader-expected.bin'),
    'process_peak_sampled_rss_kib': 0, 'terminated_by_supervisor': False,
    'timeout': False, 'rss_guard_triggered': False, 'terminal_exit_code': None,
}
save(report)
started = time.monotonic()
with (HERE/'gpu-stdout.jsonl').open('w') as stdout, (HERE/'gpu-stderr.log').open('w') as stderr:
    process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr,
                               start_new_session=True)
    report['pid'] = process.pid
    try:
        report['process_start_ticks'] = int(Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()[19])
        report['mapped_executable_sha256'] = sha(Path(f'/proc/{process.pid}/exe'))
        assert report['mapped_executable_sha256'] == report['executable_sha256']
        report['status'] = 'running'
        save(report)
        while process.poll() is None:
            elapsed = time.monotonic() - started
            rss = process_tree_rss(process.pid)
            report['process_peak_sampled_rss_kib'] = max(report['process_peak_sampled_rss_kib'], rss)
            if elapsed > TIMEOUT_SECONDS or rss > RSS_GUARD_KIB:
                report['timeout'] = elapsed > TIMEOUT_SECONDS
                report['rss_guard_triggered'] = rss > RSS_GUARD_KIB
                report['terminated_by_supervisor'] = True
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait(timeout=5)
                break
            time.sleep(0.1)
        report['terminal_exit_code'] = process.wait(timeout=5)
    except BaseException as error:
        report['runner_error'] = repr(error)
        if process.poll() is None:
            report['terminated_by_supervisor'] = True
            os.killpg(process.pid, signal.SIGKILL)
        report['terminal_exit_code'] = process.wait(timeout=5)
    finally:
        report['whole_child_wall_seconds'] = time.monotonic() - started
        report['finished_utc'] = timestamp()
        report['status'] = 'terminal'
        save(report)

lines = (HERE/'gpu-stdout.jsonl').read_text().splitlines()
observations = [json.loads(line) for line in lines if line.startswith('{')]
report['raw_observations'] = observations
report['stdout_sha256'] = sha(HERE/'gpu-stdout.jsonl')
report['stderr_sha256'] = sha(HERE/'gpu-stderr.log')
report['driver_library_unchanged'] = sha(library) == report['driver_library']['sha256']
report['bound_sources_unchanged'] = all(sha(ROOT/name) == expected for name, expected in runtime_bound_sources.items())
report['library_binding'] = 'static libraries bound at link; runtime executable hash is authoritative; later unrelated preset rebuilds excluded'
report['python_comparison'] = {'completed': False}
if (HERE/'gpu-actual.bin').exists():
    actual = (HERE/'gpu-actual.bin').read_bytes()
    expected = (HERE/'shader-expected.bin').read_bytes()
    report['actual_sha256'] = hashlib.sha256(actual).hexdigest()
    if len(actual) == len(expected) == 522804 * 4:
        aa = struct.unpack('<522804I', actual)
        ee = struct.unpack('<522804I', expected)
        mismatches = [i for i, (a, b) in enumerate(zip(aa, ee)) if a != b]
        untouched = sum(a == (b ^ 0xffffffff) for a, b in zip(aa, ee))
        report['python_comparison'] = {'completed': True, 'observed_words': len(aa),
                      'cases': 174268, 'mismatches': len(mismatches), 'first_mismatch_word_indices': mismatches[:8],
                      'untouched_complement_words': untouched,
                      'operation_counts': preparation['corpus']['operation_counts'],
                      'assisted_Add_Sub_Mul_cases': sum(preparation['corpus']['operation_counts'][str(x)] for x in [0,1,2])}
    else:
        report['python_comparison']['output_bytes'] = len(actual)
identities = [v for v in observations if 'device' in v]
report['actual_resident_driver_matches'] = (len(identities) == 1 and
      Path(identities[0]['resident_driver_path']).resolve() == library.resolve() and
      sha(Path(identities[0]['resident_driver_path'])) == report['driver_library']['sha256'])
completions = [v for v in observations if 'completed_dispatches' in v]
starts = [v for v in observations if 'dispatch_started' in v]
report['pass'] = (
    report['terminal_exit_code'] == 0 and not report['terminated_by_supervisor'] and
    report['driver_library_unchanged'] and report['actual_resident_driver_matches'] and report['bound_sources_unchanged'] and
    len(identities) == len(completions) == len(starts) == 1 and
    identities[0]['kind'] == 'software' and 'llvmpipe' in identities[0]['device'] and
    identities[0]['fp64'] == identities[0]['RTE64'] == 1 and
    starts[0]['dispatch_started'] == 1 and starts[0]['complement_initialized_words'] == 522804 and
    completions[0]['completed_dispatches'] == 1 and completions[0]['observed_words'] == 522804 and
    completions[0]['mismatches'] == completions[0]['untouched_complement_words'] == 0 and
    report['python_comparison']['completed'] and report['python_comparison']['mismatches'] == 0 and
    report['python_comparison']['untouched_complement_words'] == 0 and
    report.get('actual_sha256') == report['expected_sha256'])
save(report)
print(json.dumps({'pass': report['pass'], 'terminal_exit_code': report['terminal_exit_code'],
                  'wall_seconds': report['whole_child_wall_seconds'], 'comparison': report['python_comparison']}))
raise SystemExit(0 if report['pass'] else 1)
