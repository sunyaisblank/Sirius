from pathlib import Path
import ctypes, hashlib, importlib.util, json, os, re, signal, subprocess, sys, time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
plan = json.loads((work / 'plan.json').read_text())
spec = importlib.util.spec_from_file_location('preserved_owned_guard', work / 'owned_guard.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)
guard.ROOT = root
guard.BUILD = root / 'bin/linux-gcc'
guard.RSS_KIB = plan['outer_limits']['rss_kib']
guard.RESERVE_KIB = plan['outer_limits']['reserve_kib']
assert ctypes.CDLL(None, use_errno=True).prctl(36, 1, 0, 0, 0) == 0
for number in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
    signal.signal(number, guard.signal_stop)

source = guard.source_identity()
assert source == {'revision': plan['source_revision'], 'status': ''}
assert guard.query(['git', 'rev-parse', 'HEAD^{tree}']).decode().strip() == plan['source_tree']
assert not (work / 'science').exists()
folder = work / 'science'
folder.mkdir()
backend = root / 'bin/linux-gcc/tests/backend/sirius_backend_tests'
assert backend.stat().st_size == plan['backend_expected_bytes']
assert guard.digest(backend) == plan['backend_expected_sha256']
selected = json.loads((work / 'selected-registration.json').read_text())
assert len(selected) == 1 and selected[0]['name'] == plan['case']
assert selected[0]['command'] == [str(backend), '--gtest_filter=' + plan['case'],
                                  '--gtest_also_run_disabled_tests']
properties = {e['name']: e['value'] for e in selected[0]['properties']}
assert properties['LABELS'] == ['Correctness', 'Mandatory']
assert properties['RESOURCE_LOCK'] == ['sirius_vulkan_device']
assert 'TIMEOUT' not in properties
restoration = root / 'attestations/native-vulkan/0fc9eb6/paired-transport-restoration'
build_receipt = json.loads((restoration / 'diagnostic/build-and-arrays.json').read_text())
bindings = json.loads((restoration / 'diagnostic/linux-readonly-payloads.json').read_text())
assert build_receipt['source_revision'] == plan['source_revision'] and build_receipt['returncode'] == 0
assert bindings['source_revision'] == plan['source_revision'] and bindings['pass_']
assert bindings['actual_consumers'][str(backend.relative_to(root))]['sha256'] == guard.digest(backend)
assert guard.digest(root / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h') == build_receipt['whole_header_sha256']

env = os.environ.copy()
overrides = [key for key, value in env.items() if value and
             (key.startswith(('SIRIUS_', 'VK_', 'LP_', 'MESA_', 'GALLIVM_', 'GTEST_', 'CTEST_'))
              or key in {'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'})]
assert not overrides, ('unregistered inherited execution overrides', sorted(overrides))
env['GTEST_OUTPUT'] = 'xml:' + str(folder / 'gtest.xml')
tracked = [root / os.fsdecode(p) for p in guard.query(['git', 'ls-files', '-z']).split(b'\0') if p]
files = [*tracked, Path(__file__).resolve(), work / 'owned_guard.py', work / 'plan.json',
         work / 'selected-registration.json', backend, root / 'bin/linux-gcc/CMakeCache.txt',
         root / 'bin/linux-gcc/compile_commands.json', root / 'bin/linux-gcc/build.ninja',
         root / 'bin/linux-gcc/generated/sirius/alignment_receipt.json',
         root / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
         restoration / 'manifest.json', restoration / 'diagnostic/build-and-arrays.json',
         restoration / 'diagnostic/linux-readonly-payloads.json',
         root / plan['reused_endpoint_evidence']]
files.extend((root / 'bin/linux-gcc').rglob('CTestTestfile.cmake'))
files.extend((root / 'bin/linux-gcc').rglob('*_tests.cmake'))
files.extend((root / 'bin/linux-gcc').rglob('*_include.cmake'))
before = {str(p.relative_to(root)): guard.identity(p) for p in files}
guard.write(work / 'inputs-before.json', before)
assert guard.available_kib() >= guard.RESERVE_KIB and guard.STOP is None
argv = ['ctest', '--preset', 'linux-gcc', '--no-tests=error', '--output-on-failure',
        '-R', '^' + re.escape(plan['case']) + '$', '--output-junit', str(folder / 'ctest.xml')]
report = {'source': source, 'case': plan['case'], 'argv': argv,
          'scope': plan['scope'], 'limits': plan['outer_limits'],
          'declared_output_environment': {'GTEST_OUTPUT': env['GTEST_OUTPUT']},
          'controller': guard.proc(os.getpid()), 'boot_id': Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
          'passed': False, 'cleanup_errors': [], 'observed_births': [], 'stop_reason': None}
start = time.monotonic()
child = None
ticks = None
observed = {}
reaped = []
reason = None
code = None
peak = 0
minimum = guard.available_kib()
with (folder / 'stdout.log').open('wb') as out, (folder / 'stderr.log').open('wb') as err, (folder / 'samples.jsonl').open('w') as samples:
    try:
        child = subprocess.Popen(argv, cwd=root, env=env, stdout=out, stderr=err, start_new_session=True)
        ticks = guard.proc(child.pid)['start_ticks']
        observed[child.pid] = ticks
        report.update(child_pid=child.pid, child_start_ticks=ticks, status='running')
        guard.write(folder / 'owner.json', report)
        while child.poll() is None:
            members = guard.group_members(child.pid)
            for item in members:
                observed[item['pid']] = item['start_ticks']
            available = guard.available_kib()
            rss = sum(item['rss_kib'] for item in members) + guard.proc(os.getpid())['rss_kib']
            elapsed = time.monotonic() - start
            peak = max(peak, rss)
            minimum = min(minimum, available)
            samples.write(json.dumps({'elapsed_seconds': elapsed, 'rss_kib': rss,
                                      'available_kib': available, 'members': members}) + '\n')
            samples.flush()
            reason = guard.STOP or ('diagnostic_timeout' if elapsed >= plan['outer_limits']['seconds']
                else 'diagnostic_sampled_rss_guard' if rss > guard.RSS_KIB
                else 'diagnostic_host_reserve_guard' if available < guard.RESERVE_KIB else None)
            if reason:
                break
            time.sleep(.25)
    except BaseException as error:
        reason = reason or 'observer_exception'
        report['observer_error'] = repr(error)
    finally:
        if child is not None:
            try:
                remaining = guard.group_members(child.pid)
                if remaining:
                    reason = reason or 'owned_group_remaining'
                    guard.terminate_group(child, ticks, reaped)
            except BaseException as error:
                reason = reason or 'cleanup_observer_exception'
                report['cleanup_errors'].append(repr(error))
                # An unreaped Popen child cannot have its PID reused. Check its
                # actual kernel group/session before a bounded emergency stop.
                # This fallback retains a failed verdict if procfs observation
                # or the inherited birth-based cleanup failed.
                if child.returncode is None:
                    try:
                        assert os.getpgid(child.pid) == child.pid
                        assert os.getsid(child.pid) == child.pid
                        for number in (signal.SIGTERM, signal.SIGKILL):
                            try:
                                os.killpg(child.pid, number)
                            except ProcessLookupError:
                                pass
                            try:
                                child.wait(timeout=10)
                                break
                            except subprocess.TimeoutExpired:
                                pass
                        # Any surviving original group is still a failure; use
                        # the verified birth-based helper if it is observable.
                        if guard.group_members(child.pid):
                            guard.terminate_group(child, ticks, reaped)
                    except BaseException as fallback_error:
                        report['cleanup_errors'].append(repr(fallback_error))
            try:
                code = child.wait(timeout=10)
            except BaseException as error:
                report['cleanup_errors'].append(repr(error))
                code = child.returncode
        try:
            remaining = [] if child is None else guard.group_members(child.pid)
        except BaseException as terminal_error:
            remaining = None
            report['cleanup_errors'].append(repr(terminal_error))
        report.update(status='terminal', child_exit=code, wall_seconds=time.monotonic() - start,
                      stop_reason=reason or guard.STOP, peak_sampled_rss_kib=peak,
                      minimum_sampled_available_kib=minimum, adopted_children_reaped=reaped,
                      remaining_owned_processes=remaining)
        for pid, saved in sorted(observed.items()):
            try:
                current = guard.proc(pid)['start_ticks']
            except (FileNotFoundError, ProcessLookupError):
                current = None
            except BaseException as birth_error:
                report['cleanup_errors'].append(repr(birth_error))
                current = 'unobserved'
            report['observed_births'].append({'pid': pid, 'saved_start_ticks': saved,
                                             'current_start_ticks': current,
                                             'absent': current != 'unobserved' and current != saved})
        guard.write(folder / 'owner.json', report)

after = {str(p.relative_to(root)): guard.identity(p) for p in files}
guard.write(work / 'inputs-after.json', after)
report['source_unchanged'] = guard.source_identity() == source
report['whole_inputs_unchanged'] = before == after
guard.write(folder / 'owner.json', report)
assert code == 0 and report['stop_reason'] is None and not report['cleanup_errors']
assert report['source_unchanged'] and report['whole_inputs_unchanged'] and report['remaining_owned_processes'] == []
assert report['observed_births'] and all(e['absent'] for e in report['observed_births'])
ct = ET.parse(folder / 'ctest.xml').getroot()
assert [int(ct.get(k, '0')) for k in ('tests', 'failures', 'errors', 'skipped')] == [1, 0, 0, 0]
assert [e.get('name') for e in ct.findall('.//testcase')] == [plan['case']]
gx = ET.parse(folder / 'gtest.xml').getroot()
assert [int(gx.get(k, '0')) for k in ('tests', 'failures', 'errors', 'disabled')] == [1, 0, 0, 0]
cases = gx.findall('.//testcase')
assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
assert cases[0].get('classname') + '.' + cases[0].get('name') == plan['case'] and not cases[0].findall('skipped')
records = {e.get('name'): e.get('value') for e in cases[0].findall('./properties/property')}
names = re.findall(r'cases\.push_back\(\s*\{"([^"]+)"',
                   (root / 'tests/backend/full_path_acceptance_test.cpp').read_text().split('std::vector<Case> Cases() {', 1)[1].split('return cases;', 1)[0])
assert len(names) == plan['expected_original_witnesses']
assert set(records) == {name + '_r' + str(r) for name in names for r in range(plan['refinements'])}
report.update(failures=0, errors=0, skips=0, disabled=0, witness_records=records,
              ctest_xml=guard.identity(folder / 'ctest.xml'), gtest_xml=guard.identity(folder / 'gtest.xml'))
assert guard.STOP is None, ('controller interrupted before acceptance', guard.STOP)
report['passed'] = True
guard.write(folder / 'owner.json', report)
print(json.dumps({k: v for k, v in report.items() if k not in {'witness_records'}}, indent=2), flush=True)
