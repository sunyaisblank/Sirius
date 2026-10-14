from pathlib import Path
import ctypes, importlib.util, json, os, signal, subprocess, sys, time
sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
settings = json.loads((work / 'actions.json').read_text())
action = settings[sys.argv[1]]
plan = {'case': sys.argv[1], 'scope': action['scope'], 'outer_limits': action['limits']}
spec = importlib.util.spec_from_file_location('preserved_owned_guard', work / 'owned_guard.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)
guard.ROOT = root
assert ctypes.CDLL(None, use_errno=True).prctl(36, 1, 0, 0, 0) == 0
for number in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
    signal.signal(number, guard.signal_stop)
guard.RSS_KIB = action['limits']['rss_kib']
guard.RESERVE_KIB = action['limits']['reserve_kib']
# Ninja may create compiler process groups inside the same owned session.
# Observe and terminate the complete session rather than just its first group.
def owned_session_members(session):
    members = []
    for entry in Path('/proc').glob('[0-9]*'):
        try:
            item = guard.proc(int(entry.name))
            if item['session'] == session:
                members.append(item)
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return members

def terminate_owned_session(child, start_ticks, reaped):
    for number, grace in ((signal.SIGTERM, 10), (signal.SIGKILL, 10)):
        members = owned_session_members(child.pid)
        assert not any(e['pid'] == child.pid and e['start_ticks'] != start_ticks for e in members)
        for item in members:
            try:
                descriptor = os.pidfd_open(item['pid'])
                try:
                    current = guard.proc(item['pid'])
                    if current['start_ticks'] == item['start_ticks'] and current['session'] == child.pid:
                        signal.pidfd_send_signal(descriptor, number)
                finally:
                    os.close(descriptor)
            except (FileNotFoundError, ProcessLookupError):
                pass
        deadline = time.monotonic() + grace
        while time.monotonic() < deadline:
            if child.poll() is not None:
                # The subreaper owns adopted descendants, including other groups.
                for item in owned_session_members(child.pid):
                    try:
                        pid, status = os.waitpid(item['pid'], os.WNOHANG)
                        if pid:
                            reaped.append({'pid': pid, 'wait_status': status})
                    except ChildProcessError:
                        pass
            if child.poll() is not None and not owned_session_members(child.pid):
                return
            time.sleep(.1)
    raise RuntimeError('owned process session did not terminate within bounded cleanup grace')

guard.group_members = owned_session_members
guard.terminate_group = terminate_owned_session
source = guard.source_identity()
assert source == {'revision': action['revision'], 'status': ''}
folder = work / sys.argv[1]
folder.mkdir()
env = os.environ.copy()
allowed = ({'VK_DRIVER_FILES':'/usr/share/vulkan/icd.d/lvp_icd.json', 'SIRIUS_VULKAN_DEVICE':'0', 'MESA_SHADER_CACHE_DISABLE':'true'} if action.get('software_provider') else {})
overrides = [key for key, value in env.items() if value and allowed.get(key) != value and
             (key.startswith(('SIRIUS_', 'VK_', 'LP_', 'MESA_', 'GALLIVM_', 'GTEST_', 'CTEST_'))
              or key in {'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'})]
assert not overrides, sorted(overrides)
env.update(allowed)
if action['kind'] == 'test':
    env['GTEST_OUTPUT'] = 'xml:' + str(folder / 'gtest.xml')
else:
    env['GTEST_OUTPUT'] = ''
argv = [arg.replace('{folder}', str(folder)) for arg in action['argv']]
tracked = [root / os.fsdecode(p) for p in guard.query(['git', 'ls-files', '-z']).split(b'\0') if p]
files = [*tracked, Path(__file__).resolve(), work / 'owned_guard.py', work / 'actions.json']
files += [root / p for p in action.get('frozen_inputs', [])]
key = lambda p: str(p.relative_to(root)) if p.is_relative_to(root) else str(p)
before = {key(p): guard.identity(p) for p in files}
guard.write(folder / 'inputs-before.json', before)
report = {'source': source, 'case': plan['case'], 'argv': argv,
          'scope': plan['scope'], 'limits': plan['outer_limits'],
          'declared_output_environment': {'GTEST_OUTPUT': env['GTEST_OUTPUT'], **allowed}, 'provider_mappings': [],
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
        assert guard.STOP is None, ('interrupted before spawn', guard.STOP)
        child = subprocess.Popen(argv, cwd=root, env=env, stdout=out, stderr=err, start_new_session=True)
        ticks = guard.proc(child.pid)['start_ticks']
        observed[child.pid] = ticks
        report.update(child_pid=child.pid, child_start_ticks=ticks, status='running')
        guard.write(folder / 'owner.json', report)
        while child.poll() is None:
            members = guard.group_members(child.pid)
            for item in members:
                observed[item['pid']] = item['start_ticks']
            if allowed and not report['provider_mappings']:
                for item in members:
                    try:
                        prior = guard.proc(item['pid'])
                        assert prior['start_ticks'] == item['start_ticks'] and prior['session'] == child.pid
                        provider = '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so'
                        lines = [line for line in Path(f"/proc/{item['pid']}/maps").read_text().splitlines() if provider in line]
                        current = guard.proc(item['pid'])
                        assert current['start_ticks'] == item['start_ticks'] and current['session'] == child.pid
                        if lines:
                            report['provider_mappings'].append({'pid': item['pid'], 'start_ticks': item['start_ticks'], 'owned_session': child.pid, 'maps': lines})
                    except (FileNotFoundError, ProcessLookupError):
                        pass
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
            time.sleep(.02 if allowed else .25)
    except BaseException as error:
        reason = reason or 'observer_exception'
        report['observer_error'] = repr(error)
    finally:
        if child is not None:
            try:
                if child.poll() is not None:
                    for item in guard.group_members(child.pid):
                        if item['state'] == 'Z':
                            try:
                                pid, status = os.waitpid(item['pid'], os.WNOHANG)
                                if pid:
                                    reaped.append({'pid': pid, 'wait_status': status})
                            except ChildProcessError:
                                pass
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

after = {key(p): guard.identity(p) for p in files}
guard.write(folder / 'inputs-after.json', after)
report['source_unchanged'] = guard.source_identity() == source
report['whole_inputs_unchanged'] = before == after
guard.write(folder / 'owner.json', report)
assert report['stop_reason'] is None and not report['cleanup_errors'], report
assert report['source_unchanged'] and report['whole_inputs_unchanged']
assert report['remaining_owned_processes'] == []
assert not allowed or report['provider_mappings'], 'Actual selected software provider mapping was not observed'
assert report['observed_births'] and all(e['absent'] for e in report['observed_births'])
assert code == action['expected_exit'], (code, action['expected_exit'])
assert guard.STOP is None, ('interrupted before command acceptance', guard.STOP)
report['passed'] = True
report['meaning'] = action['scope']
guard.write(folder / 'owner.json', report)
print(json.dumps({k: report[k] for k in ('case', 'source', 'passed', 'child_exit', 'wall_seconds', 'peak_sampled_rss_kib', 'stop_reason', 'meaning')}, indent=2), flush=True)
