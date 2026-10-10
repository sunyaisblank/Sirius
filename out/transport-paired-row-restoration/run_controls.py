from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os, signal, subprocess, time, xml.etree.ElementTree as ET

r = Path.cwd(); w = Path(__file__).resolve().parent
head = json.loads((w/'source.json').read_text())['source_revision']
def git(*args): return subprocess.check_output(['git', *args], text=True).strip()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, value): p.write_text(json.dumps(value, indent=2) + '\n')
def birth(pid): return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
def memory(pid):
    pending = [pid]; seen = set(); rss = 0
    while pending:
        current = pending.pop()
        if current in seen: continue
        seen.add(current)
        try:
            folder = Path('/proc', str(current))
            status = dict(line.split(':', 1) for line in (folder/'status').read_text().splitlines() if ':' in line)
            rss += int(status.get('VmRSS', '0 kB').split()[0])
            pending.extend(map(int, (folder/'task'/str(current)/'children').read_text().split()))
        except (FileNotFoundError, ProcessLookupError): pass
    available = next(int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))
    return available, rss

def group_members(pgid):
    members = []
    for folder in Path('/proc').iterdir():
        if not folder.name.isdecimal(): continue
        try:
            fields = (folder/'stat').read_text().rsplit(')', 1)[1].split()
        except (FileNotFoundError, ProcessLookupError): continue
        if int(fields[2]) == pgid and int(fields[3]) == pgid:
            members.append({'pid':int(folder.name), 'start_ticks':fields[19]})
    return sorted(members, key=lambda entry:entry['pid'])
def signal_group(pgid, action):
    try: os.killpg(pgid, action)
    except ProcessLookupError: pass

controller_signals = []
def interrupted(signum, frame):
    controller_signals.append(signum)
signal.signal(signal.SIGTERM, interrupted)
signal.signal(signal.SIGINT, interrupted)

assert git('rev-parse', 'HEAD') == head and not git('status', '--porcelain')
gtest_overrides = {key:value for key,value in os.environ.items() if key.startswith('GTEST_') and value}
assert not gtest_overrides, ('unexpected GoogleTest environment overrides', sorted(gtest_overrides))
provider = [Path('/usr/share/vulkan/icd.d/lvp_icd.json'), Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so'), Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve()]
selected_environment = {'VK_DRIVER_FILES': str(provider[0]), 'SIRIUS_VULKAN_DEVICE': '', 'MESA_SHADER_CACHE_DISABLE': '', 'LP_NATIVE_VECTOR_WIDTH': '', 'LP_PERF': '', 'GALLIVM_PERF': ''}
backend = r/'bin/linux-gcc/tests/backend/sirius_backend_tests'
render = r/'bin/linux-gcc/src/sirius/app/sirius_render_tests'
cases = [
 ('arithmetic-selection', backend, 'RetainedComputeAdmission.FmaSelectsOnlyNativeWideProductsAndPreservesAllocation', 30, False),
 ('preparation-model', backend, 'RetainedComputeAdmission.SoftwareRendererPreparationPreservesPhysicalAccounting', 30, False),
 ('critical-transport', backend, 'RetainedComputeTest.JointRkStagesRetainCriticalIncrementsAndEmbeddedError', 90, True),
 ('schwarzschild', backend, 'RetainedComputeTest.SchwarzschildStagesPreserveIndependentFieldsAndGeneralFallback', 90, True),
 ('mixed-transport', backend, 'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal', 180, True),
 ('rejected-step', backend, 'RetainedComputeTest.RejectedStepRowsCannotExposeOldOrPartialCandidates', 90, True),
 ('coupled', backend, 'RetainedComputeTest.CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement', 90, True),
 ('shared-tracer', backend, 'RetainedComputeTest.SharedTracerCompletesDeviceIntervalsAndRetainsRollbackState', 90, True),
]
build = json.loads((w/'build-and-arrays.json').read_text())
assert build['source_revision'] == head and build['returncode'] == 0 and build['restored_all_40_arrays_exact']
payloads = json.loads((w/'linux-readonly-payloads.json').read_text())
assert payloads['source_revision'] == head and payloads['pass_']
assert sha(r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h') == build['whole_header_sha256'] == payloads['whole_header_sha256']
for relative, entry in payloads['actual_consumers'].items():
    artifact = r/relative
    assert artifact.stat().st_size == entry['bytes'] and sha(artifact) == entry['sha256']
    assert set(entry['arrays']) == set(build['arrays'])
    for name, record in entry['arrays'].items():
        assert {key: record[key] for key in ('bytes', 'sha256')} == build['arrays'][name]
        assert record['occurrences']
gate = r/'bin/linux-gcc/generated/sirius/alignment_receipt.json'
assert json.loads(gate.read_text())['source_revision'] == head
files = [Path(__file__).resolve(), w/'source.json', w/'build-and-arrays.json', w/'linux-readonly-payloads.json', gate,
 r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
 r/'bin/linux-gcc/tests/backend/retained_camera/program_fixture.h', r/'bin/linux-gcc/tests/backend/portable_binary32_reference.bin',
 *provider, *[r/p for p in payloads['actual_consumers']],
 *[r/p for p in json.loads((w/'source.json').read_text())['changed_files']],
 *sorted((r/'src/sirius/kernels').glob('retained*')),
 *sorted((r/'tests/support/retained_transport').glob('*.h')),
 r/'src/sirius/backend/retained_integrator.cpp', r/'src/sirius/backend/retained_integrator.h',
 r/'src/sirius/backend/retained_trace_executor.cpp', r/'src/sirius/backend/retained_trace_executor.h',
 r/'tests/backend/retained_dopri_test.cpp', r/'scripts/build-retained-kernels.py', r/'tests/operating_model.json', r/'scripts/verify-operating-model.py']
whole_seals = {str(p): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in files}
receipts = []
for name, executable, case, timeout, needs_provider in cases:
    maximum_rss = 4096
    folder = w/'controls'/name
    assert not folder.exists(); folder.mkdir(parents=True)
    seals = whole_seals
    env = os.environ.copy()
    for key, value in selected_environment.items(): env[key] = value
    args = [str(executable), '--gtest_filter='+case, '--gtest_output=xml:'+str(folder/'tests.xml')]
    owner = {'source_revision': head, 'source_tree': git('rev-parse', 'HEAD^{tree}'), 'case': case,
             'command': args, 'whole_input_seals': seals, 'controller_pid': os.getpid(), 'controller_start_ticks': birth(os.getpid()),
             'timeout_seconds': timeout, 'maximum_rss_mib': maximum_rss, 'minimum_available_mib': None,
             'sample_cadence_seconds': 0.25, 'peak_sampled_rss_kib': 0, 'minimum_sampled_available_kib': None,
             'selected_environment': selected_environment, 'gtest_environment_overrides': gtest_overrides,
             'scope': 'Restored original arithmetic/preparation/independent Transport/refusal/coupled/whole-state rollback controls after rejected paired-row trial. Original existing cases/settings/assertions and 30/90/180s/4096MiB guards retained. Selected llvmpipe route recorded by executed controls; no native hardware, performance/adoption/frame/full scientific or release qualification claim.'}
    started = time.monotonic(); loaded = set(); stop = None
    owner_error, cleanup_errors = None, []
    observed_births = {}
    assert not controller_signals, ('controller interrupted before spawn', controller_signals)
    with (folder/'stdout.log').open('wb') as out, (folder/'stderr.log').open('wb') as err:
        child = subprocess.Popen(args, cwd=r, env=env, stdout=out, stderr=err, start_new_session=True)
        try:
            if controller_signals:
                stop = 'controller_signal'
            owner.update({'child_pid': child.pid, 'child_start_ticks': birth(child.pid), 'started_utc': datetime.now(timezone.utc).isoformat()})
            observed_births[child.pid] = owner['child_start_ticks']
            dump(folder/'owner.json', owner)
            while child.poll() is None:
                if controller_signals:
                    stop = 'controller_signal'
                    break
                for member in group_members(child.pid):
                    observed_births[member['pid']] = member['start_ticks']
                available, rss = memory(child.pid)
                owner['peak_sampled_rss_kib'] = max(owner['peak_sampled_rss_kib'], rss)
                owner['minimum_sampled_available_kib'] = min(owner['minimum_sampled_available_kib'] or available, available)
                try:
                    for line in Path('/proc', str(child.pid), 'maps').read_text().splitlines():
                        path = line.split()[-1]
                        if path.startswith('/') and ('libvulkan_lvp.so' in path or 'libvulkan.so' in path): loaded.add(str(Path(path).resolve()))
                except (FileNotFoundError, ProcessLookupError): pass
                if time.monotonic()-started >= timeout: stop = 'case_time_bound'
                elif rss > maximum_rss*1024: stop = 'case_sampled_rss_guard'
                dump(folder/'owner.json', owner)
                if stop: break
                time.sleep(0.25)
        except BaseException as error:
            owner_error = repr(error)
            stop = stop or 'observer_exception'
        finally:
            try:
                remaining = group_members(child.pid)
                if remaining:
                    stop = stop or 'owned_group_remaining'
                    signal_group(child.pid, signal.SIGTERM)
                    deadline = time.monotonic()+10
                    while group_members(child.pid) and time.monotonic() < deadline:
                        child.poll()
                        time.sleep(0.1)
                    if group_members(child.pid): signal_group(child.pid, signal.SIGKILL)
                code = child.wait(timeout=10)
            except BaseException as error:
                cleanup_errors.append(repr(error))
                try:
                    signal_group(child.pid, signal.SIGKILL)
                    code = child.wait(timeout=10)
                except BaseException as last_error:
                    cleanup_errors.append(repr(last_error))
                    code = child.returncode
            try:
                remaining = group_members(child.pid)
            except BaseException as error:
                cleanup_errors.append(repr(error))
                remaining = None
    birth_checks = []
    for pid, saved_birth in sorted(observed_births.items()):
        try: actual_birth = birth(pid)
        except (FileNotFoundError, ProcessLookupError): actual_birth = None
        birth_checks.append({'pid':pid, 'saved_start_ticks':saved_birth, 'current_start_ticks':actual_birth, 'absent':actual_birth != saved_birth})
    if controller_signals:
        stop = 'controller_signal'
    owner.update({'controller_signals': list(controller_signals), 'returncode': code, 'elapsed_seconds': time.monotonic()-started, 'stop_reason': stop, 'sampled_loaded_provider_paths': sorted(loaded),
                  'observer_error':owner_error, 'cleanup_errors':cleanup_errors, 'saved_birth_checks':birth_checks,
                  'remaining_owned_group':remaining, 'owned_birth_absent': bool(birth_checks) and all(entry['absent'] for entry in birth_checks),
                  'owned_group_absent':remaining == [], 'terminal_visibility':'Linux procfs process group/session scan; accessible stats, no global consumer claim.'})
    owner['source_and_inputs_unchanged'] = git('rev-parse', 'HEAD') == head and not git('status', '--porcelain') and all(p.stat().st_size == x['bytes'] and sha(p) == x['sha256'] for p, x in ((Path(k), v) for k, v in seals.items()))
    dump(folder/'owner.json', owner)
    assert code == 0 and stop is None and owner_error is None and not cleanup_errors and owner['owned_birth_absent'] and owner['owned_group_absent'] and owner['source_and_inputs_unchanged'], name
    root = ET.parse(folder/'tests.xml').getroot()
    assert [int(root.get(k, '0')) for k in ['tests', 'failures', 'errors', 'disabled']] == [1, 0, 0, 0], name
    tests = root.findall('.//testcase'); assert len(tests) == 1 and tests[0].get('status') == 'run' and tests[0].get('result') == 'completed' and not tests[0].findall('skipped'), name
    properties = {p.get('name'): p.get('value') for p in tests[0].findall('./properties/property')}
    assert tests[0].get('classname')+'.'+tests[0].get('name') == case
    assert not needs_provider or {str(p.resolve()) for p in provider[1:]} <= loaded
    result = {'name': name, 'case': case, 'source_revision': head, 'returncode': code, 'failures': 0, 'errors': 0, 'skips': 0,
              'elapsed_seconds': owner['elapsed_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'], 'properties': properties,
              'owner_sha256': sha(folder/'owner.json'), 'xml_sha256': sha(folder/'tests.xml')}
    dump(folder/'result.json', result); receipts.append(result); dump(w/'controls-sequence.json', receipts)
    print(json.dumps({key: value for key, value in result.items() if key != 'properties'}), flush=True)
dump(w/'controls-complete.json', {'source_revision': head, 'returncode': 0, 'existing_tests': len(cases), 'new_test_registrations': 0, 'results': receipts, 'full_qualification': False})
