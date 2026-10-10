from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os, signal, subprocess, time, xml.etree.ElementTree as ET

r = Path.cwd(); w = Path(__file__).resolve().parent
head = json.loads((w/'test-first-source.json').read_text())['source_revision']
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

assert git('rev-parse', 'HEAD') == head and not git('status', '--porcelain')
build = json.loads((w/'baseline-build-and-arrays.json').read_text())
assert build['source_revision'] == head and build['returncode'] == 0 and build['all_40_arrays_exact']
header = r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert sha(header) == build['whole_header_sha256']
payloads = json.loads((w/'baseline-linux-readonly-payloads.json').read_text())
assert payloads['source_revision'] == head and payloads['whole_header_sha256'] == sha(header) and payloads['pass_']
assert len(build['arrays']) == 40
for relative, record in payloads['actual_consumers'].items():
    artifact = r/relative
    assert artifact.stat().st_size == record['bytes'] and sha(artifact) == record['sha256']
    assert set(record['arrays']) == set(build['arrays'])
    for name, entry in record['arrays'].items():
        assert {key:entry[key] for key in ('bytes','sha256')} == build['arrays'][name]
        assert entry['occurrences']
gate = r/'bin/linux-gcc/generated/sirius/alignment_receipt.json'
assert json.loads(gate.read_text())['source_revision'] == head
old_identity = json.loads((r/'attestations/software-vulkan/3045067/dopri-normal-guard-restoration/matched/restoration/llvmpipe/timestamps/identity.json').read_text())
old_owner = json.loads((r/'attestations/software-vulkan/3045067/dopri-normal-guard-restoration/matched/restoration/llvmpipe/timestamps/owner.json').read_text())
assert old_owner['outer_guard_seconds'] == 90 and old_owner['rss_guard_bytes'] == 4*1024**3
old = {'selected_environment': old_identity['environment']}
review = json.loads((w/'mixed-test-source-review.json').read_text())
assert review['pass']
provider = [Path('/usr/share/vulkan/icd.d/lvp_icd.json'), Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so'), Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve()]
provider_expected = {entry['path']:entry for manifest in old_identity['provider_inputs']['manifests'] for entry in [manifest['manifest'], *manifest['driver_candidates']]}
provider_expected.update({entry['path']:entry for entry in old_identity['provider_inputs']['loader_candidates']})
for path in provider:
    expected = provider_expected[str(path)]
    assert path.stat().st_size == expected['bytes'] and sha(path) == expected['sha256']
for path, expected in review['seals'].items():
    assert Path(path).stat().st_size == expected['bytes'] and sha(Path(path)) == expected['sha256']
changed = git('diff-tree', '--no-commit-id', '--name-only', '-r', 'HEAD').splitlines()
assert changed == ['scripts/verify-operating-model.py', 'tests/backend/retained_compute_test.cpp', 'tests/labels/CTestLabels.cmake', 'tests/operating_model.json']
cases = [('timestamps', r/'bin/linux-gcc/tests/backend/sirius_backend_tests', 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 90, 4096)]
receipts = []
for name, executable, case, timeout, maximum_rss in cases:
    folder = w/'baseline-timestamp-controls'/name
    assert not folder.exists(); folder.mkdir(parents=True)
    files = [executable, gate, r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h', *[r/p for p in changed]]
    if True:
        files += provider
        files += [r/'bin/linux-gcc/tests/backend/portable_binary32_reference.bin', r/'tests/support/retained_transport/reference_cases.h', r/'bin/linux-gcc/tests/backend/retained_camera/program_fixture.h', Path(__file__).resolve(), w/'mixed-test-source-review.json', w/'test-first-source.json', w/'baseline-linux-readonly-payloads.json', w/'baseline-timestamp-owner-review.json', r/'attestations/software-vulkan/3045067/dopri-normal-guard-restoration/matched/restoration/llvmpipe/timestamps/gtest.xml']
        files += sorted((r/'src/sirius/kernels').glob('retained*'))
        files += [r/'scripts/build-retained-kernels.py', r/'src/sirius/backend/retained_compute.cpp', r/'src/sirius/backend/retained_compute.h']
    seals = {str(p): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in files}
    env = os.environ.copy()
    if True:
        # Preserve the original software case's provider and cache settings.
        for key, value in old['selected_environment'].items():
            if value is None: env.pop(key, None)
            else: env[key] = value
    args = [str(executable), '--gtest_filter='+case, '--gtest_output=xml:'+str(folder/'tests.xml')]
    owner = {'source_revision': head, 'source_tree': git('rev-parse', 'HEAD^{tree}'), 'case': case,
             'command': args, 'whole_input_seals': seals, 'controller_pid': os.getpid(), 'controller_start_ticks': birth(os.getpid()),
             'timeout_seconds': timeout, 'maximum_rss_mib': maximum_rss, 'minimum_available_mib': None,
             'sample_cadence_seconds': 0.25, 'peak_sampled_rss_kib': 0, 'minimum_sampled_available_kib': None,
             'selected_environment': {key: env.get(key) for key in old['selected_environment']},
             'scope': 'One original timestamp control on accepted production at test-first revision; original90s/4GiB observer limits and numerical fixtures/modes/budgets preserved. Complete same-provider32 readbacks provide fresh software baseline. No latency retention, native, full frame or full scientific/release qualification claim.'}
    started = time.monotonic(); loaded = set(); stop = None
    owner_error, cleanup_errors = None, []
    observed_births = {}
    with (folder/'stdout.log').open('wb') as out, (folder/'stderr.log').open('wb') as err:
        child = subprocess.Popen(args, cwd=r, env=env, stdout=out, stderr=err, start_new_session=True)
        try:
            owner.update({'child_pid': child.pid, 'child_start_ticks': birth(child.pid), 'started_utc': datetime.now(timezone.utc).isoformat()})
            observed_births[child.pid] = owner['child_start_ticks']
            dump(folder/'owner.json', owner)
            while child.poll() is None:
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
    owner.update({'returncode': code, 'elapsed_seconds': time.monotonic()-started, 'stop_reason': stop, 'sampled_loaded_provider_paths': sorted(loaded),
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
    assert {str(p.resolve()) for p in provider[1:]} <= loaded
    accepted_xml = r/'attestations/software-vulkan/3045067/dopri-normal-guard-restoration/matched/restoration/llvmpipe/timestamps/gtest.xml'
    accepted = {p.get('name'):p.get('value') for p in ET.parse(accepted_xml).findall('.//property')}
    keys = [key for key in properties if '_readback_' in key]
    assert len(keys)==32 and all(properties[key]==accepted[key] for key in keys)
    assert len([key for key in properties if key.startswith('device_observation_')])==51
    result = {'name': name, 'case': case, 'source_revision': head, 'returncode': code, 'failures': 0, 'errors': 0, 'skips': 0,
              'elapsed_seconds': owner['elapsed_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'], 'properties': properties,
              'owner_sha256': sha(folder/'owner.json'), 'xml_sha256': sha(folder/'tests.xml')}
    dump(folder/'result.json', result); receipts.append(result); dump(w/'baseline-timestamp-controls-sequence.json', receipts)
    print(json.dumps({key: value for key, value in result.items() if key != 'properties'}), flush=True)
dump(w/'baseline-timestamp-controls-complete.json', {'source_revision': head, 'returncode': 0, 'existing_tests': 1, 'new_test_registrations': 0, 'results': receipts, 'full_qualification': False})
