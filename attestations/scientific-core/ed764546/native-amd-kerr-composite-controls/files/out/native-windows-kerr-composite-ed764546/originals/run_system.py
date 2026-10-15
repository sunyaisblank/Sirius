"""Own four exact original native-AMD science controls using unchanged PowerShell ownership."""
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
ACTION_NAMES = {'retained-sampler':'RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates',
                'near-conservation':'FullPathAcceptance.VulkanRetainedNearExtremalRayConservesIndependentInvariants',
                'kerr-image':'KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow',
                'kerr-precision':'KernelTrace.Fp64RungAgreesWithFp32OnKerrScene'}
ORDER = list(ACTION_NAMES)
OWNER_SECONDS = {'retained-sampler':900, 'near-conservation':28800, 'kerr-image':28800, 'kerr-precision':28800}
kind, = sys.argv[1:]
assert kind in ACTION_NAMES
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['runner_source_revision']
SOURCE = plan['source_revision']
assert plan['admitted'] is True and binding['admitted'] is True
assert SOURCE == HEAD == '23d406b8ca99695f362b88e0c5f9a706f08a7a39'
assert plan['stage'] == binding['stage'] == 'bin/windows-msvc/diagnostic-23d406b'
assert plan['windows_checkout_unc']==binding['windows_checkout_unc']==r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'
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
    assert binding['manifest_identity'] == {'bytes':16928, 'sha256':'8553aebf9a13ec3f6109bb2f82ad63d57d80cf999f17fea49f54b5383eee2ce4'}
    assert binding['artifact']['id'] == 11957202404
    gate_path = stage / 'generated/sirius/native_build_gate.json'
    assert digest(gate_path) == binding['gate_sha256']
    assert binding['files'][str(gate_path.relative_to(ROOT))] == {'bytes':gate_path.stat().st_size,'sha256':binding['gate_sha256']}
    assert gate_path.stat().st_size == 10013 and binding['gate_sha256'] == 'b962047d6e47a5dc03eb986602e681fb2c064a0dc5d728d99b94acf278709e15'
    gate = json.loads(gate_path.read_text(encoding='utf-8-sig'))
    assert gate['source']['clean'] is True and gate['source']['revision'] == SOURCE
    for role, relative in [('sirius','src/sirius/app/Release/sirius.exe'), ('sirius_backend_tests','tests/backend/Release/sirius_backend_tests.exe')]:
        product = gate['tested_artifacts'][role]
        assert product['root'] == 'build' and product['path'] == relative
        assert {key:product[key] for key in ['bytes','sha256']} == manifest['files']['bin/windows-msvc/'+relative]
    return manifest

def verify_prerequisite():
    prerequisite = binding['prerequisite']
    assert prerequisite['admitted'] is True and prerequisite['commit'] == '6ef5828147c16395ca2bff4aeca6b662215b0918'
    assert prerequisite['tag'] == 'evidence/issue-48-native-windows-inventory-smoke-23d406b'
    assert prerequisite['bundle'] == 'attestations/native-vulkan/23d406b/inventory-and-original-smoke'
    assert prerequisite['work'] == prerequisite['bundle']+'/files/out/native-windows-diagnostic-23d406b'
    assert prerequisite['manifest_identity'] == {'bytes':15510,'sha256':'e9d8ea63d256769ed117a96785a9d829a1f8a73e6e0470b2fadece8133940fd2'}
    bundle = ROOT / prerequisite['bundle']
    archive_work = ROOT / prerequisite['work']
    assert bundle.resolve(strict=True) == bundle and archive_work.resolve(strict=True) == archive_work and archive_work.is_relative_to(bundle)
    manifest_path = bundle/'manifest.json'
    assert {'bytes':manifest_path.stat().st_size,'sha256':digest(manifest_path)} == prerequisite['manifest_identity']
    assert binding['files'][str(manifest_path.relative_to(ROOT))] == prerequisite['manifest_identity']
    archive_manifest = json.loads(manifest_path.read_text())['files']
    def archived(path, accepted):
        assert path.is_relative_to(archive_work) and path.resolve(strict=True) == path
        assert {'bytes':path.stat().st_size,'sha256':digest(path)} == accepted
        assert archive_manifest[str(path.relative_to(bundle))] == accepted
        assert binding['files'][str(path.relative_to(ROOT))] == accepted
    for previous in ['inventory','smoke']:
        path = archive_work/('terminal-'+previous+'-accepted.json')
        archived(path, prerequisite['accepted'][previous])
        accepted = json.loads(path.read_text());assert accepted['pass'] is True and accepted['case'] == previous
        owner = accepted['record']['owner'];bridge = accepted['record']['bridge']
        assert owner['source_revision'] == owner['live_source_revision'] == binding['source_revision']
        assert owner['source_build_gate_sha256'] == binding['gate_sha256'] and owner['source_build_gate_unchanged'] is True
        assert owner['status'] == 'completed' and owner['returncode'] == 0 and owner['owned_process_absent'] is True
        assert bridge['accepted_owner_execution'] is True and bridge['whole_inputs_unchanged'] is True
        for relative, pin in accepted['record']['whole_inputs'].items():
            logical=PurePosixPath(relative);assert not logical.is_absolute() and str(logical)==relative and '..' not in logical.parts
            archived(archive_work/logical, pin)
        for suffix, key in [('.json','terminal_receipt_sha256'),('-bridge.json','terminal_bridge_sha256')]:
            path=archive_work/('terminal-'+previous+suffix)
            assert digest(path) == accepted[key]
            archived(path, {'bytes':path.stat().st_size,'sha256':accepted[key]})
        if previous == 'smoke':assert owner['numerical_control_passed'] is True and owner['loaded_driver_modules_recorded'] is True
    for name in ['plan.json','execution-bindings.json','finite-result.json','independent-terminal-review.json','root-terminal-admission.json']:
        path=archive_work/name;archived(path, prerequisite['receipts'][name])
    old_plan=json.loads((archive_work/'plan.json').read_text())
    old_binding=json.loads((archive_work/'execution-bindings.json').read_text())
    assert old_plan['admitted'] is True and old_binding['admitted'] is True
    assert old_plan['source_revision']==old_plan['runner_source_revision']==old_binding['source_revision']==old_binding['runner_source_revision']==binding['source_revision']
    assert old_binding['manifest_identity']==binding['manifest_identity'] and old_binding['gate_sha256']==binding['gate_sha256']
    assert len(old_binding['files'])==95  # Historical declared map, never replayed against retired lexical paths.
    finite=json.loads((archive_work/'finite-result.json').read_text())
    independent=json.loads((archive_work/'independent-terminal-review.json').read_text())
    root=json.loads((archive_work/'root-terminal-admission.json').read_text())
    assert finite['schema']=='sirius-native-windows-inventory-original-smoke-result-v1' and finite['passed'] is True
    assert independent['schema']=='sirius-native-windows-inventory-smoke-independent-terminal-review-v1' and independent['passed'] is True
    assert root['schema']=='sirius-native-windows-root-terminal-admission-v1'
    assert all(row['source_revision']==binding['source_revision'] for row in [finite,independent,root])
    assert root['actual_root_commands_once_each_closed_zero'] is True and root['independent_terminal_accepted'] is True and root['finite_inventory_and_original_smoke_passed'] is True
    assert root['independent_terminal_review']==prerequisite['receipts']['independent-terminal-review.json']
    assert root['finite_result']==prerequisite['receipts']['finite-result.json']
    original_modules=json.loads((archive_work/'native/smoke/driver-modules.json').read_text(encoding='utf-8-sig'))
    assert [row for row in original_modules if 'amdvlk' in row['name'].lower()]==[binding['expected_amd_module']]
    inventory_path=archive_work/'native/inventory/stdout.log'
    inventory=json.loads(inventory_path.read_text(encoding='utf-8-sig'))
    assert inventory['platform']=='Windows' and inventory['wsl2'] is False
    vk=inventory['backends']['vulkan'];assert vk['device_count']==1 and vk['selected_device_index']==0
    assert vk['devices'][0]['name']=='AMD Radeon 780M Graphics'
    assert vk['devices'][0]['supports_fp64'] is True and vk['devices'][0]['rounds_fp64_to_nearest'] is True
    return inventory

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
    verify_prerequisite()

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
assert plan['order'] == ORDER and set(plan['actions']) == set(ORDER)
for label in ORDER:
    action = plan['actions'][label]
    assert set(action) == {'kind','test','expected_exit','owner_seconds','bridge_outer_seconds'}
    assert action['kind'] == 'Control' and action['test'] == ACTION_NAMES[label]
    assert type(action['expected_exit']) is int and action['expected_exit'] == 0
    assert type(action['owner_seconds']) is int and action['owner_seconds'] == OWNER_SECONDS[label]
    assert type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds'] == OWNER_SECONDS[label]+480
assert plan['owner_rss_bytes'] == 4294967296 and plan['terminal_observer_seconds'] == 60
bridge_seconds = plan['actions'][kind]['bridge_outer_seconds']
assert plan['minimum_wsl_available_kib'] == 2097152

verify()
stage = ROOT / plan['stage']
output = WORK / 'native' / kind
temporary = WORK / 'temp' / kind
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
assert not (WORK / (kind+'-inputs-before.json')).exists() and not (WORK / (kind+'-inputs-after.json')).exists()
exe = stage / 'tests/backend/Release/sirius_backend_tests.exe'
case = ACTION_NAMES[kind]
seconds = plan['actions'][kind]['owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'),
        '-Stage', win(exe.parent), '-Output', win(output), '-Filter', case, '-Seconds', str(seconds),
        '-RssBytes', str(plan['owner_rss_bytes']), '-TemporaryDirectory', win(temporary), '-BuildRoot', win(stage),
        '-ExpectedGateSha256', binding['gate_sha256'], '-LiveRevision', HEAD, '-ExecutableName', exe.name,
        '-InvocationKind', 'Control']
cleanup_argv = [ps, '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'emergency_native_cleanup.ps1'),
                '-TemporaryDirectory', win(temporary), '-Output', win(output), '-Receipt', win(WORK / (kind + '-emergency-cleanup.json'))]
assert argv == binding['actions'][kind]['argv'] and cleanup_argv == binding['actions'][kind]['cleanup_argv']
for previous in ORDER[:ORDER.index(kind)]:
    prior = json.loads((WORK/('terminal-'+previous+'-accepted.json')).read_text())
    assert prior['pass'] is True and prior['case'] == previous
    assert prior['record']['bridge']['binding_sha256'] == binding_hash
    assert prior['terminal_receipt_sha256'] == digest(WORK/('terminal-'+previous+'.json'))
    assert prior['terminal_bridge_sha256'] == digest(WORK/('terminal-'+previous+'-bridge.json'))
    prior_owner = prior['record']['owner']
    assert prior_owner == json.loads((WORK/'native'/previous/'owner.json').read_text(encoding='utf-8-sig'))
    assert prior_owner['source_revision'] == SOURCE and prior_owner['numerical_control_passed'] is True and prior_owner['loaded_driver_modules_recorded'] is True
    for relative, accepted in prior['record']['whole_inputs'].items():
        logical=PurePosixPath(relative);assert not logical.is_absolute() and str(logical)==relative and '..' not in logical.parts
        path = WORK/logical
        assert {'bytes':path.stat().st_size, 'sha256':digest(path)} == accepted, relative
        predecessors[str(path.relative_to(ROOT))] = accepted
    for relative in ['terminal-'+previous+'.json','terminal-'+previous+'-bridge.json','terminal-'+previous+'-accepted.json']:
        path=WORK/relative;predecessors[str(path.relative_to(ROOT))]={'bytes':path.stat().st_size,'sha256':digest(path)}
    verify()
before = input_records()
write(WORK / (kind+'-inputs-before.json'), before)
record = {'source_revision': SOURCE, 'runner_source_revision': HEAD, 'kind': kind, 'argv': argv, 'binding_sha256': binding_hash,
          'linux_controller_pid': os.getpid(), 'linux_controller_birth': birth(os.getpid()),
          'predecessor_inputs': predecessors,
          'bridge_outer_seconds': bridge_seconds, 'windows_owner_seconds': seconds, 'minimum_wsl_available_kib': plan['minimum_wsl_available_kib'],
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
            assert time.monotonic() - started < bridge_seconds, 'bridge outer diagnostic stop'
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
    record['native_build_gate_authority_used'] = True
    record['mandatory_or_full_runtime_admitted'] = False
    assert owner['numerical_control_passed'] is True and owner['system_inventory_passed'] is False and owner['loaded_driver_modules_recorded'] is True
    xml=ET.parse(output/'gtest.xml').getroot();cases=xml.findall('.//testcase')
    assert int(xml.get('tests','-1')) == 1 and len(cases) == 1
    assert all(int(node.get(key,'0')) == 0 for node in [xml,*xml.findall('.//testsuite')] for key in ['failures','errors','skipped','disabled'])
    assert not any(xml.findall('.//'+key) for key in ['failure','error','skipped'])
    assert cases[0].get('classname')+'.'+cases[0].get('name') == case and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    properties={p.get('name'):p.get('value') for p in cases[0].findall('./properties/property')}
    if kind=='kerr-precision':assert properties.get('fp64_evidence')=='native_comparison_executed'
    if kind=='near-conservation':
        modes=int(properties['executed_product_modes']);assert modes in (1,2)
        record['executed_product_modes']=modes
        record['binary64_product_mode_executed']=modes==2
    record['loaded_provider_evidence_complete'] = owner['loaded_driver_modules_recorded'] is True
    if owner['loaded_driver_modules_recorded']:
        modules=json.loads((output/'driver-modules.json').read_text(encoding='utf-8-sig'))
        assert isinstance(modules,list) and len(modules) == 2
        loader=[row for row in modules if row['name'].lower() == 'vulkan-1.dll']
        assert len(loader) == 1 and any('amdvlk' in row['name'].lower() for row in modules)
        assert [row for row in modules if 'amdvlk' in row['name'].lower()] == [binding['expected_amd_module']]
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
