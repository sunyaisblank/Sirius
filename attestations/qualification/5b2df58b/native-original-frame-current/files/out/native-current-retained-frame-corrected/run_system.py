"""Own one original current native retained Kerr frame; no shader override."""
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
ACTION_NAMES = {
    'native-frame':'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'}
PROFILES = {'original-retained-frame':['native-frame']}
OWNER_SECONDS = {name:300 for name in ACTION_NAMES}
EVIDENCE_PROPERTIES = {}
kind, = sys.argv[1:]
assert kind in ACTION_NAMES
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['runner_source_revision']
SOURCE = plan['source_revision']
assert plan['schema'] == 'sirius-native-windows-current-retained-frame-plan-v1'
assert binding['schema'] == 'sirius-native-windows-current-retained-frame-bindings-v1'
assert plan['control_profile'] == binding['control_profile'] and plan['control_profile'] in PROFILES
ORDER = PROFILES[plan['control_profile']]
assert kind in ORDER and set(binding['actions']) == set(ORDER) and set(binding['terminal_actions']) == set(ORDER)
assert plan['admitted'] is True and binding['admitted'] is True
assert SOURCE == HEAD == '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'
assert plan['stage'] == binding['stage'] == 'bin/windows-msvc/current-original-frame-5b2df58b'
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
                       'repository':'sunyaisblank/Sirius', 'run_id':'39089567447', 'run_attempt':'1'}.items():
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
    gate_path = stage / 'generated/sirius/native_build_gate.json'
    assert digest(gate_path) == binding['gate_sha256']
    assert binding['files'][str(gate_path.relative_to(ROOT))] == {'bytes':gate_path.stat().st_size,'sha256':binding['gate_sha256']}
    gate_identity = binding['native_gate_identity']
    assert set(gate_identity) == {'bytes','sha256'}
    assert type(gate_identity['bytes']) is int and gate_identity['bytes'] > 0
    assert re.fullmatch('[0-9a-f]{64}', gate_identity['sha256'])
    assert {'bytes':gate_path.stat().st_size,'sha256':binding['gate_sha256']} == gate_identity
    gate = json.loads(gate_path.read_text(encoding='utf-8-sig'))
    assert gate['source']['clean'] is True and gate['source']['revision'] == SOURCE
    for role, relative in [('sirius','src/sirius/app/Release/sirius.exe'), ('sirius_backend_tests','tests/backend/Release/sirius_backend_tests.exe')]:
        product = gate['tested_artifacts'][role]
        assert product['root'] == 'build' and product['path'] == relative
        assert {key:product[key] for key in ['bytes','sha256']} == manifest['files']['bin/windows-msvc/'+relative]
    # Exactly one extra PE, independently bound by the unchanged native artifact manifest/Gate.
    assert len(manifest['files']) == 138
    extra_relative = 'src/sirius/app/Release/sirius_render_tests.exe'
    assert 'bin/windows-msvc/'+extra_relative not in manifest['files']
    assert gate['kind'] == 'sirius-native-build-gate' and type(gate['schema_version']) is int and gate['schema_version'] == 2 and gate['status'] == 'passed'
    assert type(gate['ctest']['executed']) is int and gate['ctest']['executed'] == 9
    assert all(type(gate['ctest'][key]) is int and gate['ctest'][key] == 0 for key in ['failures','errors','skipped'])
    extra = gate['tested_artifacts']['sirius_render_tests']
    assert set(extra) == {'root','path','bytes','sha256'} and extra['root'] == 'build' and extra['path'] == extra_relative
    extra_identity = {key:extra[key] for key in ['bytes','sha256']}
    assert type(extra_identity['bytes']) is int and extra_identity['bytes'] > 0 and re.fullmatch('[0-9a-f]{64}', extra_identity['sha256'])
    assert binding['additional_test_executable'] == {'role':'sirius_render_tests','stage_relative_path':extra_relative,'identity':extra_identity}
    extra_path = stage / extra_relative
    assert extra_path.resolve(strict=True) == extra_path and extra_path.is_file()
    assert {'bytes':extra_path.stat().st_size,'sha256':digest(extra_path)} == extra_identity == binding['files'][str(extra_path.relative_to(ROOT))]
    native_path = ROOT / binding['native_attestation_path']
    assert native_path.resolve(strict=True) == native_path and native_path.is_relative_to(ROOT)
    assert {'bytes':native_path.stat().st_size,'sha256':digest(native_path)} == binding['native_attestation_identity'] == binding['files'][str(native_path.relative_to(ROOT))]
    native = json.loads(native_path.read_text(encoding='utf-8-sig'))
    assert type(native['schema_version']) is int and native['schema_version'] == 1 and native['status'] == 'pass' and native['domains'] == ['windows-native-build']
    assert native['source_revision'] == manifest['source_revision'] and native['platform']['os'] == 'windows' and native['platform']['wsl2'] is False
    assert native['claims']['runtime_ready'] is False and native['claims']['test_estate_kind'] == 'native-build-non-render-authority'
    assert native['claims']['test_estate_passed'] is True and native['claims']['test_report'] == {'cases':9,'failures':0,'errors':0,'skipped':0}
    assert all(type(value) is int for value in native['claims']['test_report'].values())
    for key, identity in [('native_build_gate.json',gate_identity),('native-build-tested-sirius_render_tests',extra_identity)]:
        assert native['artifacts'][key] == dict(path=key, **identity)
        raw = native_path.parent / key
        assert raw.resolve(strict=True) == raw and {'bytes':raw.stat().st_size,'sha256':digest(raw)} == identity == binding['files'][str(raw.relative_to(ROOT))]
    assert binding['native_artifact']['name'] == 'sirius-windows-native-attestation' and binding['native_artifact']['repository'] == manifest['repository']
    assert binding['native_artifact']['run_id'] == manifest['run_id'] and binding['native_artifact']['run_attempt'] == manifest['run_attempt']
    assert type(binding['native_artifact']['id']) is int and binding['native_artifact']['id'] > 0
    expected_stage = {str(stage/PurePosixPath(key).relative_to('bin/windows-msvc')) for key in manifest['files']}
    expected_stage |= {str(manifest_path),str(gate_path),str(extra_path)}
    assert {str(p) for p in stage.rglob('*') if p.is_file()} == expected_stage
    assert all(not p.is_symlink() for p in stage.rglob('*'))
    return manifest

def frame_evidence(output, owner, binding):
    # Source emits this record only after successful complete retained rendering.
    text = (output/'stdout.log').read_text(encoding='utf-8-sig')
    prefix = '[Vulkan] Render evidence: '
    rows = [json.loads(line[len(prefix):]) for line in text.splitlines() if line.startswith(prefix)]
    assert len(rows) == 1
    evidence = rows[0]
    for key,value in {'schema':'sirius-vulkan-render-v1','route':'retained','source_owner':'host','device_name':'AMD Radeon 780M Graphics','metric':'Kerr','precision':'fp32'}.items():
        assert type(evidence[key]) is str and evidence[key] == value
    for key,value in {'width':128,'height':128}.items():
        assert type(evidence[key]) is int and evidence[key] == value
    assert type(evidence['device_index']) is int and evidence['device_index'] >= 0
    assert type(evidence['work_items']) is int and evidence['work_items'] > 0
    assert type(evidence['accepted_intervals']) is int and evidence['accepted_intervals'] > 0
    stage_dispatches = evidence['retained_stage_dispatches']
    assert type(stage_dispatches) is list and len(stage_dispatches) == 7
    assert all(type(count) is int and count >= 0 for count in stage_dispatches) and sum(stage_dispatches) > 0
    scene_prefix = '[Session] Scene evidence: '
    scenes = [json.loads(line[len(scene_prefix):]) for line in text.splitlines() if line.startswith(scene_prefix)]
    assert len(scenes) == 1
    scene = scenes[0]
    assert scene['schema'] == 'sirius-render-scene-v1' and scene['backend'] == 'Vulkan' and scene['metric'] == 'Kerr'
    assert type(scene['width']) is int and scene['width'] == 128 and type(scene['height']) is int and scene['height'] == 128
    assert type(scene['samples_per_pixel']) is int and scene['samples_per_pixel'] == 1
    temporary = ROOT/'renders/current-original-frame-5b2df58b/native-frame'
    expected_temp = binding['windows_checkout_unc']+'\\'+str(temporary.relative_to(ROOT)).replace('/','\\')
    for key in ['TEMP','TMP']:
        assert owner['environment'][key].replace('/','\\').casefold() == expected_temp.casefold()
    assert not (temporary/'sirius_gpu_parse.ppm').exists()
    assert '[Session] Vulkan render complete: ' in text
    return {'renderer':evidence,'scene':scene,'original_temporary_ppm_removed':True,
            'cache_policy':'explicit prebound PowerShell MESA_SHADER_CACHE_DISABLE=true inherited by Start-Process; not an observed child-environment snapshot'}

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
    assert type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds'] == 1380
assert plan['owner_rss_bytes'] == 4294967296 and plan['terminal_observer_seconds'] == 60
bridge_seconds = plan['actions'][kind]['bridge_outer_seconds']
assert plan['minimum_wsl_available_kib'] == 2097152

verify()
stage = ROOT / plan['stage']
output = WORK / 'native' / kind
temporary = ROOT / 'renders/current-original-frame-5b2df58b/native-frame'
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
assert not (WORK / (kind+'-inputs-before.json')).exists() and not (WORK / (kind+'-inputs-after.json')).exists()
exe = stage / 'src/sirius/app/Release/sirius_render_tests.exe'
case = ACTION_NAMES[kind]
seconds = plan['actions'][kind]['owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
def quote(word): return "'"+word.replace("'","''")+"'"
owner_args = ['-Stage',win(exe.parent),'-Output',win(output),'-Filter',case,'-Seconds',str(seconds),
              '-RssBytes',str(plan['owner_rss_bytes']),'-TemporaryDirectory',win(temporary),'-BuildRoot',win(stage),
              '-ExpectedGateSha256',binding['gate_sha256'],'-LiveRevision',HEAD,'-ExecutableName',exe.name,'-InvocationKind','Control']
command = "$ErrorActionPreference='Stop';foreach($v in Get-ChildItem Env:){if($v.Name -match '^(LVP_|SLANG_|MESA_)' -and !($v.Name -eq 'MESA_SHADER_CACHE_DISABLE' -and $v.Value -eq 'true')){throw 'Unexpected inherited compiler/cache override'}};$env:MESA_SHADER_CACHE_DISABLE='true';& "+quote(win(WORK/'run_native_control.ps1'))+' '+ ' '.join(word if index % 2 == 0 else quote(word) for index,word in enumerate(owner_args))+';exit $LASTEXITCODE'
argv = [ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
assert plan['child_environment'] == {'MESA_SHADER_CACHE_DISABLE':'true','TEMP':str(temporary.relative_to(ROOT)),'TMP':str(temporary.relative_to(ROOT))}
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
    if kind in EVIDENCE_PROPERTIES:
        key, expected_value = EVIDENCE_PROPERTIES[kind]
        assert properties.get(key) == expected_value
        record['required_evidence_property'] = {key:properties[key]}
    assert not [p for p in cases[0].findall('./properties/property') if p.get('name') == 'fp64_evidence']
    record['original_frame_evidence'] = frame_evidence(output, owner, binding)
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
