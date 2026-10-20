"""Own one original native packed-control action under private compile-only authority."""
from pathlib import Path, PurePosixPath, PureWindowsPath
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
ACTION_NAMES = {'control-1': 'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission', 'control-2': 'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal', 'control-3': 'RetainedComputeTest.ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation', 'control-4': 'RetainedDopriTest.CurvedTransportConnectsSamplerAndProjection', 'control-5': 'RetainedDopriTest.ConnectedIntervalsPreserveTrialsAndIndependentBudgets', 'control-6': 'RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates', 'control-7': 'RetainedComputeTest.CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement', 'control-8': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'a1': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'b1': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'b2': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'a2': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'}
ACTION_ARMS = {'control-1': 'candidate', 'control-2': 'candidate', 'control-3': 'candidate', 'control-4': 'candidate', 'control-5': 'candidate', 'control-6': 'candidate', 'control-7': 'candidate', 'control-8': 'candidate', 'a1': 'baseline', 'b1': 'candidate', 'b2': 'candidate', 'a2': 'baseline'}
PROFILES = {'timestamp': ['control-8'], 'abba': ['a1', 'b1', 'b2', 'a2'], 'remaining': ['control-1', 'control-2', 'control-3', 'control-4', 'control-5', 'control-6', 'control-7']}
OWNER_SECONDS = {name:300 for name in ACTION_NAMES}
EVIDENCE_PROPERTIES = {}
kind, = sys.argv[1:]
assert kind in ACTION_NAMES
plan = json.loads((WORK / 'plan.json').read_text())
binding_path = WORK / 'execution-bindings.json'
binding = json.loads(binding_path.read_text())
HEAD = plan['runner_source_revision']
SOURCE = plan['source_revision']
assert plan['schema'] == 'sirius-native-packed-runtime-plan-v1'
assert binding['schema'] == 'sirius-native-packed-runtime-bindings-v1'
assert plan['control_profile'] == binding['control_profile'] and plan['control_profile'] in PROFILES
ORDER = PROFILES[plan['control_profile']]
assert kind in ORDER and set(binding['actions']) == set(ORDER) and set(binding['terminal_actions']) == set(ORDER)
assert plan['admitted'] is True and binding['admitted'] is True
assert re.fullmatch('[0-9a-f]{40}',SOURCE) and re.fullmatch('[0-9a-f]{40}',HEAD)
assert plan['stage'] == binding['stage'] and plan['stage'].startswith('bin/windows-msvc/')
assert plan['canonical_CXX_source_revision'] == binding['canonical_CXX_source_revision'] == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert plan['runner_source_revision'] == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert plan['stage'] == 'bin/windows-msvc/native-coefficient-sharing-matched-current'
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
    assert stage.is_relative_to(ROOT/'bin/windows-msvc') and stage.resolve(strict=True) == stage
    manifest_path, producer_path = stage/'manifest.json', stage/'producer-result.json'
    def identity(path): return {'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    for path, name in [(manifest_path,'manifest_identity'),(producer_path,'producer_identity')]:
        assert identity(path) == binding[name] == binding['files'][str(path.relative_to(ROOT))]
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    producer = json.loads(producer_path.read_text(encoding='utf-8-sig'))
    assert manifest['schema'] == 'sirius-native-packed-ci-export-manifest-v1' and manifest['passed'] is True
    assert manifest['complete_host_or_sdk_or_tool_dso_backup'] is False
    assert manifest['source_revision'] == producer['source_revision'] == SOURCE
    assert producer['schema'] == 'sirius-native-packed-matched-ci-producer-v1'
    assert producer['producer_profile'] == 'native-coefficient-sharing-one-shader-only'
    assert producer['canonical_CXX_source_revision'] == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
    assert type(producer['expected_private_compile_objects']) is int and producer['expected_private_compile_objects'] == 5
    for key in ['passed','configure_passed','build_passed','source_unchanged','frozen_build_inputs_unchanged']:
        assert producer[key] is True, key
    for key in ['pe_or_provider_or_test_executed','mandatory_gate_executed','native_build_attestation_claimed','scientific_qualification_claimed','release_or_adoption_claimed']:
        assert producer[key] is False, key
    assert producer['failure'] is None and type(producer['private_compile_object_count']) is int and producer['private_compile_object_count'] == 5
    artifact = binding['artifact']
    assert artifact['repository'] == producer['repository'] == 'sunyaisblank/Sirius'
    assert artifact['run_id'] == producer['run_id'] and artifact['run_attempt'] == producer['run_attempt']
    assert type(artifact['run_id']) is str and artifact['run_id'].isdigit()
    assert type(artifact['run_attempt']) is str and artifact['run_attempt'].isdigit()
    assert artifact['name'] == 'sirius-native-coefficient-sharing-matched-producer' and type(artifact['id']) is int and artifact['id'] > 0
    assert isinstance(manifest['files'],dict) and manifest['files'] and len({k.casefold() for k in manifest['files']}) == len(manifest['files'])
    for key,row in manifest['files'].items():
        logical=PurePosixPath(key)
        assert isinstance(key,str) and not logical.is_absolute() and str(logical)==key and '\\' not in key and ':' not in key
        assert all(part not in ('','.','..') for part in key.split('/')) and key != 'manifest.json'
        assert type(row['bytes']) is int and row['bytes'] >= 0 and re.fullmatch('[0-9a-f]{64}',row['sha256'])
        path=stage/logical
        assert path.is_file() and path.resolve(strict=True)==path and path.is_relative_to(stage)
        assert identity(path)=={k:row[k] for k in ['bytes','sha256']}==binding['files'][str(path.relative_to(ROOT))]
        assert PureWindowsPath(row['path']) == PureWindowsPath(producer['export_root'])/PureWindowsPath(key)
    assert {str(p.relative_to(stage)).replace('\\','/') for p in stage.rglob('*') if p.is_file()} == set(manifest['files'])|{'manifest.json'}
    assert all(not p.is_symlink() for p in stage.rglob('*'))
    assert not any(p.name in ['native_build_gate.json','mandatory_gate.json'] for p in stage.rglob('*'))
    source_before=json.loads((stage/'source-before-configure.json').read_text(encoding='utf-8-sig'))
    before=json.loads((stage/'build-inputs-before.json').read_text(encoding='utf-8-sig'))
    after=json.loads((stage/'build-inputs-after.json').read_text(encoding='utf-8-sig'))
    assert isinstance(before,dict) and before and before==after
    assert isinstance(source_before,dict) and source_before and all(before.get(k)==v for k,v in source_before.items())
    commands=producer['commands']
    assert isinstance(commands,list) and [c['label'] for c in commands]==['configure','build']
    assert json.loads((stage/'commands.json').read_text(encoding='utf-8-sig')) == commands
    assert all(type(c['exit_code']) is int and c['exit_code']==0 for c in commands)
    tools=json.loads((stage/'selected-tools.json').read_text(encoding='utf-8-sig'))
    configure=commands[0]['argv'];build=commands[1]['argv']
    assert configure[:6]==[tools['cmake'],'--preset','windows-msvc','-DSIRIUS_ALIGNMENT_MODE=qualification','-DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON','-DCMAKE_CONFIGURATION_TYPES:STRING=Release'] and len(configure)==7
    option='-DCMAKE_PROJECT_Sirius_INCLUDE='
    assert configure[6].startswith(option) and PureWindowsPath(configure[6][len(option):]) == PureWindowsPath(producer['source_root'])/'out/native-coefficient-sharing-windows-source/matched-producer.cmake'
    assert build[:12]==[tools['cmake'],'--build','--preset','windows-msvc','--config','Release','--target','SiriusNativeCoefficientSharingMatchedProducer','--parallel','2','--','/verbosity:normal'] and len(build)==13
    assert build[12].startswith('/bl:') and PureWindowsPath(build[12][4:])==PureWindowsPath(producer['export_root'])/'build.binlog'
    assert all(PureWindowsPath(c['cwd'])==PureWindowsPath(producer['source_root']) for c in commands)
    cache=(stage/'metadata/config/CMakeCache.txt').read_text(encoding='utf-8-sig')
    for required in ['BUILD_TESTS:BOOL=ON','SIRIUS_MANDATORY_TESTS:BOOL=ON','SIRIUS_WERROR:BOOL=ON','SIRIUS_CONTRACT_MODE:STRING=2','SIRIUS_REQUIRE_VULKAN_RUNTIME:BOOL=ON','SIRIUS_ALIGNMENT_MODE:STRING=qualification','CMAKE_CONFIGURATION_TYPES:STRING=Release','BUILD_SHARED_LIBS:BOOL=OFF']:
        assert re.search('^'+re.escape(required)+r'\r?$',cache,re.M),required
    alignment=json.loads((stage/'metadata/config/generated/sirius/alignment_receipt.json').read_text(encoding='utf-8-sig'))
    assert alignment['source_revision']==SOURCE and alignment['source_tree_clean'] is True and alignment['alignment_mode']=='qualification'
    assert set(binding['selected_executables']) == {'baseline','candidate'}
    for arm in ['baseline','candidate']:
        selected='payload/'+arm+'/original_'+arm+'_controls.exe'
        assert binding['selected_executables'][arm] == {'relative_path':selected,'identity':{k:manifest['files'][selected][k] for k in ['bytes','sha256']}}
        for name in ['original_'+arm+'_controls.exe','original_'+arm+'_controls.map','vulkan-1.dll','msvcp140.dll','vcruntime140.dll','vcruntime140_1.dll']:
            assert 'payload/'+arm+'/'+name in manifest['files'] and manifest['files']['payload/'+arm+'/'+name]['bytes'] > 0
    arm=ACTION_ARMS[kind]
    assert plan['actions'][kind]['arm']==arm
    executable='payload/'+arm+'/original_'+arm+'_controls.exe'
    assert binding['selected_executables'][arm] == {'relative_path':executable,'identity':{k:manifest['files'][executable][k] for k in ['bytes','sha256']}}
    assert binding['native_build_gate_authority_used'] is False and binding['private_producer_authority_kind']=='private-ci-producer-result-and-manifest'
    return manifest

def verify_phase_prerequisites():
    required={'timestamp':[],'abba':['timestamp_acceptance'],'remaining':['timestamp_acceptance','benefit_acceptance']}[plan['control_profile']]
    assert set(plan['phase_prerequisites'])==set(required)
    for name in required:
        item=plan['phase_prerequisites'][name]
        assert set(item)=={'path','identity'}
        path=ROOT/item['path']
        assert path.resolve(strict=True)==path and path.is_relative_to(ROOT/'out')
        assert {'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}==item['identity']==binding['files'][item['path']]
        proof=json.loads(path.read_text(encoding='utf-8-sig'))
        assert proof['passed'] is True and proof['source_revision']==binding['source_revision']
        if name=='timestamp_acceptance': assert proof['control8_passed'] is True
        else: assert proof['benefit_passed'] is True and type(proof['checks_total']) is int and proof['checks_total']==22 and type(proof['checks_passed']) is int and proof['checks_passed']==22

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
    assert set(action) == {'kind','test','arm','expected_exit','owner_seconds','bridge_outer_seconds'}
    assert action['kind'] == 'Control' and action['test'] == ACTION_NAMES[label] and action['arm'] == ACTION_ARMS[label]
    assert type(action['expected_exit']) is int and action['expected_exit'] == 0
    assert type(action['owner_seconds']) is int and action['owner_seconds'] == OWNER_SECONDS[label]
    assert type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds'] == 1380
assert plan['owner_rss_bytes'] == 4294967296 and plan['terminal_observer_seconds'] == 60
bridge_seconds = plan['actions'][kind]['bridge_outer_seconds']
assert plan['minimum_wsl_available_kib'] == 2097152

verify_phase_prerequisites()
verify()
stage = ROOT / plan['stage']
output = WORK / 'native' / kind
temporary = ROOT / 'renders' / WORK.name / kind
record_path = WORK / (kind + '-bridge.json')
assert not output.exists() and not temporary.exists() and not record_path.exists()
assert not (WORK / (kind+'-inputs-before.json')).exists() and not (WORK / (kind+'-inputs-after.json')).exists()
arm = ACTION_ARMS[kind]
exe = stage / 'payload' / arm / ('original_'+arm+'_controls.exe')
case = ACTION_NAMES[kind]
seconds = plan['actions'][kind]['owner_seconds']
ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
def quote(word): return "'"+word.replace("'","''")+"'"
owner_args = ['-Stage',win(exe.parent),'-Output',win(output),'-Filter',case,'-Seconds',str(seconds),
              '-RssBytes',str(plan['owner_rss_bytes']),'-TemporaryDirectory',win(temporary),'-BuildRoot',win(stage),
              '-ExpectedManifestSha256',binding['manifest_identity']['sha256'],'-ExpectedProducerSha256',binding['producer_identity']['sha256'],'-Arm',arm,'-LiveRevision',HEAD,'-ExecutableName',exe.name,'-InvocationKind','Control']
command = "$ErrorActionPreference='Stop';foreach($v in Get-ChildItem Env:){if($v.Name -match '^(LVP_|SLANG_|MESA_)'){throw 'Unexpected inherited compiler/cache override'}};& "+quote(win(WORK/'run_native_control.ps1'))+' '+ ' '.join(word if index % 2 == 0 else quote(word) for index,word in enumerate(owner_args))+';exit $LASTEXITCODE'
argv = [ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
assert plan['child_environment'][kind] == {'TEMP':str(temporary.relative_to(ROOT)),'TMP':str(temporary.relative_to(ROOT))}
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
    assert owner['source_revision'] == SOURCE and owner['live_source_revision'] == HEAD and owner['source_build_gate_sha256'] == binding['producer_identity']['sha256']
    assert owner['source_private_manifest_sha256'] == binding['manifest_identity']['sha256'] and owner['source_private_manifest_unchanged'] is True
    assert owner['private_producer_authority_kind'] == 'private-ci-producer-result-and-manifest' and owner['native_build_gate_authority_used'] is False
    record['native_build_gate_authority_used'] = False
    record['private_producer_authority_used'] = True
    record['private_producer_authority_kind'] = 'private-ci-producer-result-and-manifest'
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
    record['original_control_identity'] = {'test':case,'arm':ACTION_ARMS[kind]}
    record['loaded_provider_evidence_complete'] = owner['loaded_driver_modules_recorded'] is True
    if owner['loaded_driver_modules_recorded']:
        modules=json.loads((output/'driver-modules.json').read_text(encoding='utf-8-sig'))
        assert isinstance(modules,list) and len(modules) == 2
        loader=[row for row in modules if row['name'].lower() == 'vulkan-1.dll']
        assert len(loader) == 1 and any('amdvlk' in row['name'].lower() for row in modules)
        assert [row for row in modules if 'amdvlk' in row['name'].lower()] == [binding['expected_amd_module']]
        expected_row=verify_manifest()['files']['payload/'+arm+'/vulkan-1.dll']
        expected={key:expected_row[key] for key in ['bytes','sha256']}
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
