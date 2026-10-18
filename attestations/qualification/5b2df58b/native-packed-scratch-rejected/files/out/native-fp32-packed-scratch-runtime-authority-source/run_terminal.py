"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path, PurePosixPath, PureWindowsPath
import json,os,signal,subprocess,sys,time
import hashlib,re
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'): raise RuntimeError('Python assertions must be enabled')
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
if ROOT != WORK.parents[1]: raise RuntimeError('Run from the canonical checkout')
def document(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify_manifest():
    stage = ROOT / plan['stage']
    assert stage.is_relative_to(ROOT/'bin/windows-msvc') and stage.resolve(strict=True) == stage
    manifest_path, producer_path = stage/'manifest.json', stage/'producer-result.json'
    def identity(path): return {'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    for path, name in [(manifest_path,'manifest_identity'),(producer_path,'producer_identity')]:
        assert identity(path) == bindings[name] == bindings['files'][str(path.relative_to(ROOT))]
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    producer = json.loads(producer_path.read_text(encoding='utf-8-sig'))
    assert manifest['schema'] == 'sirius-native-packed-ci-export-manifest-v1' and manifest['passed'] is True
    assert manifest['complete_host_or_sdk_or_tool_dso_backup'] is False
    assert manifest['source_revision'] == producer['source_revision'] == bindings['source_revision']
    assert producer['schema'] == 'sirius-native-packed-matched-ci-producer-v1'
    for key in ['passed','configure_passed','build_passed','source_unchanged','frozen_build_inputs_unchanged']:
        assert producer[key] is True, key
    for key in ['pe_or_provider_or_test_executed','mandatory_gate_executed','native_build_attestation_claimed','scientific_qualification_claimed','release_or_adoption_claimed']:
        assert producer[key] is False, key
    assert producer['failure'] is None and type(producer['private_compile_object_count']) is int and producer['private_compile_object_count'] == 5
    artifact = bindings['artifact']
    assert artifact['repository'] == producer['repository'] == 'sunyaisblank/Sirius'
    assert artifact['run_id'] == producer['run_id'] and artifact['run_attempt'] == producer['run_attempt']
    assert type(artifact['run_id']) is str and artifact['run_id'].isdigit()
    assert type(artifact['run_attempt']) is str and artifact['run_attempt'].isdigit()
    assert artifact['name'] == 'sirius-native-packed-matched-producer' and type(artifact['id']) is int and artifact['id'] > 0
    assert isinstance(manifest['files'],dict) and manifest['files'] and len({k.casefold() for k in manifest['files']}) == len(manifest['files'])
    for key,row in manifest['files'].items():
        logical=PurePosixPath(key)
        assert isinstance(key,str) and not logical.is_absolute() and str(logical)==key and '\\' not in key and ':' not in key
        assert all(part not in ('','.','..') for part in key.split('/')) and key != 'manifest.json'
        assert type(row['bytes']) is int and row['bytes'] >= 0 and re.fullmatch('[0-9a-f]{64}',row['sha256'])
        path=stage/logical
        assert path.is_file() and path.resolve(strict=True)==path and path.is_relative_to(stage)
        assert identity(path)=={k:row[k] for k in ['bytes','sha256']}==bindings['files'][str(path.relative_to(ROOT))]
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
    assert configure[6].startswith(option) and PureWindowsPath(configure[6][len(option):]) == PureWindowsPath(producer['source_root'])/'out/native-fp32-packed-scratch-windows-producer-source/matched-producer.cmake'
    assert build[:12]==[tools['cmake'],'--build','--preset','windows-msvc','--config','Release','--target','SiriusNativePackedMatchedProducer','--parallel','2','--','/verbosity:normal'] and len(build)==13
    assert build[12].startswith('/bl:') and PureWindowsPath(build[12][4:])==PureWindowsPath(producer['export_root'])/'build.binlog'
    assert all(PureWindowsPath(c['cwd'])==PureWindowsPath(producer['source_root']) for c in commands)
    cache=(stage/'metadata/config/CMakeCache.txt').read_text(encoding='utf-8-sig')
    for required in ['BUILD_TESTS:BOOL=ON','SIRIUS_MANDATORY_TESTS:BOOL=ON','SIRIUS_WERROR:BOOL=ON','SIRIUS_CONTRACT_MODE:STRING=2','SIRIUS_REQUIRE_VULKAN_RUNTIME:BOOL=ON','SIRIUS_ALIGNMENT_MODE:STRING=qualification','CMAKE_CONFIGURATION_TYPES:STRING=Release','BUILD_SHARED_LIBS:BOOL=OFF']:
        assert re.search('^'+re.escape(required)+r'\r?$',cache,re.M),required
    alignment=json.loads((stage/'metadata/config/generated/sirius/alignment_receipt.json').read_text(encoding='utf-8-sig'))
    assert alignment['source_revision']==bindings['source_revision'] and alignment['source_tree_clean'] is True and alignment['alignment_mode']=='qualification'
    assert set(bindings['selected_executables']) == {'baseline','candidate'}
    for arm in ['baseline','candidate']:
        selected='payload/'+arm+'/original_'+arm+'_controls.exe'
        assert bindings['selected_executables'][arm] == {'relative_path':selected,'identity':{k:manifest['files'][selected][k] for k in ['bytes','sha256']}}
        for name in ['original_'+arm+'_controls.exe','original_'+arm+'_controls.map','vulkan-1.dll','msvcp140.dll','vcruntime140.dll','vcruntime140_1.dll']:
            assert 'payload/'+arm+'/'+name in manifest['files'] and manifest['files']['payload/'+arm+'/'+name]['bytes'] > 0
    arm=ACTION_ARMS[case]
    assert plan['actions'][case]['arm']==arm
    executable='payload/'+arm+'/original_'+arm+'_controls.exe'
    assert bindings['selected_executables'][arm] == {'relative_path':executable,'identity':{k:manifest['files'][executable][k] for k in ['bytes','sha256']}}
    assert bindings['native_build_gate_authority_used'] is False and bindings['private_producer_authority_kind']=='private-ci-producer-result-and-manifest'
    return manifest

def verify_phase_prerequisites():
    required={'timestamp':[],'abba':['timestamp_acceptance'],'remaining':['timestamp_acceptance','benefit_acceptance']}[plan['control_profile']]
    assert set(plan['phase_prerequisites'])==set(required)
    for name in required:
        item=plan['phase_prerequisites'][name]
        assert set(item)=={'path','identity'}
        path=ROOT/item['path']
        assert path.resolve(strict=True)==path and path.is_relative_to(ROOT/'out')
        assert {'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}==item['identity']==bindings['files'][item['path']]
        proof=json.loads(path.read_text(encoding='utf-8-sig'))
        assert proof['passed'] is True and proof['source_revision']==bindings['source_revision']
        if name=='timestamp_acceptance': assert proof['control8_passed'] is True
        else: assert proof['benefit_passed'] is True and type(proof['checks_total']) is int and proof['checks_total']==22 and type(proof['checks_passed']) is int and proof['checks_passed']==22

def read_case(mode,label,case,bindings):
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==bindings['runner_source_revision']
 assert not subprocess.check_output(['git','status','--porcelain'])
 bridge=document(WORK/(case+'-bridge.json'));owner=document(WORK/'native'/case/'owner.json')
 assert bridge['accepted_owner_execution'] and bridge['completed'] and bridge['returncode']==0
 assert bridge['binding_sha256']==sha(WORK/'execution-bindings.json')
 assert owner['status']=='completed' and owner['returncode']==0 and not owner['outer_stop'] and not owner['cleanup_errors']
 assert owner['owned_process_absent'] and owner['source_build_gate_unchanged'] and owner['test_executable_unchanged']
 assert owner['source_build_gate_sha256']==bindings['producer_identity']['sha256'] and owner['source_private_manifest_sha256']==bindings['manifest_identity']['sha256']
 assert owner['source_private_manifest_unchanged'] is True and owner['private_producer_authority_kind']=='private-ci-producer-result-and-manifest' and owner['native_build_gate_authority_used'] is False
 assert owner['system_inventory_passed'] is False and owner['numerical_control_passed'] is True and owner['loaded_driver_modules_recorded'] is True
 verify_manifest()
 assert owner['source_revision']==bindings['source_revision'] and owner['live_source_revision']==bindings['runner_source_revision']
 for rel,item in bindings['files'].items():
  p=Path(rel) if rel.startswith('/') else ROOT/rel
  assert {'bytes':p.stat().st_size,'sha256':sha(p)}==item,rel
 assert bridge['whole_inputs_unchanged']
 assert (WORK/(case+'-inputs-before.json')).read_bytes()==(WORK/(case+'-inputs-after.json')).read_bytes()
 for rel,item in document(WORK/(case+'-inputs-before.json')).items():
  p=Path(rel) if rel.startswith('/') else ROOT/rel
  assert {'bytes':p.stat().st_size,'sha256':sha(p)}==item,rel
 for rel,item in bridge['predecessor_inputs'].items():
  path=ROOT/rel;assert {'bytes':path.stat().st_size,'sha256':sha(path)}==item
 receipts=[WORK/(case+'-bridge.json'),WORK/'native'/case/'owner.json',WORK/'native'/case/'stdout.log',WORK/'native'/case/'stderr.log',WORK/'native'/case/'memory.json',WORK/(case+'-inputs-before.json'),WORK/(case+'-inputs-after.json')]
 xml=ET.parse(WORK/'native'/case/'gtest.xml').getroot();cases=xml.findall('.//testcase')
 assert int(xml.get('tests','-1'))==1 and len(cases)==1
 assert all(int(node.get(key,'0'))==0 for node in [xml,*xml.findall('.//testsuite')] for key in ['failures','errors','skipped','disabled'])
 assert not any(xml.findall('.//'+key) for key in ['failure','error','skipped'])
 assert cases[0].get('classname')+'.'+cases[0].get('name')==ACTION_NAMES[case] and cases[0].get('status')=='run' and cases[0].get('result')=='completed'
 properties={p.get('name'):p.get('value') for p in cases[0].findall('./properties/property')}
 if case in EVIDENCE_PROPERTIES:
  key,expected_value=EVIDENCE_PROPERTIES[case]
  assert properties.get(key)==expected_value
 control_identity={'test':ACTION_NAMES[case],'arm':ACTION_ARMS[case]}
 receipts.append(WORK/'native'/case/'gtest.xml')
 if (WORK/'native'/case/'driver-modules.json').exists():
  module_path=WORK/'native'/case/'driver-modules.json';receipts.append(module_path)
  modules=document(module_path);assert isinstance(modules,list) and len(modules)==2
  loader=[row for row in modules if row['name'].lower()=='vulkan-1.dll']
  assert len(loader)==1 and any('amdvlk' in row['name'].lower() for row in modules)
  assert [row for row in modules if 'amdvlk' in row['name'].lower()]==[bindings['expected_amd_module']]
  relative='payload/'+ACTION_ARMS[case]
  expected_row=verify_manifest()['files'][relative+'/vulkan-1.dll']
  expected={key:expected_row[key] for key in ['bytes','sha256']}
  expected_path=bindings['windows_checkout_unc']+'\\'+(bindings['stage']+'/'+relative+'/vulkan-1.dll').replace('/','\\')
  assert loader[0]['path'].replace('/','\\').casefold()==expected_path.casefold()
  assert {key:loader[0][key] for key in ['bytes','sha256']}==expected
 whole_inputs={str(p.relative_to(WORK)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in receipts}
 return {'owner':owner,'bridge':bridge,'stdout_sha256':sha(WORK/'native'/case/'stdout.log'),'whole_inputs':whole_inputs,'original_control_identity':control_identity}
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
ACTION_NAMES = {'control-1': 'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission', 'control-2': 'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal', 'control-3': 'RetainedComputeTest.ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation', 'control-4': 'RetainedDopriTest.CurvedTransportConnectsSamplerAndProjection', 'control-5': 'RetainedDopriTest.ConnectedIntervalsPreserveTrialsAndIndependentBudgets', 'control-6': 'RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates', 'control-7': 'RetainedComputeTest.CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement', 'control-8': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'a1': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'b1': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'b2': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults', 'a2': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'}
ACTION_ARMS = {'control-1': 'candidate', 'control-2': 'candidate', 'control-3': 'candidate', 'control-4': 'candidate', 'control-5': 'candidate', 'control-6': 'candidate', 'control-7': 'candidate', 'control-8': 'candidate', 'a1': 'baseline', 'b1': 'candidate', 'b2': 'candidate', 'a2': 'baseline'}
PROFILES = {'timestamp': ['control-8'], 'abba': ['a1', 'b1', 'b2', 'a2'], 'remaining': ['control-1', 'control-2', 'control-3', 'control-4', 'control-5', 'control-6', 'control-7']}
OWNER_SECONDS = {name:300 for name in ACTION_NAMES}
EVIDENCE_PROPERTIES = {}
case,=sys.argv[1:];assert case in ACTION_NAMES;mode='native-packed-original-controls';label=case
binding_path=WORK/'execution-bindings.json'
bindings=document(binding_path);binding_sha256=sha(binding_path)
plan=document(WORK/'plan.json');assert plan['admitted'] is True and bindings['admitted'] is True
assert plan['schema']=='sirius-native-packed-runtime-plan-v1'
assert bindings['schema']=='sirius-native-packed-runtime-bindings-v1'
assert plan['control_profile']==bindings['control_profile'] and plan['control_profile'] in PROFILES
ORDER=PROFILES[plan['control_profile']]
assert case in ORDER and set(bindings['actions'])==set(ORDER) and set(bindings['terminal_actions'])==set(ORDER)
assert plan['source_revision']==bindings['source_revision'] and plan['runner_source_revision']==bindings['runner_source_revision']
assert all(re.fullmatch('[0-9a-f]{40}',plan[k]) for k in ['source_revision','runner_source_revision'])
assert plan['stage']==bindings['stage'] and plan['stage'].startswith('bin/windows-msvc/')
assert plan['windows_checkout_unc']==bindings['windows_checkout_unc']==r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'
assert plan['order']==ORDER and set(plan['actions'])==set(ORDER)
for action_label in ORDER:
 action=plan['actions'][action_label]
 assert set(action)=={'kind','test','arm','expected_exit','owner_seconds','bridge_outer_seconds'}
 assert action['kind']=='Control' and action['test']==ACTION_NAMES[action_label] and action['arm']==ACTION_ARMS[action_label]
 assert type(action['expected_exit']) is int and action['expected_exit']==0
 assert type(action['owner_seconds']) is int and action['owner_seconds']==OWNER_SECONDS[action_label]
 assert type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds']==1380
assert plan['owner_rss_bytes']==4294967296 and plan['minimum_wsl_available_kib']==2097152 and plan['terminal_observer_seconds']==60
verify_phase_prerequisites()
observed=read_case(mode,label,case,bindings)
prefix='terminal-'+case
receipt=WORK/(prefix+'.json');bridge=WORK/(prefix+'-bridge.json')
temp=WORK/'temp'/prefix;assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir(parents=True)
def win(path):
 path=path.resolve();assert path.is_relative_to(ROOT)
 return bindings['windows_checkout_unc']+'\\'+str(path.relative_to(ROOT)).replace('/','\\')
def quote(word):return "'"+word.replace("'","''")+"'"
audit=win(WORK/'terminal-audit.ps1');workspace=win(WORK);bootstrap=win(temp/'bootstrap.json');run=case
command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(bootstrap)+" -Encoding UTF8;& "+quote(audit)+" -Workspace "+quote(workspace)+" -RunName "+quote('native/'+run)+" -Receipt "+quote(receipt.name)+' -BuildRoot '+quote(win(ROOT/bindings['stage']))
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(WORK/(prefix+'-emergency-cleanup.json'))]
assert arguments==bindings['terminal_actions'][case]['argv'] and cleanup_args==bindings['terminal_actions'][case]['cleanup_argv']
record={'mode':mode,'label':label,'case':case,'read_case_pass':True,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(WORK/'terminal-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Selected original native-AMD packed control and private compile-only source/product/input identities; bounded Windows birth/provider/accessible-consumer audit. No hidden/global/unconditional compiler-descendant absence or full qualification claim.'}
started=time.monotonic()
with (WORK/(prefix+'.stdout')).open('wb') as out,(WORK/(prefix+'.stderr')).open('wb') as err:
 assert not controller_signals,('interrupted before spawn',controller_signals)
 child=subprocess.Popen(arguments,stdout=out,stderr=err,start_new_session=True)
 try:
  record['linux_child_pid']=child.pid;record['linux_child_birth']=Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
  while child.poll() is None:
   if controller_signals:
    record['outer_stop']=True;record['errors'].append('controller_signal');break
   remaining=60-(time.monotonic()-started)
   if remaining<=0:raise subprocess.TimeoutExpired(arguments,60)
   try:child.wait(timeout=min(0.25,remaining))
   except subprocess.TimeoutExpired:pass
  code=child.returncode
 except BaseException as error:
  record['outer_stop']=True;record['errors'].append(repr(error))
 finally:
  if controller_signals:
   record['outer_stop']=True
   if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
  if record['outer_stop'] or child.poll() is None:
   try:
    with (WORK/(prefix+'-emergency.stdout')).open('wb') as cout,(WORK/(prefix+'-emergency.stderr')).open('wb') as cerr:
     c=subprocess.run(cleanup_args,stdout=cout,stderr=cerr,timeout=60)
    assert c.returncode==0
    cleanup=document(WORK/(prefix+'-emergency-cleanup.json'));assert cleanup['recorded_births_absent'] and not cleanup['errors']
   except BaseException as error:record['errors'].append(repr(error))
  try:code=child.wait(timeout=15)
  except BaseException as error:
   record['errors'].append(repr(error))
   try:os.killpg(child.pid,signal.SIGKILL)
   except ProcessLookupError:pass
   except BaseException as error:record['errors'].append(repr(error))
   try:code=child.wait(timeout=10)
   except BaseException as error:record['errors'].append(repr(error));code=child.returncode
if controller_signals:
 record['outer_stop']=True
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
record['returncode']=code;record['elapsed_seconds']=time.monotonic()-started
try:current=Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
except (FileNotFoundError,ProcessLookupError):current=None
record['linux_child_birth_absent']=record.get('linux_child_birth') is not None and current!=record['linux_child_birth']
record['accepted']=code==0 and not record['outer_stop'] and not record['errors'] and record['linux_child_birth_absent']
if record['accepted']:
 try:
  terminal=document(receipt);assert terminal['owned_births_absent'] and terminal['all_observed_module_seals_match']
  assert terminal['actual_loaded_modules_observed'] is True
  record['loaded_provider_evidence_complete']=terminal['actual_loaded_modules_observed'] is True
  assert not terminal['consumer_scan']['matches'] and not terminal['consumer_scan']['accessible_compilers']
  record['terminal_receipt_sha256']=sha(receipt)
  record['bootstrap']=document(temp/'bootstrap.json')
  assert sha(binding_path)==binding_sha256
  assert read_case(mode,label,case,bindings)==observed
  assert not controller_signals
  record['source_binding_owner_bridge_action_post_sealed']=True
 except BaseException as error:
  record['accepted']=False;record['errors'].append('post-audit: '+repr(error))
if controller_signals:
 record['outer_stop']=True;record['accepted']=False
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
bridge.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('mode','label','case','returncode','elapsed_seconds','accepted')}),flush=True)
assert record['accepted']
result={'pass':True,'case':case,'mode':mode,'label':label,'record':observed,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'loaded_provider_evidence_complete':record['loaded_provider_evidence_complete'],'source_revision':bindings['source_revision'],'native_build_gate_authority_used':False,'private_producer_authority_used':True,'private_producer_authority_kind':'private-ci-producer-result-and-manifest','mandatory_or_full_runtime_admitted':False,'scope':'One original native-AMD control process only: exact case,completed zero-FESD XML,private producer+manifest authority,staged SDK+AMD and held-birth teardown. Full numerical/routes/timestamps/readbacks/performance remain separate readers; no frame/fullscience/UI/Mandatory/release qualification.'}
assert not controller_signals,('controller signals before acceptance publication',controller_signals)
output=WORK/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
