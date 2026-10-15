"""Audit ownership/provider teardown for one completed numerical failure; never publish a scientific pass."""
from pathlib import Path, PurePosixPath
import json,os,signal,subprocess,sys,time
import hashlib,re,math,struct
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'): raise RuntimeError('Python assertions must be enabled')
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
if ROOT != WORK.parents[1]: raise RuntimeError('Run from the canonical checkout')
def document(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_failure_context(xml_path):
    xml = ET.parse(xml_path).getroot()
    suites = xml.findall('.//testsuite')
    cases = xml.findall('.//testcase')
    assert len(suites) == len(cases) == 1
    for node in [xml, *suites]:
        assert int(node.get('tests', '-1')) == 1 and int(node.get('failures', '-1')) == 1
        assert all(int(node.get(key, '0')) == 0 for key in ['errors', 'skipped', 'disabled'])
    case = cases[0]
    assert case.get('classname')+'.'+case.get('name') == ACTION_NAMES['kerr-image']
    assert case.get('status') == 'run' and case.get('result') == 'completed'
    failures = case.findall('failure')
    assert len(failures) == 1 and len(xml.findall('.//failure')) == 1
    assert not any(xml.findall('.//'+key) for key in ['error', 'skipped'])
    message = failures[0].get('message')
    assert isinstance(message, str) and 'output.has_value()' in message
    assert 'kernel: validate trace continuation: trace termination 5 cannot publish radiance' in message
    marker = '{"schema":"sirius-trace-continuation-failure-v1",'
    assert message.count(marker) == 1
    start = message.index(marker)
    context, end = json.JSONDecoder().raw_decode(message[start:])
    raw_context = message[start:start+end]
    assert set(context) == {'schema','scalar_bytes','packet_scalar_bytes','record_bytes','byte_order',
                            'action','row','active','active_rectangle','expected_identity',
                            'physical_submissions','maximum_attempts','packet_hex',
                            'previous_record_hex','current_record_hex','device_name_bytes_hex'}
    assert context['schema'] == 'sirius-trace-continuation-failure-v1'
    for key in ['scalar_bytes','packet_scalar_bytes','record_bytes','action','row','active',
                'physical_submissions','maximum_attempts']:
        assert type(context[key]) is int
    assert context['scalar_bytes'] == context['packet_scalar_bytes'] == 4
    assert context['record_bytes'] == 272 and context['byte_order'] == 'little'
    assert context['action'] == 1 and context['active'] == 64 and 0 <= context['row'] < 64
    assert context['physical_submissions'] > 0 and context['maximum_attempts'] == 2000
    def raw(key, size):
        text = context[key]
        assert isinstance(text, str) and re.fullmatch('[0-9a-f]{'+str(size*2)+'}', text)
        return bytes.fromhex(text)
    packet_bytes = raw('packet_hex', 70*4)
    previous = raw('previous_record_hex', 272)
    current = raw('current_record_hex', 272)
    packet = struct.unpack('<70f', packet_bytes)
    words = struct.unpack('<70I', packet_bytes)
    rectangle = context['active_rectangle']
    identity = context['expected_identity']
    assert type(rectangle) is list and len(rectangle) == 4 and all(type(x) is int for x in rectangle)
    assert type(identity) is list and len(identity) == 4 and all(type(x) is int for x in identity)
    x, y, width, height = rectangle
    assert x == 0 and 0 <= y < 64 and width == 64 and height == 1
    assert tuple(packet[31:35]) == tuple(rectangle) and packet[0] == packet[1] == 64.0
    assert packet[21] == 2000.0 and packet[46] == 0.0 and packet[68] == 1.0
    assert identity == [x+context['row'] % width, y+context['row'] // width, 0, words[69]] and words[69] > 0
    previous_control = struct.unpack_from('<4I', previous, 240)
    current_control = struct.unpack_from('<4I', current, 240)
    assert list(struct.unpack_from('<4I', previous, 256)) == identity
    assert list(struct.unpack_from('<4I', current, 256)) == identity
    assert current_control[0] == 2 and current_control[1] == 5
    assert current_control[3] <= current_control[2] <= context['maximum_attempts']
    assert all(math.isfinite(x) for x in struct.unpack_from('<60f', current))
    name_hex = context['device_name_bytes_hex']
    assert isinstance(name_hex, str) and re.fullmatch('(?:[0-9a-f]{2})+', name_hex)
    assert bytes.fromhex(name_hex).decode('utf-8') == 'AMD Radeon 780M Graphics'
    rollback = {
        'physical_vectors_0_through_11_equal': previous[:192] == current[:192],
        'accepted_affine_length_last_interval_beam_seed_equal': previous[196:208] == current[196:208],
        'source_and_volume_equal': previous[208:240] == current[208:240],
        'accepted_attempt_count_equal': previous_control[3] == current_control[3],
        'one_attempt_added': current_control[2] == previous_control[2]+1,
        'previous_running': previous_control[0] == 1 and previous_control[1] == 0,
    }
    return {'schema':'sirius-native-windows-kerr-failure-context-v1',
            'raw_context_json':raw_context, 'previous_control':list(previous_control),
            'current_control':list(current_control), 'rollback_comparisons':rollback,
            'rollback_comparison_is_observation_not_acceptance':True,
            'numerical_test_passed':False, 'radiance_publication_admitted':False}

def verify_manifest():
    stage = ROOT / plan['stage']
    manifest_path = stage / 'generated/sirius/windows-diagnostic-build-products.json'
    assert {'bytes': manifest_path.stat().st_size, 'sha256': sha(manifest_path)} == bindings['manifest_identity']
    assert bindings['files'][str(manifest_path.relative_to(ROOT))] == bindings['manifest_identity']
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    for key, value in {'schema':'sirius-windows-diagnostic-build-products-v1', 'source_revision':bindings['source_revision'],
                       'repository':'sunyaisblank/Sirius', 'run_id':'38711720394', 'run_attempt':'1'}.items():
        assert type(manifest[key]) is str and manifest[key] == value
    assert manifest['diagnostic_build_only'] is True and manifest['qualification_claimed'] is False
    assert bindings['artifact']['repository'] == manifest['repository'] and bindings['artifact']['run_id'] == manifest['run_id']
    assert bindings['artifact']['run_attempt'] == manifest['run_attempt'] and bindings['artifact']['name'] == 'sirius-windows-diagnostic-build-products'
    assert type(bindings['artifact']['id']) is int and bindings['artifact']['id'] > 0
    assert isinstance(manifest['files'], dict) and manifest['files']
    for key, row in manifest['files'].items():
        logical = PurePosixPath(key)
        assert isinstance(key, str) and key.startswith('bin/windows-msvc/') and '\\' not in key and ':' not in key
        assert not logical.is_absolute() and str(logical) == key and all(part not in ('', '.', '..') for part in key.split('/'))
        assert set(row) == {'bytes', 'sha256'} and type(row['bytes']) is int and row['bytes'] >= 0 and re.fullmatch('[0-9a-f]{64}', row['sha256'])
        relative = logical.relative_to('bin/windows-msvc')
        path = stage / relative
        assert path.is_file() and path.resolve(strict=True) == path and path.is_relative_to(stage)
        assert bindings['files'][str(path.relative_to(ROOT))] == row
    for relative in ['src/sirius/app/Release/sirius.exe', 'tests/backend/Release/sirius_backend_tests.exe',
                     'src/sirius/app/Release/vulkan-1.dll', 'tests/backend/Release/vulkan-1.dll']:
        assert 'bin/windows-msvc/' + relative in manifest['files']
    alignment_path = stage / 'src/sirius/app/Release/resources/model/alignment_receipt.json'
    alignment = json.loads(alignment_path.read_text(encoding='utf-8-sig'))
    assert alignment['source_revision'] == bindings['source_revision'] and alignment['source_tree_clean'] is True and alignment['alignment_mode'] == 'qualification'
    assert bindings['files'][str(alignment_path.relative_to(ROOT))] == manifest['files']['bin/windows-msvc/src/sirius/app/Release/resources/model/alignment_receipt.json']
    assert type(bindings['manifest_identity']['bytes']) is int and bindings['manifest_identity']['bytes'] > 0
    assert re.fullmatch('[0-9a-f]{64}', bindings['manifest_identity']['sha256'])
    gate_path = stage / 'generated/sirius/native_build_gate.json'
    assert sha(gate_path) == bindings['gate_sha256']
    assert bindings['files'][str(gate_path.relative_to(ROOT))] == {'bytes':gate_path.stat().st_size,'sha256':bindings['gate_sha256']}
    assert {'bytes':gate_path.stat().st_size,'sha256':bindings['gate_sha256']} == bindings['gate_identity']
    assert type(bindings['gate_identity']['bytes']) is int and bindings['gate_identity']['bytes'] > 0
    assert re.fullmatch('[0-9a-f]{64}', bindings['gate_sha256'])
    gate = json.loads(gate_path.read_text(encoding='utf-8-sig'))
    assert gate['source']['clean'] is True and gate['source']['revision'] == bindings['source_revision']
    for role, relative in [('sirius','src/sirius/app/Release/sirius.exe'), ('sirius_backend_tests','tests/backend/Release/sirius_backend_tests.exe')]:
        product = gate['tested_artifacts'][role]
        assert product['root'] == 'build' and product['path'] == relative
        assert {key:product[key] for key in ['bytes','sha256']} == manifest['files']['bin/windows-msvc/'+relative]
    return manifest

def read_case(mode,label,case,bindings):
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==bindings['runner_source_revision']
 assert not subprocess.check_output(['git','status','--porcelain'])
 bridge=document(WORK/(case+'-bridge.json'));owner=document(WORK/'native'/case/'owner.json')
 assert bridge['accepted_owner_execution'] and bridge['completed'] and bridge['returncode']==1
 assert bridge['binding_sha256']==sha(WORK/'execution-bindings.json')
 assert owner['status']=='completed' and owner['returncode']==1 and not owner['outer_stop'] and not owner['cleanup_errors']
 assert owner['owned_process_absent'] and owner['source_build_gate_unchanged'] and owner['test_executable_unchanged']
 assert owner['source_build_gate_sha256']==bindings['gate_sha256']
 assert owner['system_inventory_passed'] is False and owner['numerical_control_passed'] is False and owner['loaded_driver_modules_recorded'] is True
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
 failure_context=read_failure_context(WORK/'native'/case/'gtest.xml')
 assert failure_context==bridge['failure_context'] and bridge['numerical_test_passed'] is False and bridge['expected_completed_failure_captured'] is True
 receipts.append(WORK/'native'/case/'gtest.xml')
 if (WORK/'native'/case/'driver-modules.json').exists():
  module_path=WORK/'native'/case/'driver-modules.json';receipts.append(module_path)
  modules=document(module_path);assert isinstance(modules,list) and len(modules)==2
  loader=[row for row in modules if row['name'].lower()=='vulkan-1.dll']
  assert len(loader)==1 and any('amdvlk' in row['name'].lower() for row in modules)
  assert [row for row in modules if 'amdvlk' in row['name'].lower()]==[bindings['expected_amd_module']]
  relative='tests/backend/Release'
  expected=verify_manifest()['files']['bin/windows-msvc/'+relative+'/vulkan-1.dll']
  expected_path=bindings['windows_checkout_unc']+'\\'+(bindings['stage']+'/'+relative+'/vulkan-1.dll').replace('/','\\')
  assert loader[0]['path'].replace('/','\\').casefold()==expected_path.casefold()
  assert {key:loader[0][key] for key in ['bytes','sha256']}==expected
 whole_inputs={str(p.relative_to(WORK)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in receipts}
 return {'owner':owner,'bridge':bridge,'stdout_sha256':sha(WORK/'native'/case/'stdout.log'),'whole_inputs':whole_inputs,'failure_context':failure_context}
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
ACTION_NAMES = {'kerr-image':'KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow'}
ORDER = list(ACTION_NAMES)
OWNER_SECONDS = {'kerr-image':28800}
case,=sys.argv[1:];assert case in ACTION_NAMES;mode='native-failure-context';label=case
binding_path=WORK/'execution-bindings.json'
bindings=document(binding_path);binding_sha256=sha(binding_path)
plan=document(WORK/'plan.json');assert plan['admitted'] is True and bindings['admitted'] is True
assert plan['source_revision']==plan['runner_source_revision']==bindings['source_revision']==bindings['runner_source_revision']=='02ee51f85e45e8342923bae1d9f26d50f1d39e6d'
assert plan['stage']==bindings['stage']=='bin/windows-msvc/diagnostic-02ee51f'
assert plan['windows_checkout_unc']==bindings['windows_checkout_unc']==r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'
assert plan['order']==ORDER and set(plan['actions'])==set(ORDER)
for action_label in ORDER:
 action=plan['actions'][action_label]
 assert set(action)=={'kind','test','expected_exit','owner_seconds','bridge_outer_seconds'}
 assert action['kind']=='Control' and action['test']==ACTION_NAMES[action_label]
 assert type(action['expected_exit']) is int and action['expected_exit']==1
 assert type(action['owner_seconds']) is int and action['owner_seconds']==OWNER_SECONDS[action_label]
 assert type(action['bridge_outer_seconds']) is int and action['bridge_outer_seconds']==OWNER_SECONDS[action_label]+480
assert plan['owner_rss_bytes']==4294967296 and plan['minimum_wsl_available_kib']==2097152 and plan['terminal_observer_seconds']==60
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
record={'mode':mode,'label':label,'case':case,'read_case_pass':True,'numerical_test_passed':False,'expected_child_exit':1,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(WORK/'terminal-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Expected failed original Kerr capture and immutable prerequisite/source/input identities; bounded Windows birth/provider/accessible-consumer audit. No hidden/global/unconditional compiler-descendant absence or full qualification claim.'}
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
result={'pass':False,'audit_completed':True,'ownership_record_valid':True,'failure_context_complete':True,'numerical_test_passed':False,'case':case,'mode':mode,'label':label,'record':observed,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'loaded_provider_evidence_complete':record['loaded_provider_evidence_complete'],'source_revision':bindings['source_revision'],'native_build_gate_authority_used':True,'mandatory_or_full_runtime_admitted':False,'scope':'One expected completed original Kerr MinimumStep failure with exact host context and actual SDKloader+AMD. Ownership/provider/context audit accepted; numerical test remains failed. No image publication, independent-physics repair, precision/UI/fullMandatory/fullscience/domain/release claim.'}
assert not controller_signals,('controller signals before acceptance publication',controller_signals)
output=WORK/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
