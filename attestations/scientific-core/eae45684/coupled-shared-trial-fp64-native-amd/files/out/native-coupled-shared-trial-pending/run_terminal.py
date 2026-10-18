"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path, PurePosixPath
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
    manifest_path = stage / 'generated/sirius/windows-diagnostic-build-products.json'
    assert {'bytes': manifest_path.stat().st_size, 'sha256': sha(manifest_path)} == bindings['manifest_identity']
    assert bindings['files'][str(manifest_path.relative_to(ROOT))] == bindings['manifest_identity']
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    for key, value in {'schema':'sirius-windows-diagnostic-build-products-v1', 'source_revision':bindings['source_revision'],
                       'repository':'sunyaisblank/Sirius', 'run_id':'39082394780', 'run_attempt':'1'}.items():
        assert type(manifest[key]) is str and manifest[key] == value
    assert manifest['diagnostic_build_only'] is True and manifest['qualification_claimed'] is False
    assert bindings['artifact']['repository'] == manifest['repository'] and bindings['artifact']['run_id'] == manifest['run_id']
    assert bindings['artifact']['run_attempt'] == manifest['run_attempt'] and bindings['artifact']['name'] == 'sirius-windows-diagnostic-build-products'
    assert type(bindings['artifact']['id']) is int and bindings['artifact']['id'] > 0
    assert isinstance(manifest['files'], dict) and manifest['files']
    overlay_path = stage / 'generated/sirius/private-coupled-shader-overlay.json'
    overlay_identity = bindings['overlay_identity']
    assert set(overlay_identity) == {'bytes','sha256'} and type(overlay_identity['bytes']) is int and overlay_identity['bytes'] > 0
    assert re.fullmatch('[0-9a-f]{64}', overlay_identity['sha256'])
    assert overlay_path.is_file() and overlay_path.resolve(strict=True) == overlay_path
    assert {'bytes':overlay_path.stat().st_size,'sha256':sha(overlay_path)} == overlay_identity
    assert bindings['files'][str(overlay_path.relative_to(ROOT))] == overlay_identity
    overlay = json.loads(overlay_path.read_text(encoding='utf-8-sig'))
    assert set(overlay) == {'schema','private_shader_experiment','source_revision','stage','base_manifest_identity','base_native_gate_identity','override','provenance'}
    assert overlay['schema'] == 'sirius-private-coupled-shader-overlay-v1' and overlay['private_shader_experiment'] is True
    assert overlay['source_revision'] == bindings['source_revision'] and overlay['stage'] == plan['stage']
    assert overlay['base_manifest_identity'] == bindings['manifest_identity'] and overlay['base_native_gate_identity'] == bindings['native_gate_identity']
    override_key = 'bin/windows-msvc/tests/backend/Release/resources/kernels/coupled_probe_fp64.spv'
    old_row = {'bytes': 1225888, 'sha256': 'ee60ee290ca0f18c3031ca456a731e88568e6e268b5b032bc0cb840df999a667'}
    new_row = {'bytes': 878700, 'sha256': 'fbdd17e7fb3cd1e541bc7e51ce9ca407462f54482e570cb2adf2c2dd8fd71767'}
    assert overlay['override'] == {'key':override_key,'old':old_row,'new':new_row}
    assert manifest['files'][override_key] == old_row
    provenance_paths = {'candidate_entry': 'out/coupled-shared-trial-body/coupled_probe.slang', 'candidate_module': 'out/coupled-shared-trial-body/gr_coupled.slang', 'source_contract': 'out/coupled-shared-trial-body/source-contracts.json', 'source_review': 'out/coupled-shared-trial-body/independent-source-review.json', 'compiler_product': 'out/coupled-shared-trial-screen/candidate-fp64/coupled_probe_fp64.spv', 'compiler_selection': 'out/coupled-shared-trial-screen/root-owner-and-source-selection.json', 'compiler_posttools': 'out/coupled-shared-trial-screen/root-posttools.json', 'compiler_static': 'out/coupled-shared-trial-screen/root-static-interface-and-footprint.json', 'independent_emitted_review': 'out/coupled-shared-trial-screen/independent-emitted-review.json'}
    assert overlay['provenance'] == bindings['shader_provenance'] and set(overlay['provenance']) == set(provenance_paths)
    for label, relative in provenance_paths.items():
        item = overlay['provenance'][label]
        assert set(item) == {'path','identity'} and item['path'] == relative
        row = item['identity']
        assert set(row) == {'bytes','sha256'} and type(row['bytes']) is int and row['bytes'] > 0 and re.fullmatch('[0-9a-f]{64}',row['sha256'])
        path = ROOT / relative
        assert path.is_file() and path.resolve(strict=True) == path
        assert bindings['files'][relative] == row and {'bytes':path.stat().st_size,'sha256':sha(path)} == row
    assert overlay['provenance']['compiler_product']['identity'] == new_row
    assert overlay['provenance']['candidate_entry']['identity'] == {'bytes':6259,'sha256':'748058bbfecfab730d2e652d21c5c10eadd9414c14410829d3180e5965774ddc'}
    assert overlay['provenance']['candidate_module']['identity'] == {'bytes':13076,'sha256':'8542f1092485d59c6e8f2ca9640a5ce8732935bfa5c7661def3bd99efff053ce'}
    for key, row in manifest['files'].items():
        logical = PurePosixPath(key)
        assert isinstance(key, str) and key.startswith('bin/windows-msvc/') and '\\' not in key and ':' not in key
        assert not logical.is_absolute() and str(logical) == key and all(part not in ('', '.', '..') for part in key.split('/'))
        assert set(row) == {'bytes', 'sha256'} and type(row['bytes']) is int and row['bytes'] >= 0 and re.fullmatch('[0-9a-f]{64}', row['sha256'])
        relative = logical.relative_to('bin/windows-msvc')
        path = stage / relative
        assert path.is_file() and path.resolve(strict=True) == path and path.is_relative_to(stage)
        expected = new_row if key == override_key else row
        assert bindings['files'][str(path.relative_to(ROOT))] == expected
        assert {'bytes':path.stat().st_size,'sha256':sha(path)} == expected
    for relative in ['src/sirius/app/Release/sirius.exe', 'tests/backend/Release/sirius_backend_tests.exe',
                     'src/sirius/app/Release/vulkan-1.dll', 'tests/backend/Release/vulkan-1.dll']:
        assert 'bin/windows-msvc/' + relative in manifest['files']
    alignment_path = stage / 'src/sirius/app/Release/resources/model/alignment_receipt.json'
    alignment = json.loads(alignment_path.read_text(encoding='utf-8-sig'))
    assert alignment['source_revision'] == bindings['source_revision'] and alignment['source_tree_clean'] is True and alignment['alignment_mode'] == 'qualification'
    assert bindings['files'][str(alignment_path.relative_to(ROOT))] == manifest['files']['bin/windows-msvc/src/sirius/app/Release/resources/model/alignment_receipt.json']
    gate_path = stage / 'generated/sirius/native_build_gate.json'
    assert sha(gate_path) == bindings['gate_sha256']
    assert bindings['files'][str(gate_path.relative_to(ROOT))] == {'bytes':gate_path.stat().st_size,'sha256':bindings['gate_sha256']}
    gate_identity = bindings['native_gate_identity']
    assert set(gate_identity) == {'bytes','sha256'}
    assert type(gate_identity['bytes']) is int and gate_identity['bytes'] > 0
    assert re.fullmatch('[0-9a-f]{64}', gate_identity['sha256'])
    assert {'bytes':gate_path.stat().st_size,'sha256':bindings['gate_sha256']} == gate_identity
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
 assert bridge['accepted_owner_execution'] and bridge['completed'] and bridge['returncode']==0
 assert bridge['binding_sha256']==sha(WORK/'execution-bindings.json')
 assert owner['status']=='completed' and owner['returncode']==0 and not owner['outer_stop'] and not owner['cleanup_errors']
 assert owner['owned_process_absent'] and owner['source_build_gate_unchanged'] and owner['test_executable_unchanged']
 assert owner['source_build_gate_sha256']==bindings['gate_sha256']
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
 assert [p.get('value') for p in cases[0].findall('./properties/property') if p.get('name')=='fp64_evidence']==['actual_four_column_joint_rk4_and_independent_scalar_neighbours']
 assert [p.get('value') for p in cases[0].findall('./properties/property') if p.get('name')=='coupled_supports_fp64']==['true']
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
 return {'owner':owner,'bridge':bridge,'stdout_sha256':sha(WORK/'native'/case/'stdout.log'),'whole_inputs':whole_inputs}
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
ACTION_NAMES = {
    'coupled-fp64':'KernelParity.FourColumnRK4MapMatchesIndependentNeighboursFp64'}
PROFILES = {'private-coupled-shared-trial-fp64':['coupled-fp64']}
OWNER_SECONDS = {name:900 for name in ACTION_NAMES}
EVIDENCE_PROPERTIES = {'coupled-fp64':('fp64_evidence','actual_four_column_joint_rk4_and_independent_scalar_neighbours')}
case,=sys.argv[1:];assert case in ACTION_NAMES;mode='native-private-coupled-shared-trial';label=case
binding_path=WORK/'execution-bindings.json'
bindings=document(binding_path);binding_sha256=sha(binding_path)
plan=document(WORK/'plan.json');assert plan['admitted'] is True and bindings['admitted'] is True
assert plan['schema']=='sirius-native-windows-private-coupled-shared-trial-plan-v1'
assert bindings['schema']=='sirius-native-windows-private-coupled-shared-trial-bindings-v1'
assert plan['control_profile']==bindings['control_profile'] and plan['control_profile'] in PROFILES
ORDER=PROFILES[plan['control_profile']]
assert case in ORDER and set(bindings['actions'])==set(ORDER) and set(bindings['terminal_actions'])==set(ORDER)
assert plan['source_revision']==plan['runner_source_revision']==bindings['source_revision']==bindings['runner_source_revision']=='eae45684ee6f2d1b1baf9673a53ec9c349f936e1'
assert plan['stage']==bindings['stage']=='bin/windows-msvc/coupled-shared-trial-eae45684'
assert plan['windows_checkout_unc']==bindings['windows_checkout_unc']==r'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius'
assert plan['order']==ORDER and set(plan['actions'])==set(ORDER)
for action_label in ORDER:
 action=plan['actions'][action_label]
 assert set(action)=={'kind','test','expected_exit','owner_seconds','bridge_outer_seconds'}
 assert action['kind']=='Control' and action['test']==ACTION_NAMES[action_label]
 assert type(action['expected_exit']) is int and action['expected_exit']==0
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
record={'private_shader_experiment':True,'CI_shader_qualification_claimed':False,'mode':mode,'label':label,'case':case,'read_case_pass':True,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(WORK/'terminal-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'PRIVATE SHADER EXPERIMENT with unchanged original native-AMD fp64 joint-map assertions and eae host; bounded Windows birth/provider/accessible-consumer audit. No hidden/global/unconditional compiler-descendant absence or full qualification claim.'}
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
result={'pass':True,'case':case,'mode':mode,'label':label,'record':observed,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'loaded_provider_evidence_complete':record['loaded_provider_evidence_complete'],'source_revision':bindings['source_revision'],'native_build_gate_authority_used':True,'native_gate_authority_scope':'Unchanged CI host and BASE manifest only; private coupled shader excluded.','private_shader_experiment':True,'CI_shader_qualification_claimed':False,'mandatory_or_full_runtime_admitted':False,'scope':'PRIVATE SHADER EXPERIMENT running one unchanged original native-AMD fp64 joint-map control: actual_four_column_joint_rk4_and_independent_scalar_neighbours and coupled_supports_fp64=true, stagedSDKloader+AMD, complete zero-FESD XML and original assertions. Unsupported FP64 decline is rejected. Limited accessible census; no full trace/source-map/fullscience/Mandatory/performance/release claim.'}
assert not controller_signals,('controller signals before acceptance publication',controller_signals)
output=WORK/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
