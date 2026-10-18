"""ROOT-only negative shutdown audit AFTER the original owner helper is closed nonpass.
Uses the exact bound terminal/cleanup argv; never reads positive frame evidence or launches rendering.
"""
from pathlib import Path
import hashlib,json,os,signal,subprocess,sys,time
sys.dont_write_bytecode=True
assert not sys.flags.optimize and len(sys.argv)==1
WORK=Path(__file__).resolve().parent;ROOT=Path.cwd();assert ROOT==WORK.parents[1]
def document(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def identity(path):return {'bytes':path.stat().st_size,'sha256':sha(path)}
binding_path=WORK/'execution-bindings.json';bindings=document(binding_path);binding_sha256=sha(binding_path)
plan=document(WORK/'plan.json');case=label='native-frame';mode='native-frame-negative-terminal-audit'
assert bindings['admitted'] is True and plan['admitted'] is True
assert bindings['source_revision']==bindings['runner_source_revision']==plan['source_revision']==plan['runner_source_revision']=='2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert plan['order']==['native-frame'] and plan['terminal_observer_seconds']==60
runtime_path=WORK/(case+'-bridge.json');owner_path=WORK/'native'/case/'owner.json'
runtime=document(runtime_path);owner=document(owner_path)
assert type(runtime['returncode']) is int and runtime['returncode']!=0
assert runtime['binding_sha256']==binding_sha256 and runtime['whole_inputs_unchanged'] is True
assert owner['numerical_control_passed'] is False and owner['owned_process_absent'] is True
assert owner['end_utc'] and owner['source_revision']==bindings['source_revision'] and owner['live_source_revision']==bindings['runner_source_revision']
assert owner['source_build_gate_sha256']==bindings['gate_sha256'] and owner['source_build_gate_unchanged'] and owner['test_executable_unchanged']
before=document(WORK/(case+'-inputs-before.json'));assert before==document(WORK/(case+'-inputs-after.json'))
raw=[runtime_path,owner_path,WORK/(case+'-inputs-before.json'),WORK/(case+'-inputs-after.json')]
for name in ['stdout.log','stderr.log','memory.json','driver-modules.json']:
 path=WORK/'native'/case/name
 if path.exists():raw.append(path)
snapshot={str(path.relative_to(ROOT)):identity(path) for path in raw}
script_identity=identity(Path(__file__));observer_identity=identity(WORK/'run_terminal.py')
def verify():
 assert sha(binding_path)==binding_sha256
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==bindings['runner_source_revision']
 assert not subprocess.check_output(['git','status','--porcelain'])
 for mapping in (bindings['files'],before,snapshot):
  for name,expected in mapping.items():
   path=Path(name) if name.startswith('/') else ROOT/name
   assert identity(path)==expected,name
 assert identity(Path(__file__))==script_identity and identity(WORK/'run_terminal.py')==observer_identity
verify()
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
prefix='terminal-'+case;receipt=WORK/(prefix+'.json');temp=WORK/'temp'/prefix
bridge=WORK/(prefix+'-failure-audit-bridge.json')
assert not temp.exists() and not receipt.exists() and not bridge.exists()
assert not (WORK/(prefix+'.stdout')).exists() and not (WORK/(prefix+'.stderr')).exists()
arguments=bindings['terminal_actions'][case]['argv'];cleanup_args=bindings['terminal_actions'][case]['cleanup_argv']
assert arguments[:6]==['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command']
assert len(arguments)==7 and 'terminal-audit.ps1' in arguments[-1] and 'run_native_control.ps1' not in arguments[-1]
assert 'emergency_native_cleanup.ps1' in ' '.join(cleanup_args)
temp.mkdir(parents=True)
record={'schema':'sirius-native-constant-frame-nonpass-readonly-terminal-audit-v1','mode':mode,'label':label,'case':case,'numerical_admitted':False,'numerical_passed':False,'failure_audit_only':True,'read_case_pass':False,'observer_sha256':script_identity['sha256'],'positive_observer_source':observer_identity,'audit_script_sha256':sha(WORK/'terminal-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'controller_pid':os.getpid(),'controller_birth':Path('/proc',str(os.getpid()),'stat').read_text().rsplit(')',1)[1].split()[19],'binding':identity(binding_path),'closed_nonpass_bridge':identity(runtime_path),'qualification_claimed':False}
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
  terminal=document(receipt)
  assert terminal['owned_births_absent'] and terminal['all_observed_module_seals_match']
  assert terminal['actual_loaded_modules_observed'] is True
  assert terminal['consumer_scan']['matches']==terminal['consumer_scan']['accessible_compilers']==[]
  assert len(terminal['births'])==2 and {b['kind'] for b in terminal['births']}=={'owner','test'}
  for birth in terminal['births']:
   assert birth['pid']==owner[birth['kind']+'_pid'] and birth['recorded_start_utc']==owner[birth['kind']+'_start_time_utc'] and birth['owned_birth_absent']
  modules=document(WORK/'native'/case/'driver-modules.json');assert len(modules)==2
  assert [m for m in modules if 'amdvlk' in m['name'].lower()]==[bindings['expected_amd_module']]
  loader,=[m for m in modules if m['name'].lower()=='vulkan-1.dll']
  relative='src/sirius/app/Release/vulkan-1.dll';expected=bindings['files'][bindings['stage']+'/'+relative]
  assert {k:loader[k] for k in ['bytes','sha256']}==expected
  assert loader['path'].replace('/','\\').casefold()==(bindings['windows_checkout_unc']+'\\'+(bindings['stage']+'/'+relative).replace('/','\\')).casefold()
  assert len(terminal['loaded_module_post_seals'])==2
  for seal in terminal['loaded_module_post_seals']:
   module,=[m for m in modules if m['path']==seal['path']]
   assert seal['unchanged'] and seal['expected_sha256']==seal['post_sha256']==module['sha256'] and seal['expected_bytes']==seal['post_bytes']==module['bytes']
  assert terminal['owned_executable_post_seal']['unchanged'] and terminal['owned_executable_post_seal']['sha256']==owner['test_executable_sha256']
  verify();assert not controller_signals
  record['terminal_receipt']=identity(receipt);record['bootstrap']=document(temp/'bootstrap.json')
  record['source_binding_owner_bridge_action_post_sealed']=True
 except BaseException as error:
  record['accepted']=False;record['errors'].append('post-audit: '+repr(error))
if controller_signals:
 record['outer_stop']=True;record['accepted']=False
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
record['accepted_failure_shutdown_audit']=record['accepted']
record['scope']='Closed original nonpass teardown only;exactbound60sPS audit/emergency/birth/module/input postseals and accessible-census scope. No numerical success,renderer launch,retry,source cleanup or qualification authority.'
with bridge.open('x',encoding='utf-8') as stream:json.dump(record,stream,indent=2);stream.write('\n')
print(json.dumps({k:record[k] for k in ['returncode','elapsed_seconds','accepted_failure_shutdown_audit','numerical_admitted']}),flush=True)
assert record['accepted_failure_shutdown_audit']
