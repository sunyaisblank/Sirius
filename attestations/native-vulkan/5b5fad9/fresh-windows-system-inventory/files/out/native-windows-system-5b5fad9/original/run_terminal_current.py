"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path
import json,os,signal,subprocess,sys,time
import hashlib
sys.dont_write_bytecode=True
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def document(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_case(mode,label,case,bindings):
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==bindings['source_revision']
 assert not subprocess.check_output(['git','status','--porcelain'])
 bridge=document(WORK/(case+'-bridge.json'));owner=document(WORK/'native'/case/'owner.json')
 assert bridge['accepted_owner_execution'] and bridge['completed'] and bridge['returncode']==0
 assert bridge['binding_sha256']==sha(WORK/'execution-bindings.json')
 assert owner['status']=='completed' and owner['returncode']==0 and not owner['outer_stop'] and not owner['cleanup_errors']
 assert owner['owned_process_absent'] and owner['numerical_control_passed'] and owner['loaded_driver_modules_recorded']
 assert owner['source_revision']==owner['live_source_revision']==bindings['source_revision']
 for rel,item in bindings['files'].items():
  p=Path(rel) if rel.startswith('/') else ROOT/rel
  assert {'bytes':p.stat().st_size,'sha256':sha(p)}==item,rel
 receipts=[WORK/(case+'-bridge.json'),WORK/'native'/case/'owner.json',WORK/'native'/case/'gtest.xml',WORK/'native'/case/'driver-modules.json']
 whole_inputs={str(p.relative_to(WORK)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in receipts}
 return {'owner':owner,'bridge':bridge,'gtest_sha256':sha(WORK/'native'/case/'gtest.xml'),'whole_inputs':whole_inputs}
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
case,=sys.argv[1:];assert case in ('science','render');mode='current';label=case
binding_path=WORK/'execution-bindings.json'
bindings=document(binding_path);binding_sha256=sha(binding_path)
observed=read_case(mode,label,case,bindings)
prefix='terminal-'+case
receipt=WORK/(prefix+'.json');bridge=WORK/(prefix+'-bridge.json')
temp=WORK/'temp'/prefix;assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir(parents=True)
def win(path):return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()
def quote(word):return "'"+word.replace("'","''")+"'"
audit=win(WORK/'terminal-audit.ps1');workspace=win(WORK);bootstrap=win(temp/'bootstrap.json');run=case
command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(bootstrap)+" -Encoding UTF8;& "+quote(audit)+" -Workspace "+quote(workspace)+" -RunName "+quote('native/'+run)+" -Receipt "+quote(receipt.name)
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(WORK/(prefix+'-emergency-cleanup.json'))]
record={'mode':mode,'label':label,'case':case,'read_case_pass':True,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(WORK/'terminal-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Actual normal owner completion and exact scientific XML checked before bounded read-only Windows birth/provider/accessible-consumer audit. Hidden fields/handles/global or unconditional descendant absence not asserted.'}
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
  terminal=document(receipt);assert terminal['owned_births_absent'] and terminal['all_two_module_seals_match']
  assert not terminal['consumer_scan']['matches'] and not terminal['consumer_scan']['accessible_compilers']
  record['terminal_receipt_sha256']=sha(receipt)
  record['bootstrap']=document(temp/'bootstrap.json')
  assert sha(binding_path)==binding_sha256
  assert read_case(mode,label,case,bindings)==observed
  assert not controller_signals
  record['source_binding_owner_bridge_xml_post_sealed']=True
 except BaseException as error:
  record['accepted']=False;record['errors'].append('post-audit: '+repr(error))
if controller_signals:
 record['outer_stop']=True;record['accepted']=False
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
bridge.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('mode','label','case','returncode','elapsed_seconds','accepted')}),flush=True)
assert record['accepted']
result={'pass':True,'case':case,'mode':mode,'label':label,'record':observed,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'scope':'One finite numerical/native owner/provider acceptance with limited accessible terminal census, no performance/full-frame/full qualification claim.'}
assert not controller_signals,('controller signals before acceptance publication',controller_signals)
output=WORK/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
