"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path
import json,os,signal,subprocess,sys,time
from native_contracts import ROOT,WORK,read_case,document,sha
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
mode,label,case=sys.argv[1:]
binding_path=WORK/('baseline-native-bindings.json' if mode=='baseline' else 'execution-bindings.json')
bindings=document(binding_path)
binding_sha256=sha(binding_path)
observed=read_case(mode,label,case,bindings)
prefix='terminal-'+mode+'-'+label+'-'+case
receipt=WORK/(prefix+'.json');bridge=WORK/(prefix+'-bridge.json')
temp=WORK/'temp'/prefix;assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir(parents=True)
def win(path):return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()
def quote(word):return "'"+word.replace("'","''")+"'"
audit=win(WORK/'baseline-terminal-audit.ps1');workspace=win(WORK);bootstrap=win(temp/'bootstrap.json');run=mode+'/'+label+'/'+case
command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(bootstrap)+" -Encoding UTF8;& "+quote(audit)+" -Workspace "+quote(workspace)+" -RunName "+quote('native/'+run)+" -Receipt "+quote(receipt.name)
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(WORK/(prefix+'-emergency-cleanup.json'))]
record={'mode':mode,'label':label,'case':case,'read_case_pass':True,'execution_bindings_sha256':binding_sha256,'whole_input_source_seals_checked_before_audit':True,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(WORK/'baseline-terminal-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Actual normal owner completion, exact scientific XML and whole-input/clean-source seals checked before and after bounded read-only Windows birth/provider/accessible-consumer audit. Hidden fields/handles/global or unconditional descendant absence not asserted.'}
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
  assert sha(binding_path)==binding_sha256
  assert read_case(mode,label,case,bindings)==observed
  assert sha(binding_path)==binding_sha256
  record['whole_input_source_seals_checked_after_audit']=True
  record['terminal_receipt_sha256']=sha(receipt)
  record['bootstrap']=document(temp/'bootstrap.json')
 except Exception as error:
  record['accepted']=False;record['errors'].append(repr(error))
 # A completed interop bridge normally establishes auditor exit. Keep its exact
 # birth for the next/final census, without claiming hidden consumer ownership.
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
