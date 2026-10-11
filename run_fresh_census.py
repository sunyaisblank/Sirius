"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path
import json,os,signal,subprocess,sys,time
from native_contracts import ROOT,WORK,document,sha
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):
 signal.signal(signum,interrupted)
assert len(sys.argv)==1
source=document(WORK/'restoration/source.json');head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{u}'])
assert document(WORK/'native-preservation-reviewed.json')['accepted']
assert document(WORK/'software-preservation-reviewed.json')['accepted']
assert document(WORK/'native-terminal-independent-review.json')['pass']
def cleanup_input_seals():
 seals={}; native=None
 for name,review_name in (('native-protection.json','native-preservation-reviewed.json'),('software-protection.json','software-preservation-reviewed.json')):
  receipt_value=document(WORK/name);review=document(WORK/review_name)
  assert review['accepted'] and review['manifest_sha256']==receipt_value['manifest']['sha256']
  archive=ROOT/receipt_value['archive'];manifest_path=archive/'manifest.json'
  assert manifest_path.stat().st_size==receipt_value['manifest']['bytes'] and sha(manifest_path)==review['manifest_sha256']
  for path in (WORK/name,WORK/review_name,manifest_path):
   seals[str(path.relative_to(ROOT))]={'bytes':path.stat().st_size,'sha256':sha(path)}
  if name=='native-protection.json':native=document(manifest_path)
 records={v['origin']:v for v in native['diagnostics']}
 owners=sorted((WORK/'native').rglob('owner.json'))
 auditors=sorted(p for p in (WORK/'temp').rglob('bootstrap.json') if p.relative_to(WORK/'temp').parts[0].startswith('terminal-'))
 assert len(owners)==12 and len(auditors)==12
 births=set()
 for path in [*owners,*auditors]:
  origin=str(path.relative_to(ROOT));value={'bytes':path.stat().st_size,'sha256':sha(path)}
  assert value=={k:records[origin][k] for k in ('bytes','sha256')},origin
  seals[origin]=value;data=document(path)
  for kind in (('owner','test') if path in owners else ('owner',)):
   number=data[kind+'_pid'];birth=data[kind+'_start_time_utc']
   assert type(number) is int and number>0 and birth
   assert (number,birth) not in births
   births.add((number,birth))
 assert len(births)==36
 return seals,sorted(births)
input_seals,expected_births=cleanup_input_seals()
mode,label,case='finished','restored','trial'
binding_path=WORK/'execution-bindings.json'
bindings=document(binding_path)
binding_sha256=sha(binding_path)
observer_sha256=sha(Path(__file__))
audit_sha256=sha(WORK/'fresh_native_consumers.ps1')
prefix='fresh-native-consumers'
receipt=WORK/(prefix+'.json');bridge=WORK/(prefix+'-bridge.json')
temp=WORK/'temp'/prefix;assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir(parents=True)
def win(path):return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()
def quote(word):return "'"+word.replace("'","''")+"'"
audit=win(WORK/'fresh_native_consumers.ps1');workspace=win(WORK);bootstrap=win(temp/'bootstrap.json');run=mode+'/'+label+'/'+case
command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(bootstrap)+" -Encoding UTF8;& "+quote(audit)+" -Trial "+quote(workspace)+" -Receipt "+quote(win(receipt))
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(WORK/(prefix+'-emergency-cleanup.json'))]
record={'mode':mode,'label':label,'case':case,'historical_recovery_review_pass':True,'restoration_revision':head,'linux_owner_pid':os.getpid(),'linux_owner_start_ticks':Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],'execution_bindings_sha256':binding_sha256,'historical_archive_manifests_checked_before_audit':True,'consumed_input_seals':input_seals,'expected_Windows_births':expected_births,'observer_sha256':observer_sha256,'audit_script_sha256':audit_sha256,'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Restored clean pushed source and historical archive manifests checked before bounded read-only all recorded Windows birth and accessible scoped-consumer audit. Hidden fields/handles/global or unconditional descendant absence not asserted.'}
started=time.monotonic()
child=None;code=None;spawn_attempted=False
with (WORK/(prefix+'.stdout')).open('wb') as out,(WORK/(prefix+'.stderr')).open('wb') as err:
 try:
  assert not controller_signals,('interrupted before spawn',controller_signals)
  spawn_attempted=True;record['spawn_attempted']=True
  child=subprocess.Popen(arguments,stdout=out,stderr=err,start_new_session=True)
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
  if spawn_attempted and (record['outer_stop'] or child is None or child.poll() is None):
   try:
    with (WORK/(prefix+'-emergency.stdout')).open('wb') as cout,(WORK/(prefix+'-emergency.stderr')).open('wb') as cerr:
     c=subprocess.run(cleanup_args,stdout=cout,stderr=cerr,timeout=60)
    assert c.returncode==0
    cleanup=document(WORK/(prefix+'-emergency-cleanup.json'));assert cleanup['recorded_births_absent'] and not cleanup['errors']
   except BaseException as error:record['errors'].append(repr(error))
  if child is not None:
   try:code=child.wait(timeout=15)
   except BaseException as error:
    record['errors'].append(repr(error))
    try:
     saved=record.get('linux_child_birth')
     current=Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
     assert saved is not None and current==saved and os.getpgid(child.pid)==child.pid, 'Unconfirmed Linux birth/process group; signal declined'
     os.killpg(child.pid,signal.SIGKILL)
    except (FileNotFoundError,ProcessLookupError):pass
    except BaseException as error:record['errors'].append(repr(error))
    # Reap independently even if identity verification or signalling failed.
    try:code=child.wait(timeout=10)
    except BaseException as error:record['errors'].append(repr(error));code=child.returncode
if controller_signals:
 record['outer_stop']=True
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
record['returncode']=code;record['elapsed_seconds']=time.monotonic()-started
current=None
if child is not None:
 try:current=Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):current=None
record['linux_child_birth_absent']=record.get('linux_child_birth') is not None and current!=record['linux_child_birth']
record['accepted']=code==0 and not record['outer_stop'] and not record['errors'] and record['linux_child_birth_absent']
if record['accepted']:
 try:
  terminal=document(receipt);assert terminal['pass'] and len(terminal['known_births'])==36 and all(v['owned_birth_absent'] for v in terminal['known_births'])
  assert sorted((v['pid'],v['recorded_start_utc']) for v in terminal['known_births'])==expected_births
  assert not terminal['matches'] and set(terminal['targets'])=={'fma-dense-dp-review','native-7341f20','native-b728680'}
  assert sha(binding_path)==binding_sha256
  assert sha(Path(__file__))==observer_sha256 and sha(WORK/'fresh_native_consumers.ps1')==audit_sha256
  assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head and not subprocess.check_output(['git','status','--porcelain'])
  assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{u}'])
  assert cleanup_input_seals()==(input_seals,expected_births)
  assert sha(binding_path)==binding_sha256
  record['restoration_source_and_controller_seals_checked_after_audit']=True
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
result={'pass':True,'case':case,'mode':mode,'label':label,'restoration_revision':head,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'scope':'Fresh36 recorded native/terminal births and accessible scoped-consumer census, restored-source and controller continuity, historical archived inputs; own observer normal interop completion. No deletion or broader qualification claim.'}
assert not controller_signals,('controller signals before acceptance publication',controller_signals)
output=WORK/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
