"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path
import json,os,signal,subprocess,sys,time
sys.path.insert(0,str(Path.cwd()/'out/transport-paired-row-integration'))
from native_contracts import ROOT,WORK,read_case,document,sha
REST=Path(__file__).resolve().parent
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
mode,label,case='matched','a2','timestamps'
binding_path=WORK/('baseline-native-bindings.json' if mode=='baseline' else 'execution-bindings.json')
bindings=document(binding_path)
observed=read_case(mode,label,case,bindings)
prefix='final-native-cleanup'
receipt=WORK/(prefix+'.json');bridge=WORK/(prefix+'-bridge.json')
temp=REST/'temp'/prefix;assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir(parents=True)
def win(path):return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()
def quote(word):return "'"+word.replace("'","''")+"'"
audit=win(REST/'final_native_terminal.ps1');workspace=win(WORK);bootstrap=win(temp/'bootstrap.json');run=mode+'/'+label+'/'+case
command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(bootstrap)+" -Encoding UTF8;& "+quote(audit)+" -Workspace "+quote(workspace)+" -BaselineStage "+quote(win(ROOT/bindings['producers']['baseline']['stage']))+" -CandidateStage "+quote(win(ROOT/bindings['producers']['candidate']['stage']))+" -Receipt "+quote(win(receipt))
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(WORK/(prefix+'-emergency-cleanup.json'))]
record={'mode':mode,'label':label,'case':case,'read_case_pass':True,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(REST/'final_native_terminal.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Historical accepted latest control receipts plus fresh all45Windowsbirth/provider/product/accessiblereference audit beforefinishedtrialcleanup. Hiddenfields/handles/global/unconditionaldescendantabsence not asserted.'}

head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert head==document(REST/'source.json')['source_revision'] and not subprocess.check_output(['git','status','--porcelain'])
protection=document(WORK/'archive-protection.json');archive=ROOT/protection['archive'];manifest=document(archive/'manifest.json')
assert sha(archive/'manifest.json')==protection['manifest']['sha256']
seal_files={}
for r in manifest['payloads']:
 p=ROOT/r['origin']
 if p.is_relative_to(WORK/'native') or p.is_relative_to(WORK/'temp') or p.is_relative_to(ROOT/'bin/windows-msvc'):
  assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'];seal_files[str(p)]={k:r[k] for k in ('bytes','sha256')}
for path,expected in bindings['files'].items():
 if path.startswith('/mnt/c/'):
  p=Path(path);assert p.stat().st_size==expected['bytes'] and sha(p)==expected['sha256']
parser=WORK/'native_contracts.py'
parser_origin=str(parser.relative_to(ROOT))
expected=next(q for q in manifest['payloads'] if q['origin']==parser_origin)
assert parser.stat().st_size==expected['bytes'] and sha(parser)==expected['sha256']
for p in [Path(__file__),REST/'final_native_terminal.ps1',REST/'source.json',WORK/'archive-protection.json',archive/'manifest.json',parser,WORK/'execution-bindings.json',WORK/'emergency_native_cleanup.ps1',*map(Path,[q for q in bindings['files'] if q.startswith('/mnt/c/')])]:
 seal_files[str(p)]={'bytes':p.stat().st_size,'sha256':sha(p)}
def verify_seals():
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head and not subprocess.check_output(['git','status','--porcelain'])
 for path,r in seal_files.items():assert Path(path).stat().st_size==r['bytes'] and sha(Path(path))==r['sha256'],path
verify_seals()

record['input_seals']=seal_files;record['restored_source_revision']=head
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
 terminal=document(receipt);assert terminal['pass'] and len(terminal['births'])==45 and all(x['owned_birth_absent'] for x in terminal['births']) and len(terminal['post_seals'])==4
 assert not terminal['census']['matches'] and not terminal['census']['accessible_compilers']
 record['terminal_receipt_sha256']=sha(receipt)
 record['bootstrap']=document(temp/'bootstrap.json')
 # A completed interop bridge normally establishes auditor exit. Keep its exact
 # birth for the next/final census, without claiming hidden consumer ownership.
if controller_signals:
 record['outer_stop']=True;record['accepted']=False
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
verify_seals();record['source_and_input_seals_exact_at_end']=True
if controller_signals:
 record['outer_stop']=True;record['accepted']=False
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
bridge.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('mode','label','case','returncode','elapsed_seconds','accepted')}),flush=True)
assert record['accepted']
result={'pass':True,'case':case,'mode':mode,'label':label,'record':observed,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'scope':'Fresh final45recordedbirth/provider/product/accessible-reference audit for historical finished native trial and stages; no new scientific execution, performance or qualification claim.'}
assert not controller_signals,('controller signals before acceptance publication',controller_signals)
output=WORK/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
