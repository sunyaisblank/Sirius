"""Bounded read-only accessible consumer audit for three obsolete native stages."""
from pathlib import Path
import json,os,signal,subprocess,sys,time,hashlib
sys.dont_write_bytecode=True
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def document(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def identity(p):return {'bytes':p.stat().st_size,'sha256':sha(p)}
plan=document(WORK/'plan.json')
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==plan['source_revision']
assert not subprocess.check_output(['git','status','--porcelain'])
for entries in plan['files'].values():
 for path,item in entries.items():
  assert identity(ROOT/path)==item['original']
  r=item['retained_official_export'];assert identity(ROOT/r['path'])=={k:r[k] for k in ['bytes','sha256']}
def win(path):return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()
def quote(word):return "'"+word.replace("'","''")+"'"
roots=[win(ROOT/rel) for rel in plan['cleanup_roots']]
(WORK/'windows-roots.json').write_text(json.dumps(roots,indent=2)+'\n')
inputs=[Path(__file__),WORK/'consumer-audit.ps1',WORK/'emergency_native_cleanup.ps1',WORK/'plan.json',WORK/'windows-roots.json']
input_seals={str(p):identity(p) for p in inputs}
receipt=WORK/'windows-consumers.json';bridge=WORK/'observer-bridge.json';temp=WORK/'temp'
assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir()
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'consumer-audit.ps1'),'-Workspace',win(WORK)]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(WORK/'observer-emergency-cleanup.json')]
record={'kind':'obsolete native stage accessible consumer audit','deadline_seconds':60,'outer_stop':False,'errors':[],'scientific_pass':False,'qualification_claimed':False,'input_seals':input_seals}
prefix='observer'
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted);signal.signal(signal.SIGHUP,interrupted)
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
started=time.monotonic();child=None;code=None;record['spawn_attempted']=False
with (WORK/(prefix+'.stdout')).open('wb') as out,(WORK/(prefix+'.stderr')).open('wb') as err:
 assert not controller_signals,('interrupted before spawn',controller_signals)
 try:
  record['spawn_attempted']=True
  child=subprocess.Popen(arguments,stdout=out,stderr=err,start_new_session=True)
  record['linux_child_pid']=child.pid;record['linux_child_birth']=birth(child.pid)
  assert record['linux_child_birth'] is not None,'Owned Linux birth unavailable'
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
  if record['spawn_attempted'] and (record['outer_stop'] or child is None or child.poll() is None):
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
     expected=record.get('linux_child_birth')
     assert expected is not None and birth(child.pid)==expected,'Owned birth unavailable or changed; no group signal'
     assert os.getpgid(child.pid)==child.pid,'Owned group changed; no group signal'
     os.killpg(child.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    except BaseException as error:record['errors'].append(repr(error))
    try:code=child.wait(timeout=10)
    except BaseException as error:record['errors'].append(repr(error));code=child.returncode

if controller_signals:
 record['outer_stop']=True
 if 'controller_signal' not in record['errors']:record['errors'].append('controller_signal')
record['controller_signals']=list(controller_signals)
record['returncode']=code;record['elapsed_seconds']=time.monotonic()-started
current=birth(child.pid) if child is not None else None
record['linux_child_birth_absent']=record.get('linux_child_birth') is not None and current!=record['linux_child_birth']
record['accepted']=code==0 and not record['outer_stop'] and not record['errors'] and record['linux_child_birth_absent']
if record['accepted']:
 try:
  census=document(receipt);assert census['pass'] and not census['matches']
  assert census['roots']==roots
  assert {str(p):identity(p) for p in inputs}==input_seals
  for entries in plan['files'].values():
   for path,item in entries.items():
    assert identity(ROOT/path)==item['original']
    r=item['retained_official_export'];assert identity(ROOT/r['path'])=={k:r[k] for k in ['bytes','sha256']}
  record['bootstrap']=document(temp/'bootstrap.json')
  record['receipt_sha256']=sha(receipt)
  record['inputs_and_retention_post_sealed']=True
 except BaseException as error:record['accepted']=False;record['errors'].append(repr(error))
if controller_signals:
 record['outer_stop']=True;record['accepted']=False
record['controller_signals']=list(controller_signals)
bridge.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ['returncode','elapsed_seconds','accepted']}),flush=True)
assert record['accepted'] and not controller_signals
