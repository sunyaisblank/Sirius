"""Bounded help inventory with historical exact-birth emergency recovery."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def win(p):return subprocess.check_output(['wslpath','-w',str(p.resolve())],text=True).strip()
def git(*a):return subprocess.check_output(['git',*a],text=True).strip()
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
assert len(sys.argv)==2 and sys.argv[1] in ('profiler-trace-ui',)
action=sys.argv[1];head=read(WORK/'source.json')['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
review=read(WORK/'trace-ui-source-review.json');assert review['accepted']
extraction=read(WORK/'tools-extraction.json')
paths=[Path(__file__),WORK/'profiler_trace_ui.ps1',WORK/'emergency_native_cleanup.ps1',WORK/'source.json',WORK/'tools-extraction.json',WORK/'trace-ui-source-review.json']
gui=read(WORK/'tools-gui-extraction.json');paths += [WORK/'tools-gui-extraction.json']
trace=WORK/'native-capture-early-console/original-legacy-te.rgp';acquisition=read(WORK/'native-capture-early-console-bridge.json');assert acquisition['accepted']
assert trace.stat().st_size==acquisition['original_control_and_acquired_trace']['trace_bytes'] and sha(trace)==acquisition['original_control_and_acquired_trace']['trace_sha256']
paths += [trace,WORK/'native-capture-early-console-bridge.json',WORK/'prepare_trace_ui.py']
paths += [ROOT/v['path'] for v in extraction['files']+gui['files']]
seals={str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in paths}
for name,v in review['controller_seals'].items():assert seals[str((WORK/name).relative_to(ROOT))]==v
for v in extraction['files']+gui['files']:assert seals[v['path']]=={k:v[k] for k in ('bytes','sha256')}
tool='RadeonGPUProfiler.exe'
matches=[ROOT/v['path'] for v in extraction['files'] if Path(v['path']).name==tool];assert len(matches)==1
exe=matches[0];tools=exe.parent;output=WORK/action;assert not output.exists()
prefix=WORK/(action+'-bridge');assert not prefix.with_suffix('.json').exists()
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
args=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'profiler_trace_ui.ps1'),'-Tools',win(tools),'-Output',win(output),'-Action',action,'-ExpectedHash',sha(exe),'-LiveRevision',head,'-Trace',win(trace),'-TraceHash',sha(trace)]
def quote(v):return "'"+v.replace("'","''")+"'"
cleanup_bootstrap=WORK/(action+'-cleanup-helper-bootstrap.json')
cleanup_command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');scope='Exact emergency helper observer birth before scoped recovery'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(win(cleanup_bootstrap))+" -Encoding UTF8;& "+quote(win(WORK/'emergency_native_cleanup.ps1'))+" -TemporaryDirectory "+quote(win(output/'temp'))+" -Output "+quote(win(output))+" -Receipt "+quote(win(WORK/(action+'-emergency.json')))
cleanup=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',cleanup_command]
signals=[]
for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(signum,lambda signum,frame:signals.append(signum))
record={'action':action,'source_revision':head,'linux_owner_pid':os.getpid(),'linux_owner_start_ticks':birth(os.getpid()),'input_seals':seals,'arguments':args,'deadline_seconds':90,'errors':[],'outer_stop':False,'scope':'Owned GUI window observation only. Intentional owned stop after observation; no service/science/device capture; recorded tool/controller ownership and finite accessible developer-process check.'}
def reap_owned(c,ticks,timeout):
 try:return c.wait(timeout=timeout)
 except subprocess.TimeoutExpired:
  # Signal only the saved birth in its independently created process group.
  try:
   assert ticks and birth(c.pid)==ticks and os.getpgid(c.pid)==c.pid
   os.killpg(c.pid,signal.SIGKILL)
  except (FileNotFoundError,ProcessLookupError):pass
  return c.wait(timeout=10)
def stop_cleanup_observer():
 # The leaf observer records its own birth before targeting the saved helper.
 # A leaf timeout is retained as unresolved ownership, never accepted.
 leaf_bootstrap=WORK/(action+'-cleanup-stop-bootstrap.json')
 leaf_receipt=WORK/(action+'-cleanup-stop.json')
 command="$ErrorActionPreference='Stop';$self=Get-Process -Id $PID;$selfHandle=$self.Handle;@{owner_pid=$PID;owner_start_time_utc=$self.StartTime.ToUniversalTime().ToString('o')}|ConvertTo-Json|Set-Content -LiteralPath "+quote(win(leaf_bootstrap))+" -Encoding UTF8;$b=Get-Content -Raw -LiteralPath "+quote(win(cleanup_bootstrap))+"|ConvertFrom-Json;$p=Get-Process -Id ([int]$b.owner_pid) -ErrorAction SilentlyContinue;if($p){$h=$p.Handle;if($p.StartTime.ToUniversalTime().ToString('o') -eq $b.owner_start_time_utc){$p.Kill();if(!$p.WaitForExit(5000)){throw 'Emergency helper did not stop'}}};$p=Get-Process -Id ([int]$b.owner_pid) -ErrorAction SilentlyContinue;$absent=(!$p -or $p.StartTime.ToUniversalTime().ToString('o') -ne $b.owner_start_time_utc);@{pass=$absent;target=$b}|ConvertTo-Json -Depth 6|Set-Content -LiteralPath "+quote(win(leaf_receipt))+" -Encoding UTF8;if(!$absent){exit 1}"
 args=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
 c=None;ticks=None
 try:
  with (WORK/(action+'-cleanup-stop.stdout')).open('xb') as out,(WORK/(action+'-cleanup-stop.stderr')).open('xb') as err:
   c=subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True);ticks=birth(c.pid)
   record['cleanup_stop_observer']={'pid':c.pid,'start_ticks':ticks,'bootstrap':str(leaf_bootstrap.relative_to(ROOT)),'arguments':args}
   code=c.wait(timeout=20);assert code==0 and read(leaf_receipt)['pass']
 finally:
  if c is not None:
   code=reap_owned(c,ticks,5)
   record['cleanup_stop_observer'].update(returncode=code,linux_birth_absent=birth(c.pid)!=ticks)
def recover_windows():
 c=None;ticks=None
 try:
  with (WORK/(action+'-emergency.stdout')).open('xb') as out,(WORK/(action+'-emergency.stderr')).open('xb') as err:
   c=subprocess.Popen(cleanup,stdout=out,stderr=err,start_new_session=True);ticks=birth(c.pid)
   record['emergency_observer']={'pid':c.pid,'start_ticks':ticks,'bootstrap':str(cleanup_bootstrap.relative_to(ROOT)),'arguments':cleanup}
   try:code=c.wait(timeout=60)
   except subprocess.TimeoutExpired:
    record['errors'].append('Emergency helper exceeded60s; tool ownership remains invalidated')
    assert cleanup_bootstrap.exists(),'Emergency helper exact Windows birth unavailable'
    stop_cleanup_observer()
    raise
  assert code==0
  receipt=read(WORK/(action+'-emergency.json'));assert receipt['recorded_births_absent'] and not receipt['errors']
 finally:
  if c is not None:
   code=reap_owned(c,ticks,5)
   record['emergency_observer'].update(returncode=code,linux_birth_absent=birth(c.pid)!=ticks)
child=None;code=None;attempted=False;start=time.monotonic()
try:
 assert not signals
 with prefix.with_suffix('.stdout').open('xb') as out,prefix.with_suffix('.stderr').open('xb') as err:
  attempted=True;child=subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True)
  record['linux_child_pid']=child.pid;record['linux_child_start_ticks']=birth(child.pid)
  while child.poll() is None:
   if signals:raise RuntimeError(('controller_signals',signals))
   remaining=90-(time.monotonic()-start)
   if remaining<=0:raise subprocess.TimeoutExpired(args,90)
   try:child.wait(timeout=min(.25,remaining))
   except subprocess.TimeoutExpired:pass
  code=child.returncode
except BaseException as e:record['errors'].append(repr(e));record['outer_stop']=True
finally:
 if signals:record['outer_stop']=True
 if attempted and (record['outer_stop'] or child is None or child.poll() is None or child.returncode!=0):
  try:
   recover_windows()
  except BaseException as e:record['errors'].append(repr(e))
 if child is not None:
  try:code=child.wait(timeout=15)
  except BaseException as e:
   record['errors'].append(repr(e))
   try:
    saved=record.get('linux_child_start_ticks');assert saved and birth(child.pid)==saved and os.getpgid(child.pid)==child.pid
    os.killpg(child.pid,signal.SIGKILL)
   except (FileNotFoundError,ProcessLookupError):pass
   except BaseException as e:record['errors'].append(repr(e))
   try:code=child.wait(timeout=10)
   except BaseException as e:record['errors'].append(repr(e))
record['controller_signals']=signals;record['returncode']=code;record['elapsed_seconds']=time.monotonic()-start
record['linux_child_birth_absent']=child is not None and birth(child.pid)!=record['linux_child_start_ticks']
record['accepted']=code==0 and not record['outer_stop'] and not record['errors'] and not signals and record['linux_child_birth_absent']
if record['accepted']:
 try:
  o=read(output/'owner.json');assert o['status']=='observed-and-stopped' and o['intentional_owned_stop'] and o['window_observation'] and not o['outer_stop'] and not o['cleanup_errors'] and o['owned_process_absent'] and o['tool_unchanged'] and o['developer_census_completed'] and not o['accessible_remaining_developer_processes']
  assert o['test_executable_sha256']==sha(exe) and o['trace_unchanged'] and o['trace_sha256']==sha(trace)
  saved=[v for v in o['window_observation'] if v['print_window_succeeded'] and v['image_path']]
  assert saved
  for v in saved:
   p=Path(subprocess.check_output(['wslpath','-u',v['image_path']],text=True).strip())
   assert p.resolve().parent==output.resolve() and p.stat().st_size==v['image_bytes'] and sha(p)==v['image_sha256']
  assert read(output/'signature.json')['status']=='Valid'
  assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
  assert {str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in paths}==seals
  record['Windows_owner']=o
  record['stdout_sha256']=sha(output/'stdout.log');record['stderr_sha256']=sha(output/'stderr.log')
 except BaseException as e:
  record['accepted']=False;record['errors'].append(repr(e))
  try:recover_windows()
  except BaseException as e:record['errors'].append(repr(e))
if signals:record['outer_stop']=True;record['accepted']=False
record['controller_signals']=list(signals)
prefix.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('action','returncode','elapsed_seconds','accepted')}),flush=True)
assert record['accepted']
