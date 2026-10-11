"""Root-owned original numerical control with shader-clock sidecars; bounded bridge and scoped recovery."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();HERE=Path(__file__).resolve().parent
assert len(sys.argv)==2 and sys.argv[1] in ('device-timestamps','mixed-transport','factored-endpoints')
action=sys.argv[1];WORK=HERE/action
assert not WORK.exists(), 'Fresh per-case task working directory required'
WORK.mkdir()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def seal_key(p):return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p.resolve())
def win(p):return subprocess.check_output(['wslpath','-w',str(p.resolve())],text=True,timeout=10).strip()
def git(*a):return subprocess.check_output(['git',*a],text=True,timeout=20).strip()
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
plan=read(HERE/'execution-plan.json');head=plan['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
review=read(HERE/'controller-review.json');assert review['accepted']
paths=[Path(__file__),HERE/'execution-plan.json',HERE/'controller-review.json',HERE/'run_shader_clock_control.ps1',HERE/'emergency_native_cleanup.ps1',HERE/'run_native_retirement.py',HERE/'verify_native_retirement.ps1',HERE/'source/source-review-final.json']
paths += [Path(p) if Path(p).is_absolute() else ROOT/p for p in plan['input_seals']]
paths += [Path(p['linux_path']) for p in plan['providers']]
seals={seal_key(p):seal(p) for p in paths}
for name,v in review['controller_seals'].items():assert seal(HERE/name)==v
for name,v in plan['input_seals'].items():assert seals[name]==v
for p in plan['providers']:assert seal(Path(p['linux_path']))=={k:p[k] for k in ('bytes','sha256')}
output=WORK/'run';temporary=WORK/'owner-temp';prefix=WORK/'native-control-bridge'
assert not output.exists() and not temporary.exists() and not prefix.with_suffix('.json').exists()
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
filter=plan['cases'][action]
args=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(HERE/'run_shader_clock_control.ps1'),'-Stage',win(ROOT/plan['executable_directory']),'-Output',win(output),'-Filter',filter,'-Expected',filter,'-Seconds','320','-RssBytes','4294967296','-TemporaryDirectory',win(temporary),'-BuildRoot',win(ROOT/plan['build_root']),'-ExpectedGateSha256',plan['gate_sha256'],'-LiveRevision',head]
def quote(v):return "'"+v.replace("'","''")+"'"
cleanup_bootstrap=WORK/'cleanup-helper-bootstrap.json'
cleanup_command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');scope='Exact emergency helper observer birth before scoped recovery'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(win(cleanup_bootstrap))+" -Encoding UTF8;& "+quote(win(HERE/'emergency_native_cleanup.ps1'))+" -TemporaryDirectory "+quote(win(temporary))+" -Output "+quote(win(output))+" -Receipt "+quote(win(WORK/'emergency.json'))
cleanup=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',cleanup_command]
signals=[]
for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(signum,lambda signum,frame:signals.append(signum))
record={'action':action,'source_revision':head,'linux_owner_pid':os.getpid(),'linux_owner_start_ticks':birth(os.getpid()),'input_seals':seals,'arguments':args,'deadline_seconds':400,'errors':[],'outer_stop':False,'linux_child_pid':None,'linux_child_start_ticks':None,'scope':plan['purpose']}
def observed_absent(pid,ticks):
 if pid is None or ticks is None:return False
 try:return birth(pid)!=ticks
 except BaseException as e:record['errors'].append('birth observation: '+repr(e));return False
def reap_owned(c,ticks,timeout):
 try:return c.wait(timeout=timeout)
 except subprocess.TimeoutExpired:
  try:
   assert ticks and birth(c.pid)==ticks and os.getpgid(c.pid)==c.pid
   os.killpg(c.pid,signal.SIGKILL)
  except (FileNotFoundError,ProcessLookupError):pass
  finally:
   # Every post-launch exit path attempts to reap the known Popen handle.
   c.wait(timeout=10)
  return c.returncode
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
   c=subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True)
   record['cleanup_stop_observer']={'pid':c.pid,'start_ticks':ticks,'bootstrap':str(leaf_bootstrap.relative_to(ROOT)),'arguments':args}
   ticks=birth(c.pid);record['cleanup_stop_observer']['start_ticks']=ticks;assert ticks,'Linux observer birth unavailable'
   code=c.wait(timeout=20);assert code==0 and read(leaf_receipt)['pass']
 finally:
  if c is not None:
   try:code=reap_owned(c,ticks,5)
   except BaseException as e:record['errors'].append('observer retirement: '+repr(e));code=c.returncode
   record['cleanup_stop_observer'].update(returncode=code,linux_birth_absent=observed_absent(c.pid,ticks))
def recover_windows():
 c=None;ticks=None
 try:
  with (WORK/(action+'-emergency.stdout')).open('xb') as out,(WORK/(action+'-emergency.stderr')).open('xb') as err:
   c=subprocess.Popen(cleanup,stdout=out,stderr=err,start_new_session=True)
   record['emergency_observer']={'pid':c.pid,'start_ticks':ticks,'bootstrap':str(cleanup_bootstrap.relative_to(ROOT)),'arguments':cleanup}
   ticks=birth(c.pid);record['emergency_observer']['start_ticks']=ticks;assert ticks,'Linux observer birth unavailable'
   try:code=c.wait(timeout=60)
   except subprocess.TimeoutExpired:
    record['errors'].append('Emergency helper exceeded60s; tool ownership remains invalidated')
    assert cleanup_bootstrap.exists(),'Emergency helper exact Windows birth unavailable'
    stop_cleanup_observer()
    raise
  assert code==0
  receipt=read(WORK/'emergency.json');assert receipt['recorded_births_absent'] and not receipt['errors']
 finally:
  if c is not None:
   try:code=reap_owned(c,ticks,5)
   except BaseException as e:record['errors'].append('observer retirement: '+repr(e));code=c.returncode
   record['emergency_observer'].update(returncode=code,linux_birth_absent=observed_absent(c.pid,ticks))
child=None;code=None;attempted=False;start=time.monotonic()
try:
 assert not signals
 with prefix.with_suffix('.stdout').open('xb') as out,prefix.with_suffix('.stderr').open('xb') as err:
  attempted=True;child=subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True)
  record['linux_child_pid']=child.pid;record['linux_child_start_ticks']=birth(child.pid);assert record['linux_child_start_ticks'],'Linux native bridge birth unavailable'
  while child.poll() is None:
   if signals:raise RuntimeError(('controller_signals',signals))
   remaining=400-(time.monotonic()-start)
   if remaining<=0:raise subprocess.TimeoutExpired(args,400)
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
record['linux_child_birth_absent']=child is not None and observed_absent(child.pid,record.get('linux_child_start_ticks'))
record['accepted']=code==0 and not record['outer_stop'] and not record['errors'] and not signals and record['linux_child_birth_absent']

if record['accepted']:
 try:
  owner=read(output/'owner.json');bootstrap=read(temporary/'bootstrap.json')
  assert owner['status']=='completed' and owner['returncode']==0 and not owner['outer_stop'] and owner['numerical_control_passed'] and owner['owned_process_absent'] and not owner['cleanup_errors']
  assert owner['outer_guard_seconds']==320 and owner['rss_guard_bytes']==4294967296 and owner['peak_sampled_rss_bytes']<=4294967296
  assert owner['owner_pid']==bootstrap['owner_pid'] and owner['owner_start_time_utc']==bootstrap['owner_start_time_utc']
  assert owner['source_revision']==head==owner['live_source_revision'] and owner['source_build_gate_sha256']==plan['gate_sha256'] and owner['source_build_gate_unchanged']
  assert owner['test_executable_unchanged'] and owner['loaded_driver_modules_recorded']
  modules=read(output/'driver-modules.json');assert len(modules)==2
  for module in modules:
   expected=next(p for p in plan['providers'] if p['name'].lower()==module['name'].lower())
   assert module['path'].lower()==expected['path'].lower() and module['bytes']==expected['bytes'] and module['sha256']==expected['sha256']
  assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
  assert {seal_key(p):seal(p) for p in paths}==seals
  record['Windows_owner']=owner
  record['outputs']={seal_key(p):seal(p) for p in (output/'owner.json',output/'driver-modules.json',output/'stdout.log',output/'stderr.log',output/'gtest.xml',temporary/'bootstrap.json')}
 except BaseException as e:
  record['accepted']=False;record['errors'].append(repr(e))
  try:recover_windows()
  except BaseException as e:record['errors'].append(repr(e))

if signals:record['outer_stop']=True;record['accepted']=False
record['controller_signals']=list(signals)
record['linux_child_group_absent']=False
try:
 if record.get('linux_child_pid') is not None:
  os.killpg(record['linux_child_pid'],0)
except ProcessLookupError:record['linux_child_group_absent']=True
except BaseException as e:record['errors'].append('group observation: '+repr(e))
record['accepted']=record['accepted'] and record['linux_child_group_absent']
if signals:record['outer_stop']=True;record['accepted']=False
record['controller_signals']=list(signals)
prefix.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('action','returncode','elapsed_seconds','accepted') if k in record}),flush=True)
assert record['accepted']
