"""Root-owned finite Windows retirement census for the completed native control."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();HERE=Path(__file__).resolve().parent
assert len(sys.argv)==1
WORK=HERE;action='cleanup-retirement';signals=[]
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
def win(p):return subprocess.check_output(['wslpath','-w',str(p.resolve())],text=True,timeout=10).strip()
def quote(v):return "'"+v.replace("'","''")+"'"
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
cleanup_bootstrap=WORK/'cleanup-retirement-bootstrap.json'
record={'errors':[],'scope':'Fresh exact native control births + accessible scoped consumer census only.'}
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

args=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(HERE/'verify_final_retirement.ps1'),'-Stage',win(WORK),'-Receipt',win(WORK/'cleanup-native-retirement.json'),'-Bootstrap',win(cleanup_bootstrap),'-NativeStage',win(ROOT/read(HERE/'source.json')['stage_root'])]
input_paths=(Path(__file__),HERE/'verify_final_retirement.ps1',HERE/'source.json',*[WORK/case/path for case in ('device-timestamps','mixed-transport','factored-endpoints') for path in ('run/owner.json','owner-temp/bootstrap.json')])
def seal_inputs():return {str(p):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in input_paths}
record['arguments']=args;record['input_seals']=seal_inputs()
for s in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(s,lambda n,f:signals.append(n))
c=None;ticks=None;start=time.monotonic()
try:
 with (WORK/'cleanup-retirement.stdout').open('xb') as out,(WORK/'cleanup-retirement.stderr').open('xb') as err:
  c=subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True);record.update(pid=c.pid,start_ticks=None)
  ticks=birth(c.pid);record['start_ticks']=ticks;assert ticks
  while c.poll() is None:
   if signals or time.monotonic()-start>40:raise RuntimeError('Retirement observer exceeded40s or received signal')
   try:c.wait(timeout=.1)
   except subprocess.TimeoutExpired:pass
  assert c.returncode==0 and read(WORK/'cleanup-native-retirement.json')['pass']
except BaseException as e:
 record['errors'].append(repr(e))
 if c is not None and c.poll() is None:
  try:stop_cleanup_observer()
  except BaseException as e:record['errors'].append('Windows observer retirement: '+repr(e))
finally:
 if c is not None:
  try:record['returncode']=reap_owned(c,ticks,5)
  except BaseException as e:record['errors'].append('Linux observer retirement: '+repr(e))
  record['birth_absent']=observed_absent(c.pid,ticks)
  record['group_absent']=False
  try:os.killpg(c.pid,0)
  except ProcessLookupError:record['group_absent']=True
  except BaseException as e:record['errors'].append('group observation: '+repr(e))
try:
 record['inputs_unchanged']=seal_inputs()==record['input_seals']
 if not record['inputs_unchanged']:record['errors'].append('Retirement inputs changed during observer')
except BaseException as e:
 record['inputs_unchanged']=False;record['errors'].append('post-seal verification: '+repr(e))
record['elapsed_seconds']=time.monotonic()-start;record['signals']=signals
record['accepted']=not record['errors'] and not signals and record.get('returncode')==0 and record.get('birth_absent') and record.get('group_absent')
(WORK/'cleanup-retirement-run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record));assert record['accepted']
