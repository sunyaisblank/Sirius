"""Derive the finite root bridge from the previously reviewed birth-recovery owner."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
w=Path(__file__).resolve().parent
prior=(w/'prior-run-native-spm.py').read_text()
helpers=prior[prior.index('def reap_owned('):prior.index('child=None;code=None;attempted=False;')]
body=prior[prior.index('child=None;code=None;attempted=False;'):prior.index("if record['accepted']:")]
body=body.replace('remaining=300-', 'remaining=120-').replace('TimeoutExpired(args,300)','TimeoutExpired(args,120)')
prefix=r'''"""Root-owned read-only native feature admission; bounded bridge and scoped recovery."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def seal_key(p):return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p.resolve())
def win(p):return subprocess.check_output(['wslpath','-w',str(p.resolve())],text=True,timeout=10).strip()
def git(*a):return subprocess.check_output(['git',*a],text=True,timeout=20).strip()
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
plan=read(WORK/'plan.json');head=plan['source_revision'];action='clock-query'
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
review=read(WORK/'source-review.json');assert review['accepted']
paths=[Path(__file__),WORK/'plan.json',WORK/'source-review.json',WORK/'clock-owner.ps1',WORK/'clock-child.ps1',WORK/'native_clock_features.cs',WORK/'emergency_native_cleanup.ps1',WORK/'prepare_clock_owner.py',WORK/'prepare_clock_bridge.py',WORK/'owner-preparation.json']
paths += [Path(p) if Path(p).is_absolute() else ROOT/p for p in plan['input_seals']]
paths += [Path(p['linux_path']) for p in plan['providers']]
seals={seal_key(p):seal(p) for p in paths}
for name,v in review['controller_seals'].items():assert seal(WORK/name)==v
for name,v in plan['input_seals'].items():assert seals[name]==v
for p in plan['providers']:assert seal(Path(p['linux_path']))=={k:p[k] for k in ('bytes','sha256')}
output=WORK/'run';temporary=WORK/'owner-temp';prefix=WORK/'clock-query-bridge'
assert not output.exists() and not temporary.exists() and not (WORK/'native-temp').exists() and not prefix.with_suffix('.json').exists()
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
args=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'clock-owner.ps1'),'-Stage',win(WORK),'-Output',win(output),'-Filter','native-shader-clock-features','-Expected','native-shader-clock-features','-Seconds','45','-RssBytes','536870912','-TemporaryDirectory',win(temporary)]
def quote(v):return "'"+v.replace("'","''")+"'"
cleanup_bootstrap=WORK/(action+'-cleanup-helper-bootstrap.json')
cleanup_command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');scope='Exact emergency helper observer birth before scoped recovery'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(win(cleanup_bootstrap))+" -Encoding UTF8;& "+quote(win(WORK/'emergency_native_cleanup.ps1'))+" -TemporaryDirectory "+quote(win(temporary))+" -CompilerDirectory "+quote(win(WORK/'native-temp'))+" -Output "+quote(win(output))+" -Receipt "+quote(win(WORK/(action+'-emergency.json')))
cleanup=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',cleanup_command]
signals=[]
for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(signum,lambda signum,frame:signals.append(signum))
record={'action':action,'source_revision':head,'linux_owner_pid':os.getpid(),'linux_owner_start_ticks':birth(os.getpid()),'input_seals':seals,'arguments':args,'deadline_seconds':120,'errors':[],'outer_stop':False,'scope':plan['purpose']}
'''
suffix=r'''
if record['accepted']:
 try:
  owner=read(output/'owner.json');result=read(WORK/'clock-features.json');bootstrap=read(temporary/'bootstrap.json');child_birth=read(WORK/'child-bootstrap.json')
  assert owner['status']=='completed' and owner['returncode']==0 and not owner['outer_stop'] and owner['screen_execution_completed'] and owner['owned_process_absent'] and not owner['cleanup_errors']
  assert owner['outer_guard_seconds']==45 and owner['rss_guard_bytes']==536870912 and owner['peak_sampled_rss_bytes']<=536870912
  assert owner['owner_pid']==bootstrap['owner_pid'] and owner['owner_start_time_utc']==bootstrap['owner_start_time_utc']
  assert owner['test_pid']==child_birth['child_pid'] and owner['test_start_time_utc']==child_birth['child_start_time_utc']
  assert result['completed'] and result['device']=='AMD Radeon 780M Graphics' and result['driver']=='26.8.1 (LLPC)'
  assert result['logical_devices_created']==result['pipelines_created']==result['dispatches']==0
  assert result['extension_count']==len(result['extensions'])
  assert owner['test_executable_unchanged'] and owner['loaded_driver_modules_recorded']
  assert len(result['loaded_modules'])==2
  for module in result['loaded_modules']:
   expected=next(p for p in plan['providers'] if p['name'].lower()==module['name'].lower())
   assert module['path'].lower()==expected['path'].lower() and module['bytes']==expected['bytes'] and module['sha256']==expected['sha256']
  assert (output/'stderr.log').stat().st_size==0 and prefix.with_suffix('.stderr').stat().st_size==0
  assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
  assert {seal_key(p):seal(p) for p in paths}==seals
  record['Windows_owner']=owner;record['clock_features']=result['clock_features'];record['outputs']={seal_key(p):seal(p) for p in (output/'owner.json',output/'driver-modules.json',output/'stdout.log',output/'stderr.log',WORK/'clock-features.json',temporary/'bootstrap.json',WORK/'child-bootstrap.json')}
 except BaseException as e:
  record['accepted']=False;record['errors'].append(repr(e))
  try:recover_windows()
  except BaseException as e:record['errors'].append(repr(e))
if signals:record['outer_stop']=True;record['accepted']=False
record['controller_signals']=list(signals)
try:os.killpg(record['linux_child_pid'],0);record['linux_child_group_absent']=False
except ProcessLookupError:record['linux_child_group_absent']=True
record['accepted']=record['accepted'] and record['linux_child_group_absent']
prefix.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('action','returncode','elapsed_seconds','accepted','clock_features') if k in record}),flush=True)
assert record['accepted']
'''
p=w/'run_clock_query.py';assert not p.exists();p.write_text(prefix+helpers+body+suffix)
import ast
ast.parse(p.read_text())
print(p.stat().st_size)
