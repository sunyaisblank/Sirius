from pathlib import Path
import json
root=Path('out/native-shader-clock-execution')
source=Path('attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/run_placement_query.py').read_text()
start=source.index('ROOT=Path.cwd();'); end=source.index('def observed_absent(')
setup='''ROOT=Path.cwd();HERE=Path(__file__).resolve().parent
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
'''
source=source[:start]+setup+source[end:]
source=source.replace("WORK/(action+'-emergency.json')", "WORK/'emergency.json'")
source=source.replace('remaining=120-(time.monotonic()-start)', 'remaining=400-(time.monotonic()-start)').replace('subprocess.TimeoutExpired(args,120)', 'subprocess.TimeoutExpired(args,400)')
a=source.index("if record['accepted']:\n try:");b=source.index('\nif signals:record',a)
validation='''if record['accepted']:
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
'''
source=source[:a]+validation+source[b:]
source=source.replace("('action','returncode','elapsed_seconds','accepted','clock_placement')", "('action','returncode','elapsed_seconds','accepted')")
(root/'run_native_control.py').write_text(source)
# Reuse the accepted finite retirement observer for a named case. Root owns all invocations.
p=root/'run_native_retirement.py';s=p.read_text();s=s.replace("ROOT=Path.cwd();WORK=Path(__file__).resolve().parent;action='retirement';signals=[]", "ROOT=Path.cwd();HERE=Path(__file__).resolve().parent\nassert len(sys.argv)==2 and sys.argv[1] in ('device-timestamps','mixed-transport','factored-endpoints')\nWORK=HERE/sys.argv[1];action='retirement';signals=[]")
s=s.replace("WORK/'verify_native_retirement.ps1'", "HERE/'verify_native_retirement.ps1'")
s=s.replace('completed query', 'completed native control').replace('native query births', 'native control births')
p.write_text(s)
p=root/'verify_native_retirement.ps1';s=p.read_text().replace('finite-native-clock-query-retirement','finite-native-clock-control-retirement').replace('native query birth','native control birth');p.write_text(s)
print('prepared controller sources only')
