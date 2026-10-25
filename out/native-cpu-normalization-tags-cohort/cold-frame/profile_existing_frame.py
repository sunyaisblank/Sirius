"""One existing3s USER-PC capture during the sole original frame, never a render."""
from pathlib import Path
import hashlib,json,os,resource,subprocess,time
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'frame';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
def proc(pid):
 v=(Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()
 return {'pid':pid,'state':v[0],'ppid':int(v[1]),'session':int(v[3]),'start_ticks':v[19]}
plan=J(W/'profile-plan.json');exe=R/plan['target_executable'];sampler=R/plan['sampler']
assert pin(sampler)==plan['sampler_identity'] and not (F/'all-threads-observation.json').exists()
record={'status':'started','plan':pin(W/'profile-plan.json'),'source':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'observer_started':False,'errors':[]}
p=W/'profile-controller.json'
with p.open('x') as f:json.dump(record,f)
def save():p.write_text(json.dumps(record,indent=2)+'\n')
began=time.monotonic()
try:
 while time.monotonic()-began<plan['wait_for_first_tile_seconds']:
  if not (F/'owner.json').is_file():time.sleep(.1);continue
  owner=J(F/'owner.json');assert owner['source']=={'revision':record['source'],'status':''}
  if owner['status']=='terminal':raise RuntimeError('original frame terminated before profile trigger')
  stdout=F/'stdout.log'
  if not stdout.exists() or b'Tile' not in stdout.read_bytes():time.sleep(.1);continue
  matches=[]
  for path in Path('/proc').iterdir():
   if not path.name.isdecimal():continue
   try:
    target=proc(int(path.name))
    if target['session']==owner['child_pid'] and target['state'] not in ('Z','X') and os.path.samefile(path/'exe',exe):matches.append(target)
   except (FileNotFoundError,PermissionError,ProcessLookupError):continue
  if len(matches)==1:break
  assert not matches,'more than one matching owned renderer';time.sleep(.1)
 else:raise TimeoutError('single predeclared first-tile trigger window expired; no retry')
 target=matches[0];assert proc(owner['child_pid'])['start_ticks']==owner['child_start_ticks']
 argv=['/usr/bin/python3','-B',str(sampler),str(target['pid']),target['start_ticks'],str(target['session']),str(F)]
 record.update(observer_started=True,target=target,argv=argv,trigger_delay_seconds=time.monotonic()-began);save()
 started=time.monotonic();code=subprocess.run(argv,cwd=R).returncode
 usage=resource.getrusage(resource.RUSAGE_CHILDREN);raw=J(F/'all-threads-observation.json')
 record.update(observer_exit=code,observer_wall_seconds=time.monotonic()-started,observer_peak_rss_kib=usage.ru_maxrss,observer_user_cpu_seconds=usage.ru_utime,observer_system_cpu_seconds=usage.ru_stime,raw_observation=pin(F/'all-threads-observation.json'),raw_status=raw['status'],sample_count=raw.get('sample_count',0),status='complete' if code==0 and raw['status']=='complete' else 'incomplete')
except BaseException as e:
 record.update(status='error');record['errors'].append(type(e).__name__+': '+str(e))
record['scope']='One fixed3s existing USER-PC observer during the original frame. Original command/settings/cache/governors/output assertions and300s4GiB2GiB limits stay unchanged. Renderer wall time includes observer interference; separate observer footprint reported, not subtracted. No uninstrumented speed, kernel/off-CPU/callchain or whole-frame-causality claim. No render or capture retry.'
save();print(json.dumps({k:record.get(k) for k in ('status','observer_exit','sample_count','observer_wall_seconds','observer_peak_rss_kib','errors')}));raise SystemExit(0 if record['status']=='complete' else 1)
