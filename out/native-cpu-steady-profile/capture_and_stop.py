"""Run original frame command, sample once in its body, then stop owned child."""
from pathlib import Path
import hashlib,json,os,re,resource,signal,subprocess,time
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'observe'
J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
def proc(pid):
 v=(Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()
 return dict(pid=pid,state=v[0],ppid=int(v[1]),session=int(v[3]),start_ticks=v[19])
def stop(saved):
 try:
  fd=os.pidfd_open(saved['pid'])
  try:
   now=proc(saved['pid']);assert now['start_ticks']==saved['start_ticks'] and now['session']==saved['session']
   signal.pidfd_send_signal(fd,signal.SIGTERM)
  finally:os.close(fd)
 except (FileNotFoundError,ProcessLookupError):pass
plan=J(W/'plan.json');exe=R/plan['executable'];sampler=R/plan['sampler']
assert pin(exe)==plan['executable_identity'] and pin(sampler)==plan['sampler_identity']
assert not (F/'all-threads-observation.json').exists()
result=dict(status='started',source=plan['source'],observer_started=False,frame_qualification=False,errors=[])
output=W/'result.json'
with output.open('x') as f:json.dump(result,f)
child=None;target=None;child_identity=None
try:
 child=subprocess.Popen(plan['original_argv'],cwd=R)
 child_identity=proc(child.pid);began=time.monotonic()
 while time.monotonic()-began<plan['wait_seconds']:
  assert child.poll() is None,'original frame exited before diagnostic trigger'
  if not (F/'owner.json').exists():time.sleep(.05);continue
  owner=J(F/'owner.json')
  assert owner['source']==dict(revision=plan['source'],status='') and owner['status']!='terminal'
  progress=re.findall(rb'\bTile\s+(\d+)\s*/\s*(\d+)',(F/'stdout.log').read_bytes())
  eligible=[(int(d),int(t)) for d,t in progress if int(t)==plan['trigger_total_tiles'] and plan['trigger_completed_tiles']<=int(d)<=int(t)]
  if not eligible:time.sleep(.05);continue
  matches=[]
  for p in Path('/proc').iterdir():
   if not p.name.isdecimal():continue
   try:
    v=proc(int(p.name))
    if v['session']==owner['child_pid'] and v['state'] not in ('Z','X') and os.path.samefile(p/'exe',exe):matches.append(v)
   except (FileNotFoundError,ProcessLookupError,PermissionError):continue
  assert len(matches)==1,'expected exactly one owned renderer'
  target=matches[0];break
 else:raise TimeoutError('single30s128-completed-tile trigger expired; no retry')
 result.update(target=target,trigger_progress=eligible[-1],trigger_delay_seconds=time.monotonic()-began,observer_started=True)
 output.write_text(json.dumps(result,indent=2)+'\n')
 argv=['/usr/bin/python3','-B',str(sampler),str(target['pid']),target['start_ticks'],str(target['session']),str(F)]
 started=time.monotonic();code=subprocess.run(argv,cwd=R).returncode
 usage=resource.getrusage(resource.RUSAGE_CHILDREN);raw=J(F/'all-threads-observation.json')
 result.update(observer_argv=argv,observer_exit=code,observer_wall_seconds=time.monotonic()-started,
  observer_peak_rss_kib=usage.ru_maxrss,observer_user_cpu_seconds=usage.ru_utime,observer_system_cpu_seconds=usage.ru_stime,
  raw_observation=pin(F/'all-threads-observation.json'),sample_count=raw.get('sample_count',0),raw_status=raw['status'])
 assert code==0 and raw['status']=='complete','single capture incomplete; no retry'
 result['status']='complete'
except BaseException as e:
 result['status']='error';result['errors'].append(type(e).__name__+': '+str(e))
finally:
 if child is not None:
  try:
   if target is not None:stop(target)
   try:child.wait(timeout=5)
   except subprocess.TimeoutExpired:stop(child_identity);child.wait(timeout=5)
   result['original_command_exit']=child.returncode
   result['intentional_stop']=True
  except BaseException as e:
   result['status']='error';result['errors'].append('cleanup '+type(e).__name__+': '+str(e))
result['scope']=plan['scope'];output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result.get(k) for k in ('status','trigger_progress','sample_count','observer_wall_seconds','original_command_exit','errors')}),flush=True)
raise SystemExit(0 if result['status']=='complete' else 1)
