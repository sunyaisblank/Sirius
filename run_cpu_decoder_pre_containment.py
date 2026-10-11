"""Own, limit, reap and audit the single saved-byte native decoder child."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def birth(pid):
 s=Path('/proc/'+str(pid)+'/stat').read_text();f=s[s.rfind(')')+2:].split()
 return {'pid':pid,'start_ticks':int(f[19]),'process_group':int(f[2]),'parent_pid':int(f[1])}
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def census(group):
 found=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   b=birth(int(p.name))
   if b['process_group']==group:found.append(b)
  except (FileNotFoundError,ProcessLookupError):pass
 return found
owner=birth(os.getpid());probe=W/'decode_saved_sqtt.py';inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),probe,W/'plan.json')}
plan=json.loads((W/'plan.json').read_text())
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==plan['source_revision']
assert not subprocess.check_output(['git','status','--porcelain'],text=True).strip()
assert not (W/'decoder-run.json').exists()
start=time.monotonic();peak=0;reason=None;child=None;known={};returncode=None
with (W/'decoder-stdout.txt').open('xb') as stdout,(W/'decoder-stderr.txt').open('xb') as stderr:
 try:
  child=subprocess.Popen([sys.executable,'-B',str(probe)],cwd=ROOT,stdin=subprocess.PIPE,stdout=stdout,stderr=stderr,start_new_session=True)
  b=birth(child.pid);assert b['parent_pid']==owner['pid'] and b['process_group']==child.pid
  known[(b['pid'],b['start_ticks'])]=b
  (W/'decoder-owner.json').write_text(json.dumps({'owner':owner,'child':b,'input_seals':inputs,'wall_seconds':45,'rss_bytes':512*1024**2},indent=2)+'\n')
  child.stdin.write(b'START\n');child.stdin.flush();child.stdin.close()
  while child.poll() is None:
   for p in census(child.pid):known[(p['pid'],p['start_ticks'])]=p
   status=Path('/proc/'+str(child.pid)+'/status')
   try:
    rss=next((int(line.split()[1])*1024 for line in status.read_text().splitlines() if line.startswith('VmRSS:')),0);peak=max(peak,rss)
   except FileNotFoundError:rss=0
   if time.monotonic()-start>45 or rss>512*1024**2:
    reason='wall_limit' if time.monotonic()-start>45 else 'rss_limit';break
   time.sleep(.05)
 finally:
  if child is not None:
   if child.poll() is None:
    os.killpg(child.pid,signal.SIGTERM)
    try:child.wait(timeout=2)
    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL)
   returncode=child.wait(timeout=3)
   remaining=census(child.pid)
   if remaining:
    os.killpg(child.pid,signal.SIGKILL);time.sleep(.1);remaining=census(child.pid)
   assert not remaining
closures=[]
for b in known.values():
 try:current=birth(b['pid']);absent=current['start_ticks']!=b['start_ticks']
 except (FileNotFoundError,ProcessLookupError):current=None;absent=True
 closures.append({'birth':b,'current':current,'recorded_birth_absent':absent});assert absent
for p,v in inputs.items():assert seal(ROOT/p)==v
result={'owner':owner,'input_seals':inputs,'recorded_children':list(known.values()),'recorded_child_closure':closures,'final_owned_group':remaining,'returncode':returncode,'limit_reason':reason,'peak_sampled_rss_bytes':peak,'elapsed_seconds':time.monotonic()-start,'outputs':{p.name:seal(p) for p in W.glob('decoder-*') if p.is_file()},'scope':'Recorded Linux child births and final owned-group closure only; no global process, hidden task or scheduled-job absence claim.'}
(W/'decoder-run.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'returncode':returncode,'limit_reason':reason,'recorded_children':len(known),'closure':all(v['recorded_birth_absent'] for v in closures),'peak_sampled_rss_bytes':peak,'elapsed_seconds':result['elapsed_seconds']}))
assert returncode==0 and reason is None
