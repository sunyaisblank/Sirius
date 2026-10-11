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
owner=birth(os.getpid());probe=W/'decode_with_saved_isa.py';inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),probe,W/'plan.json')}
pending_signals=[]
def defer_signal(sig,frame):pending_signals.append(sig)
signal.signal(signal.SIGTERM,defer_signal);signal.signal(signal.SIGINT,defer_signal)
plan=json.loads((W/'plan.json').read_text())
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==plan['source_revision']
assert not subprocess.check_output(['git','status','--porcelain'],text=True).strip()
assert not (W/'isa-run.json').exists()
start=time.monotonic();peak=0;reason=None;child=None;known={};returncode=None;errors=[];remaining=None
def signal_group(sig):
 try:os.killpg(child.pid,sig)
 except ProcessLookupError:pass
 except BaseException as e:errors.append({'action':'signal_group','signal':sig,'error':repr(e)})
with (W/'isa-stdout.txt').open('xb') as stdout,(W/'isa-stderr.txt').open('xb') as stderr:
 try:
  if pending_signals:reason='owner_signal';raise InterruptedError('Interrupted before spawn')
  child=subprocess.Popen([sys.executable,'-B',str(probe)],cwd=ROOT,stdin=subprocess.PIPE,stdout=stdout,stderr=stderr,start_new_session=True)
  b=birth(child.pid);assert b['parent_pid']==owner['pid'] and b['process_group']==child.pid
  known[(b['pid'],b['start_ticks'])]=b
  (W/'isa-owner.json').write_text(json.dumps({'owner':owner,'child':b,'input_seals':inputs,'wall_seconds':45,'rss_bytes':512*1024**2},indent=2)+'\n')
  if pending_signals:reason='owner_signal';raise InterruptedError('Interrupted before START')
  child.stdin.write(b'START\n');child.stdin.flush();child.stdin.close()
  while child.poll() is None:
   if pending_signals:reason='owner_signal';break
   for p in census(child.pid):known[(p['pid'],p['start_ticks'])]=p
   status=Path('/proc/'+str(child.pid)+'/status')
   try:
    rss=next((int(line.split()[1])*1024 for line in status.read_text().splitlines() if line.startswith('VmRSS:')),0);peak=max(peak,rss)
   except FileNotFoundError:rss=0
   if time.monotonic()-start>45 or rss>512*1024**2:
    reason='wall_limit' if time.monotonic()-start>45 else 'rss_limit';break
   time.sleep(.05)
 except BaseException as e:errors.append({'action':'run','error':repr(e)})
 finally:
  if child is not None:
   try:
    if child.poll() is None:signal_group(signal.SIGTERM)
    try:returncode=child.wait(timeout=2)
    except subprocess.TimeoutExpired:
     signal_group(signal.SIGKILL)
     try:returncode=child.wait(timeout=3)
     except BaseException as e:errors.append({'action':'leader_wait_after_kill','error':repr(e)})
   except BaseException as e:errors.append({'action':'leader_wait','error':repr(e)})
   try:remaining=census(child.pid)
   except BaseException as e:errors.append({'action':'group_census','error':repr(e)})
   if remaining:
    signal_group(signal.SIGKILL);time.sleep(.1)
   try:returncode=child.wait(timeout=1)
   except BaseException as e:errors.append({'action':'final_leader_wait','error':repr(e)})
   try:remaining=census(child.pid)
   except BaseException as e:remaining=None;errors.append({'action':'final_group_census','error':repr(e)})
closures=[]
for b in known.values():
 try:current=birth(b['pid']);absent=current['start_ticks']!=b['start_ticks']
 except (FileNotFoundError,ProcessLookupError):current=None;absent=True
 except BaseException as e:current=None;absent=False;errors.append({'action':'birth_closure','birth':b,'error':repr(e)})
 closures.append({'birth':b,'current':current,'recorded_birth_absent':absent})
for p,v in inputs.items():assert seal(ROOT/p)==v
if pending_signals:reason='owner_signal'
result={'owner':owner,'input_seals':inputs,'recorded_children':list(known.values()),'recorded_child_closure':closures,'final_owned_group':remaining,'returncode':returncode,'limit_reason':reason,'deferred_signals':pending_signals,'errors':errors,'peak_sampled_rss_bytes':peak,'elapsed_seconds':time.monotonic()-start,'outputs':{p.name:seal(p) for p in W.glob('isa-*') if p.is_file()},'scope':'Recorded Linux child births and final owned-group closure only; no global process, hidden task or scheduled-job absence claim.'}
(W/'isa-run.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'returncode':returncode,'limit_reason':reason,'recorded_children':len(known),'closure':all(v['recorded_birth_absent'] for v in closures),'peak_sampled_rss_bytes':peak,'elapsed_seconds':result['elapsed_seconds']}))
assert returncode==0 and reason is None and not pending_signals and not errors and remaining==[] and all(v['recorded_birth_absent'] for v in closures)
