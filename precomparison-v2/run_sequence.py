"""Launch only the four frozen serial owner actions, once."""
from pathlib import Path
import hashlib,json,signal,subprocess,sys
sys.dont_write_bytecode=True
W=Path(__file__).resolve().parent; R=W.parents[1]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
freeze=json.loads((W/'retention-gate-frozen.json').read_text())
assert freeze['sequence']==['a1','b1','b2','a2']
assert freeze['launcher_sha256']==sha(__file__)
assert freeze['reader_sha256']==sha(W/'software_contracts.py')
assert freeze['comparator_sha256']==sha(W/'compare_software.py')
assert freeze['actions_sha256']==sha(W/'actions.json')
assert not any((W/label).exists() for label in freeze['sequence'])
path=W/'software-sequence.json'; assert not path.exists()
record={'labels':freeze['sequence'],'freeze_sha256':sha(W/'retention-gate-frozen.json'),
    'launcher_sha256':sha(__file__),'results':[],'completed':False}
def write():
    temporary=path.with_suffix('.tmp'); temporary.write_text(json.dumps(record,indent=2)+'\n'); temporary.replace(path)
child=None; stop=None
def interrupted(number,_):
    global stop
    stop=number
    # This direct child is unreaped, so its PID cannot yet be reused. The
    # reviewed owner handles the signal and contains its actual device session.
    if child is not None and child.returncode is None: child.send_signal(signal.SIGTERM)
for number in [signal.SIGINT,signal.SIGTERM,signal.SIGHUP]:signal.signal(number,interrupted)
write()
for label in record['labels']:
    assert stop is None
    child=subprocess.Popen([sys.executable,'-B',str(W/'run_owned.py'),label],cwd=R)
    code=child.wait()
    assert code==0 and stop is None,(label,code,stop)
    record['results'].append({'label':label,'owner_sha256':sha(W/label/'owner.json'),'exit':code}); write()
record['completed']=True; write()
print(json.dumps({'completed':True,'sequence':record['labels'],'results':len(record['results'])}))
