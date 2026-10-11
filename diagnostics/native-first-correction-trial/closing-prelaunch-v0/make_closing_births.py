"""Join recorded scientific and audit births; no launches, stops or process census."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
assert WORK == ROOT/'out/native-division-first-correction-screen'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
review = read(WORK/'closing-retirement-source-review.json')
assert review['accepted'] and review['inventory_reader_seal'] == seal(Path(__file__))
inputs = {}; windows = []; linux = []
def get(p):
    inputs[str(p.relative_to(ROOT))] = seal(p)
    return read(p)
def native(folder, label):
    owner = get(folder/'run/owner.json'); bootstrap = get(folder/'owner-temp/bootstrap.json')
    assert owner['owner_pid'] == bootstrap['owner_pid'] and owner['owner_start_time_utc'] == bootstrap['owner_start_time_utc']
    assert owner['status'] == 'completed' and owner['returncode'] == 0 and not owner['outer_stop']
    for kind, pid, birth in [('owner',owner['owner_pid'],owner['owner_start_time_utc']),('child',owner['test_pid'],owner['test_start_time_utc'])]:
        windows.append({'category':'scientific','kind':label+' '+kind,'pid':int(pid),'birth':birth})
    audit = get(folder/'retirement-bootstrap.json'); retired = get(folder/'native-retirement.json')
    assert retired['pass'] and not retired['matches'] and not retired['errors']
    assert retired['observer_pid'] == audit['owner_pid'] and retired['observer_birth'] == audit['owner_start_time_utc']
    windows.append({'category':'audit','kind':label+' retirement observer','pid':int(audit['owner_pid']),'birth':audit['owner_start_time_utc']})
    bridge = get(folder/('placement-query-bridge.json' if label == 'private-query' else 'native-control-bridge.json'))
    rr = get(folder/'retirement-run.json')
    assert bridge['accepted'] and rr['accepted'] and rr['inputs_unchanged'] and not rr['errors']
    for kind, pid, ticks in [('owner',bridge['linux_owner_pid'],bridge['linux_owner_start_ticks']),('bridge',bridge['linux_child_pid'],bridge['linux_child_start_ticks']),('retirement',rr['pid'],rr['start_ticks'])]:
        linux.append({'kind':label+' '+kind,'pid':int(pid),'start_ticks':str(ticks)})
native(WORK/'placement', 'private-query')
actions = get(WORK/'runtime/native-cases.json')['actions']; assert len(actions) == 16
for action in actions: native(WORK/'runtime'/action, action)
for name in ('compiler.json','runtime/source-build.json','runtime/product-inspection.json'):
    report = get(WORK/name); assert report['accepted'] and not report['errors'] and not report['signals']
    for run in report['runs']:
        assert run['accepted'] and run['returncode'] == 0 and run['owned_birth_absent'] and run['owned_group_absent']
        linux.append({'kind':name+' '+run['name'],'pid':int(run['pid']),'start_ticks':str(run['start_ticks'])})
assert len(windows) == 51 and len({(t['pid'],t['birth']) for t in windows}) == 51
assert len(linux) == 62 and len({(t['pid'],t['start_ticks']) for t in linux}) == 62
inputs[str(Path(__file__).relative_to(ROOT))] = seal(Path(__file__))
inputs[str((WORK/'closing-retirement-source-review.json').relative_to(ROOT))] = seal(WORK/'closing-retirement-source-review.json')
for path, expected in inputs.items(): assert seal(ROOT/path) == expected, path
out = WORK/'closing-births.json'; assert not out.exists()
out.write_text(json.dumps({'accepted': True,'windows_targets':windows,'linux_targets':linux,'input_seals':inputs,
    'scope':'Exact recorded34 scientific and17 audit Windows births and62 Linux births. No current absence, group census, source-epoch revalidation, numerical or full qualification claim.'}, indent=2)+'\n')
print(json.dumps({'windows_scientific':34,'windows_audit':17,'linux_births':62,'inputs':len(inputs),'output':seal(out)}))
