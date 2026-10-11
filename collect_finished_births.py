"""Collect recorded process births from finished diagnostic receipts only."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
WORK = Path(__file__).resolve().parent
windows = {}
linux = {}
def collect(node, origin):
    if isinstance(node, dict):
        for prefix in ('owner', 'test', 'helper', 'observer'):
            pid, date = node.get(prefix + '_pid'), node.get(prefix + '_start_time_utc')
            if type(pid) is int and pid > 0 and isinstance(date, str):
                windows.setdefault((pid, date), set()).add(origin)
        pid, date = node.get('pid'), node.get('birth')
        if type(pid) is int and pid > 0 and isinstance(date, str) and 'T' in date:
            windows.setdefault((pid, date), set()).add(origin)
        for prefix in ('linux_owner', 'linux_child'):
            pid = node.get(prefix + '_pid')
            ticks = node.get(prefix + '_start_ticks')
            if type(pid) is int and pid > 0 and isinstance(ticks, str) and ticks.isdecimal():
                linux.setdefault((pid, ticks), set()).add(origin)
        pid, ticks = node.get('pid'), node.get('start_ticks')
        if type(pid) is int and pid > 0 and isinstance(ticks, str) and ticks.isdecimal():
            linux.setdefault((pid, ticks), set()).add(origin)
        for child in node.values():
            collect(child, origin)
    elif isinstance(node, list):
        for child in node:
            collect(child, origin)
for path in sorted(WORK.rglob('*.json')):
    if path.name == 'finished-births.json':
        continue
    collect(json.loads(path.read_text(encoding='utf-8-sig')), str(path.relative_to(WORK)))
assert windows and linux
record = {'windows': [{'pid': p, 'birth': b, 'origins': sorted(v)}
                      for (p, b), v in sorted(windows.items())],
          'linux': [{'pid': p, 'start_ticks': b, 'origins': sorted(v)}
                    for (p, b), v in sorted(linux.items())],
          'scope': 'Recorded receipts only; no inference about unrecorded or hidden consumers.'}
destination = WORK / 'finished-births.json'
assert not destination.exists()
destination.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'windows_births': len(windows), 'linux_births': len(linux)}))
