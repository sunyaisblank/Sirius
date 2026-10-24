"""One original C8 under passive observation; no repeat, build, or warmup."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
label, output = sys.argv[1:]
assert label == 'control-8'
folder = Path(output)
assert folder == work / label
action = json.loads((work / 'actions.json').read_text())[label]
for name, expected in action['expected_identities'].items():
    path = Path(name) if name.startswith('/') else root / name
    data = path.read_bytes()
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected, name
commands = folder / 'observe.gdb'
commands.write_text('set pagination off\nset confirm off\nset debuginfod enabled off\n'
    'set auto-load off\nset startup-with-shell off\nset disable-randomization off\n'
    'source ' + str(work / 'native_observer.py') + '\nrun\n')
argv = ['/usr/bin/gdb', '-nx', '-nh', '-q', '-batch', '-x', str(commands), '--args',
        action['executable'], '--gtest_filter=' + action['case'], '--gtest_color=no']
result = subprocess.run(argv, cwd=root, env=os.environ.copy())
assert result.returncode == 0, result.returncode
record = json.loads((folder / 'native-observation.json').read_text())
assert not record['failures'] and record['inferior_exit'] == 0, record
assert len(record['calls']) == 3
assert {row['function'] for row in record['calls']} == {
    'sirius_native_' + slug + '_try_packet' for slug in ('transport', 'endpoint', 'dopri_phase')}
assert all(row['returned'] and row['handled'] for row in record['calls'])
assert (folder / 'gtest.xml').is_file()
print(json.dumps({'original_C8_child_exit': 0, 'first_native_packets_handled': 3,
                  'timing_is_not_a_benefit_trial': True}))
