"""Execute one unchanged original control, once, after its original cohort."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
current = work.parents[1]
label, output = sys.argv[1:]
assert label in ['control-' + str(i) for i in range(1, 8)]
assert Path(output) == work / label
action = json.loads((work / 'actions.json').read_text())[label]
benefit = json.loads((current / 'benefit/readers/benefit.json').read_text())
decimal = json.loads((current / 'benefit/readers/decimal-agreement.json').read_text())
assert benefit['audit_passed'] and benefit['all22_numeric_pass']
assert decimal['agrees_on_all22'] and decimal['all22_pass']
for previous in action['prerequisites']:
    folder = (current / 'control-8' if previous == 'control-8' else
              current / 'benefit/runtime' / previous if previous in ('a1', 'b1', 'b2', 'a2') else
              work / previous)
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
for name, expected in action['expected_identities'].items():
    data = (root / name).read_bytes()
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected, name
original_argv = [action['executable'], '--gtest_filter=' + action['case'], '--gtest_color=no']
if label == 'control-4':
    commands = work / label / 'observe.gdb'
    commands.write_text('set pagination off\nset confirm off\nset debuginfod enabled off\n'
        'set auto-load off\nset startup-with-shell off\nset disable-randomization off\n'
        'set print thread-events off\nsource ' + str(work / 'control4_observer.py') +
        '\nrun\npython resume_after_sample()\n')
    subprocess.run(['/usr/bin/gdb', '-nx', '-nh', '-q', '-batch', '-x', str(commands),
                    '--args', *original_argv], cwd=root, check=True)
    observation = json.loads((work / label / 'route-observation.json').read_text())
    assert not observation['failures'] and observation['inferior_exit'] == 0, observation
else:
    subprocess.run(original_argv, cwd=root, check=True)
