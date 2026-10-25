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
              current / 'controls/runtime' / previous if previous in ('control-1', 'control-2', 'control-3') else
              work / previous)
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
for name, expected in action['expected_identities'].items():
    data = (root / name).read_bytes()
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected, name
original_argv = [action['executable'], '--gtest_filter=' + action['case'], '--gtest_color=no']
subprocess.run(original_argv, cwd=root, check=True)
