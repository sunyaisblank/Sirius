"""One uninstrumented original ABBA slot on the exact same matched host."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
label, output = sys.argv[1:]
assert label in ('a1', 'b1', 'b2', 'a2')
folder = Path(output)
assert folder == work / label
action = json.loads((work / 'actions.json').read_text())[label]
control = root / 'out/native-cpu-active-prefix/control-8-acceptance.json'
assert json.loads(control.read_text())['passed']
for previous in action['prerequisites']:
    owner_path = (root / 'out/native-cpu-active-prefix/control-8/owner.json'
                  if previous == 'control-8' else work / previous / 'owner.json')
    owner = json.loads(owner_path.read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
for name, expected in action['expected_identities'].items():
    path = root / name
    data = path.read_bytes()
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected, name
subprocess.run([action['executable'], '--gtest_filter=' + action['case'], '--gtest_color=no'],
               cwd=root, check=True)
