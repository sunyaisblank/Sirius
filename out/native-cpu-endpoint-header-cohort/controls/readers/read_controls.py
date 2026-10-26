from pathlib import Path
import sys
import json
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from control_authority import ao

# The original prepare-time timing gates now belong at final admission.
current = Path(__file__).resolve().parents[2]
benefit = current / 'benefit'
accepted = json.loads((benefit / 'readers/benefit.json').read_text())
decimal = json.loads((benefit / 'readers/decimal-agreement.json').read_text())
assert accepted['audit_passed'] and accepted['all22_numeric_pass'] and accepted['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass']
for label in ('a1', 'b1', 'b2', 'a2'):
    owner = json.loads((benefit / 'runtime' / label / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
exec(compile(ao('readers/read_controls.py'), '/home/astra/.project/Sirius/out/native-cpu-v4-cohort/controls/readers/read_controls.py', "exec"))
