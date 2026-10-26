"""Keep each original control invocation; move only timing prerequisites later."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from control_authority import ao, load

ao('runtime/run_case.py')
source = load('runtime/run_case.py')
before = """benefit = json.loads((current / 'benefit/readers/benefit.json').read_text())
decimal = json.loads((current / 'benefit/readers/decimal-agreement.json').read_text())
assert benefit['audit_passed'] and benefit['all22_numeric_pass']
assert decimal['agrees_on_all22'] and decimal['all22_pass']"""
after = """control8 = json.loads((current / 'control-8-acceptance.json').read_text())
assert control8['passed'] and control8['optional_math']['all_selected_optional_math_observed']
number = int(label.removeprefix('control-'))
assert action['prerequisites'] == ['control-8'] + ['control-' + str(n) for n in range(1, number)]"""
assert source.count(before) == 1
source = source.replace(before, after, 1)
exec(compile(source, 'out/native-cpu-local-identities-cohort/controls/runtime/run_case.py', 'exec'))
