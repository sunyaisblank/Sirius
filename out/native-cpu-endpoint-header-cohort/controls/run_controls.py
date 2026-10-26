"""Run the original seven once, fail fast; defer final acceptance until fresh ABBA."""
from control_authority import ao, load

ao('run_controls.py')
source = load('run_controls.py')
before = """accepted = json.loads((work.parent / 'benefit/readers/benefit.json').read_text())
decimal = json.loads((work.parent / 'benefit/readers/decimal-agreement.json').read_text())
assert accepted['audit_passed'] and accepted['all22_numeric_pass']
assert decimal['agrees_on_all22'] and decimal['all22_pass']"""
after = """accepted = json.loads((work.parent / 'control-8-acceptance.json').read_text())
assert accepted['passed'] and accepted['optional_math']['all_selected_optional_math_observed']"""
assert source.count(before) == 1
source = source.replace(before, after, 1)
before = "subprocess.run(['/usr/bin/python3', '-B', str(work / 'readers/read_controls.py')], check=True)"
assert source.count(before) == 1
source = source.replace(before,
    "print(json.dumps({'original_controls_completed': 7, 'final_all8_plus22_reader_pending': True}))", 1)
exec(compile(source, 'out/native-cpu-local-identities-cohort/controls/run_controls.py', 'exec'))
