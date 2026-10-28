"""Original seven once in order after current AX ABBA/Decimal; stop on failure."""
from control_authority import ao, load

ao('run_controls.py')
exec(compile(load('run_controls.py'), 'out/native-cpu-local-identities-cohort/controls/run_controls.py', 'exec'))
