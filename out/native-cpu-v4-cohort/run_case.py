"""Run the exact original C8 authority; separately read bounded optional coverage."""
from pathlib import Path
from rebind import load

source = load('run_case.py')
exec(compile(source, str(Path.cwd() / 'out/native-cpu-local-identities-cohort/run_case.py'), 'exec'))
from optional_math import check
optional = check(root, action, record)
print(json.dumps({'optional_math': optional}))
