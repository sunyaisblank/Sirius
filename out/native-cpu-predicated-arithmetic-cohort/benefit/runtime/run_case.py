"""Reuse original ABBA cases, counts, limits and numerical/Decimal gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/runtime/run_case.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '25acbe181fe3df9238dcf51e211c789672ec761d05f0cae8275e7b14deefb7fb'
source = basis.read_text()
assert 'out/native-cpu-active-prefix' in source
source = source.replace('out/native-cpu-active-prefix', 'out/native-cpu-predicated-arithmetic-cohort')
exec(compile(source, str(basis), 'exec'))
