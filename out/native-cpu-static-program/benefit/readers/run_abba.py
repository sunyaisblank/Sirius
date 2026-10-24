"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/run_abba.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'c363a719e652fcb289bded3cf70124f2f1e0fc30960f42bcdb9301d1fe1cc3df'
source = basis.read_text()
changes = [('out/native-cpu-active-prefix', 'out/native-cpu-static-program')]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
