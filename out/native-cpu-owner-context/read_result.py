"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/read_result.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '05926a8e76d0786ac344b811c6b170378c79ec7104beb1981549d8c69639e87d'
source = basis.read_text()
changes = [('providerB', 'providerC')]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
