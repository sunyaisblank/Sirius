"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/runtime/owned_guard.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6fc32fd7fbed8e24370e929829a70b82aa53e54668fb38ded68ace2d1e9a6037'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
