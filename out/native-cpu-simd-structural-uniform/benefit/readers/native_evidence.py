"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/native_evidence.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5bde8a6e69d6e7f7ed0d581787aef33b4696be9c9d9fc718c31e067eefe82097'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
