"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/read_linked_native.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'f6ad126b91333496b30f084f08bbdce415d75cd99cb46b4b9230548197dfb45e'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
