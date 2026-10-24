"""Reuse sealed original authority with the new task-local file context."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-slang-bridge/build_batch.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'b1fec31561d7ae127be4a0ec46a4d036141d4cd80cf14409c36a9b18c921237f'
source = basis.read_text()
changes = {}
for before, after in changes.items():
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
