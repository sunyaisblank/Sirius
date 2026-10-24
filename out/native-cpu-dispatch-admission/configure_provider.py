"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/configure_provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6eae5bac51f3059bf9d75d95d9dd4db77a56d7bedd8301efe4bd9cb128f2de07'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
