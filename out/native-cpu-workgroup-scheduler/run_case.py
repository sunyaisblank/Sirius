"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/run_case.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'a5b98807d45afbc56435e553367e8baed9d7ce6628b0d1c4735713f0494748a9'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
