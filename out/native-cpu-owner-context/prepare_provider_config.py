"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare_provider_config.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '333cf237b0cf71129475b4c3e3c1c6eee9891c6e969e8f04891de062607ff81d'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
