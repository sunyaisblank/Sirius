"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare_provider_build.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'e66f9a4ce9ca148e53c72d322c9e96e60c71b0de4b32541f3e0b6b6ffe473446'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
