"""Reuse sealed original authority with the new task-local file context."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-slang-bridge/read_linked_native.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '9a94c2b45471807ad74b28060edc60e19f33d42e7b9f94ba6126065d79c6f064'
source = basis.read_text()
changes = {"batch = json.loads((W / 'batch-generation.json').read_text())": "batch = json.loads((R / 'out/native-cpu-slang-bridge/batch-generation.json').read_text())"}
for before, after in changes.items():
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
