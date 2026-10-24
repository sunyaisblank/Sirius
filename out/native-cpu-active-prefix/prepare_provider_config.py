"""Reuse sealed original authority with the new task-local file context."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_config.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5ce6ffb5e7ff923cf0a34696ec088132915e90b57205b9da4cbaad7471c04762'
source = basis.read_text()
changes = {"'typed-three-stage-batch'": "'compile'"}
for before, after in changes.items():
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
