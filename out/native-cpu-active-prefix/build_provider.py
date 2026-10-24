"""Reuse sealed original authority with the new task-local file context."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-slang-bridge/build_provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'c3d87afecb5ecb975c4dc49d8cf7a30f54e07cc553418d57e6f967b66b3b6325'
source = basis.read_text()
changes = {'native-cpu-slang-bridge/provider': 'native-cpu-active-prefix/provider'}
for before, after in changes.items():
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
