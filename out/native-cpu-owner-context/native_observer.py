"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/native_observer.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'c38d87c7cea09dd50b53750982987ab5bfdd46859b53d0d61598a9e7c66377fa'
source = basis.read_text()
changes = [('out/native-cpu-active-prefix', 'out/native-cpu-owner-context')]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
