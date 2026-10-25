"""Reuse sealed original acceptance; current artifact-path binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '1b44f0678dbfc7ee455b3d49e170be884b294882896fd1ac3706ff5ebded20ab'
source = basis.read_text()
changes = []
for old, new in changes:
    assert old in source, old
    source = source.replace(old, new)
exec(compile(source, str(basis), 'exec'))
