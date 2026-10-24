"""Reuse the sealed original authority; path binding changes only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '1b44f0678dbfc7ee455b3d49e170be884b294882896fd1ac3706ff5ebded20ab'
source = basis.read_text()
changes = []
for before, after in changes:
    assert source.count(before) >= 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
