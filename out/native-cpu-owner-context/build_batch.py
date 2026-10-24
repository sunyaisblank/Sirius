"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/build_batch.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6b92bd64e5714d73dabb873f9f413db7a79a2f144117a8a3043d105ce3a0bf6d'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
