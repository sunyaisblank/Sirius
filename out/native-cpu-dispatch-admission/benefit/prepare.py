"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'ba3c35750f9e1b0b8e876b084d1f3f738ab34d1080709080dd5d4dc76915a440'
source = basis.read_text()
changes = [('current B', 'current F'), ('native B', 'native F')]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
