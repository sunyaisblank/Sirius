"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/numeric_gates.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'f0c6b9fc00567c71cf738a0cc2b43797dc7d3b05aae3040ee349aaecc11962a7'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
