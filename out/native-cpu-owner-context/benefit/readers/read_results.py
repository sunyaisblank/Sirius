"""Reuse sealed B authority with task-local context and unchanged gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/read_results.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '1371748434db9de22af3555de6e2b27b25e9bd25f7928d1937f8f0202a7ace92'
source = basis.read_text()
changes = [('out/native-cpu-active-prefix', 'out/native-cpu-owner-context')]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
