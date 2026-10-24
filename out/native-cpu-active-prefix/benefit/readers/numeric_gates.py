"""Reuse the sealed original authority; path binding changes only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-arithmetic-sharing/readers/numeric_gates.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '051b1d0a83c5669099ca00b0496302d662d001b0f341d291a7cf89f68e3abdf1'
source = basis.read_text()
changes = []
for before, after in changes:
    assert source.count(before) >= 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
