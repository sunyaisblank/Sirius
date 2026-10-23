"""Reuse a sealed original authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cooperative-transport-intergroup/readers/run_abba.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'fc334cb14225acda8afbf4fdf74f5bf1c75cdcedb3d48778f8c00c0090c3978d'
_source = _basis.read_text()
assert 'out/native-cooperative-transport-intergroup' in _source
_source = _source.replace('out/native-cooperative-transport-intergroup', 'out/native-cpu-arithmetic-sharing')
exec(compile(_source, str(_basis), "exec"))
