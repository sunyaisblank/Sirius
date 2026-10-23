"""Reuse a sealed original authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cooperative-transport-intergroup/readers/read_results.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '00a896fb01881d3fa3c90ddc0b1d0421f8bcf688a02e2c1752965cb74f41e118'
_source = _basis.read_text()
assert 'out/native-cooperative-transport-intergroup' in _source
_source = _source.replace('out/native-cooperative-transport-intergroup', 'out/lavapipe-owner-machine-codegen')
exec(compile(_source, str(_basis), "exec"))
