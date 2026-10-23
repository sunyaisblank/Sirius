"""Reuse a sealed original authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cooperative-transport-intergroup/readers/decimal_check.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '95661dbb3e82d8c02f15b9404865a07a6398ce252a9a5871abe4d0a3b63091a5'
_source = _basis.read_text()
assert 'out/native-cooperative-transport-intergroup' in _source
_source = _source.replace('out/native-cooperative-transport-intergroup', 'out/lavapipe-owner-bulk-zero')
exec(compile(_source, str(_basis), "exec"))
