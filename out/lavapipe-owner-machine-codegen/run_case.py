"""Reuse a sealed original authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cooperative-transport-intergroup/run_case.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '32e6e275147046fd27be75b2f022b01431d36bc909357e519a9de32ce23d32fe'
_source = _basis.read_text()
exec(compile(_source, str(_basis), "exec"))
