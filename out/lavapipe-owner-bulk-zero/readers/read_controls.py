"""Reuse a sealed original authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cooperative-transport-intergroup/readers/read_controls.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '1c8ce9f9870a2e36bc885e035d108b2cd20180eb0d0b1be1244777f359c9e33e'
_source = _basis.read_text()
assert 'out/native-cooperative-transport-intergroup/compiler-result.json' in _source
_source = _source.replace('out/native-cooperative-transport-intergroup/compiler-result.json', 'out/lavapipe-owner-bulk-zero/compiler-result.json')
assert 'out/lavapipe-intergroup-row-owner/ordinary/endpoint' in _source
_source = _source.replace('out/lavapipe-intergroup-row-owner/ordinary/endpoint', 'out/lavapipe-owner-bulk-zero/ordinary/endpoint')
exec(compile(_source, str(_basis), "exec"))
