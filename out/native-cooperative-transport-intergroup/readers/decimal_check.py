"""Reuse the original acceptance authority with current exact artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-acqrel-ordering/readers/decimal_check.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'b9da8e9a9b551e831ad3afdc54095061add9aa84bd1e70fd78436b97217d64f1'
_source = _basis.read_text()
assert 'out/lavapipe-acqrel-ordering' in _source
_source = _source.replace('out/lavapipe-acqrel-ordering', 'out/native-cooperative-transport-intergroup')
exec(compile(_source, str(_basis), "exec"))
