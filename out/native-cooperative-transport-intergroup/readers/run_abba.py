"""Reuse the original acceptance authority with current exact artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-acqrel-ordering/readers/run_abba.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '32273d8ad5c85669a8c548c24921bcba3131182113232680adc74a74a1ebaf90'
_source = _basis.read_text()
assert 'out/lavapipe-acqrel-ordering' in _source
_source = _source.replace('out/lavapipe-acqrel-ordering', 'out/native-cooperative-transport-intergroup')
exec(compile(_source, str(_basis), "exec"))
