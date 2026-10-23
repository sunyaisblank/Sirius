"""Reuse the original acceptance authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-acqrel-ordering/readers/read_controls.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '4db9e2f8cc4581924b95ab08c70a3465f0229d900e7da4e20aeea2fc7b82ef6d'
_source = _basis.read_text()
assert 'out/native-original-input-ownership/compiler-result.json' in _source
_source = _source.replace('out/native-original-input-ownership/compiler-result.json', 'out/native-cpu-row-ownership/compiler-result.json')
exec(compile(_source, str(_basis), "exec"))
