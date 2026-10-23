"""Use original predicates with the sealed regular module catalogue."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-acqrel-ordering/readers/read_results.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'e4792eac2c0505f3870fb1c3b8b0095221fedcfca340b78d9aa990398a6b8974'
_source = _basis.read_text()
assert 'out/lavapipe-acqrel-ordering' in _source
_source = _source.replace('out/lavapipe-acqrel-ordering', 'out/lavapipe-owner-machine-codegen')
assert 'out/native-original-input-ownership/compiler-result.json' in _source
_source = _source.replace('out/native-original-input-ownership/compiler-result.json', 'out/native-cooperative-transport-intergroup/compiler-result.json')
exec(compile(_source, str(_basis), 'exec'))
