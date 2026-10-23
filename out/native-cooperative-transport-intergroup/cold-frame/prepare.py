"""Reuse the original gated link/frame protocol with the composed backend."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cpu-row-ownership/cold-frame/prepare.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '5c702abf6582dba1a79b47628576db51c96bd6544f8bac46136f562de7f34f6f'
_source = _basis.read_text()
assert _source.count("B = C / 'native-cpu-row-ownership'") == 1
_source = _source.replace("B = C / 'native-cpu-row-ownership'", "B = C / 'native-cooperative-transport-intergroup'")
_source = _source.replace("SOURCE = 'bfa095baf0ff9592cf648cfd49a73d5f6814c05c'", "SOURCE = 'b1048e57c613ff47c959b704932cc0ecaaadaaf1'")
_source = _source.replace('renders/native-cpu-row-bfa095ba/cold-frame-tmp', 'renders/native-cooperative-transport-intergroup-b1048e57/cold-frame-tmp')
exec(compile(_source, str(_basis), 'exec'))
