"""Reuse original gated frame link/runbook with unchanged mixed backend."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-cpu-row-ownership/cold-frame/prepare.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '5c702abf6582dba1a79b47628576db51c96bd6544f8bac46136f562de7f34f6f'
_source = _basis.read_text()
assert "B = C / 'native-cpu-row-ownership'" in _source
_source = _source.replace("B = C / 'native-cpu-row-ownership'", "B = C / 'lavapipe-owner-machine-codegen'")
assert "SOURCE = 'bfa095baf0ff9592cf648cfd49a73d5f6814c05c'" in _source
_source = _source.replace("SOURCE = 'bfa095baf0ff9592cf648cfd49a73d5f6814c05c'", "SOURCE = json.loads((P / 'runtime/actions.json').read_text())['control-8']['revision']")
assert "B / 'matched/libsirius_backend.a'" in _source
_source = _source.replace("B / 'matched/libsirius_backend.a'", "C / 'native-cooperative-transport-intergroup/matched/libsirius_backend.a'")
assert 'renders/native-cpu-row-bfa095ba/cold-frame-tmp' in _source
_source = _source.replace('renders/native-cpu-row-bfa095ba/cold-frame-tmp', 'renders/lavapipe-owner-machine-codegen/cold-frame-tmp')
exec(compile(_source, str(_basis), 'exec'))
