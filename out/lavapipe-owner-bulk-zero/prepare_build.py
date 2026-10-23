"""Reuse the sealed incremental provider builder with the scoped codegen patch."""
from pathlib import Path
import hashlib, json
_root = Path.cwd()
_workspace = Path(__file__).resolve().parent
_mesa = _root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
_review = json.loads((_workspace / 'review-admission.json').read_text())
_current = {str(p.relative_to(_mesa)): {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
    for p in _mesa.rglob('*') if p.is_file()}
_canonical = json.dumps(_current, sort_keys=True, separators=(',', ':')).encode()
assert {'bytes': len(_canonical), 'sha256': hashlib.sha256(_canonical).hexdigest()} == _review['full_current_source_identity']
_basis = Path.cwd() / 'out/lavapipe-intergroup-row-owner/prepare_build.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '45caa976920ec3c020d5387c4a7486f9ebdb11a7f5d15842c1ecef82b61c01be'
_source = _basis.read_text()
assert "'prototype.patch'" in _source
_source = _source.replace("'prototype.patch'", "'bulk-zero.patch'")
assert "'prototype_patch'" in _source
_source = _source.replace("'prototype_patch'", "'prototype_patch'")
assert 'One incremental private provider build from independently reviewed intergroup row-owner source.' in _source
_source = _source.replace('One incremental private provider build from independently reviewed intergroup row-owner source.', 'One incremental private provider build with equivalent bulk scratch initialization for exact intergroup owners and restored Default machine optimization.')
exec(compile(_source, str(_basis), "exec"))
