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
_source = _source.replace("'prototype.patch'", "'machine-codegen.patch'")
assert "'prototype_patch'" in _source
_source = _source.replace("'prototype_patch'", "'machine_codegen_patch'")
assert 'One incremental private provider build from independently reviewed intergroup row-owner source.' in _source
_source = _source.replace('One incremental private provider build from independently reviewed intergroup row-owner source.', 'One incremental private provider build of the independently reviewed post-IR machine CodeGenNone policy for exact intergroup owners.')
exec(compile(_source, str(_basis), "exec"))
