"""Reuse the sealed incremental provider builder with the scoped codegen patch."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-intergroup-row-owner/build_incremental.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'd53eabc6275fa6d664ad2180208d405d131c2b72bb3596d030b41d81c70a8a47'
_source = _basis.read_text()
assert 'bin/linux-gcc/lavapipe-intergroup-row-owner/provider' in _source
_source = _source.replace('bin/linux-gcc/lavapipe-intergroup-row-owner/provider', 'bin/linux-gcc/lavapipe-owner-machine-codegen/provider')
assert 'One incremental private provider build with unchanged original configuration and sealed reused artifacts.' in _source
_source = _source.replace('One incremental private provider build with unchanged original configuration and sealed reused artifacts.', 'One incremental scoped post-IR machine-codegen provider build with unchanged original configuration and sealed reused artifacts.')
exec(compile(_source, str(_basis), "exec"))
