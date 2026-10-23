"""Stage the one new provider; reuse original configured build and product seals."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-owner-bulk-zero/build_incremental.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='67fa6a2a279f2b166fa4a7a5562cfc4727fd2c46008819c8f6e60bbc37fc8e9e'
_source=_basis.read_text();_old='bin/linux-gcc/lavapipe-owner-bulk-zero/provider'
assert _source.count(_old)==1
_source=_source.replace(_old,'bin/linux-gcc/lavapipe-transport-bulk-zero/provider')
exec(compile(_source,str(_basis),'exec'))
