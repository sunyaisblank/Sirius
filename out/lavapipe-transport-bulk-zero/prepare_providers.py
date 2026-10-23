"""Bind the actual provider and retain explicit current-versus-reused code scope."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-intergroup-row-owner/prepare_providers.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='89d0bd88a904da0bcb11c64f1f23fc63c4f1b4391118a80e3e4a2ff7fe3ee24b'
_source=_basis.read_text();_old='all-three current code/features/properties/MXCSR'
assert _source.count(_old)==1
_source=_source.replace(_old,'current Transport code/features/properties/MXCSR, authenticated unchanged E/DP source/prior-code composition')
exec(compile(_source,str(_basis),'exec'))
