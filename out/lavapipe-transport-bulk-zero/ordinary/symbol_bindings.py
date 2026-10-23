"""Reuse exact ELF symbol binding; add only the exact Transport clear extent."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-owner-bulk-zero/ordinary/symbol_bindings.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='30075994e5d0de3f4fa9db381044a4c77ba6d4c4268e6060b495023b874da738'
_source=_basis.read_text();_old="{'Endpoint': 97920, 'DopriPhase': 100640}"
assert _source.count(_old)==1
_source=_source.replace(_old,"{'Transport': 73440, 'Endpoint': 97920, 'DopriPhase': 100640}")
exec(compile(_source,str(_basis),'exec'))
