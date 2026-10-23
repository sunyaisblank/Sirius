"""Bind actual new provider only after its complete code gate."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-intergroup-row-owner/prepare_providers.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '89d0bd88a904da0bcb11c64f1f23fc63c4f1b4391118a80e3e4a2ff7fe3ee24b'
exec(compile(_basis.read_text(), str(_basis), 'exec'))
