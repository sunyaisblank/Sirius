"""Reuse G's actual-code checks on H's separately sealed objects."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-simd-rows/inspect_objects.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '2f4572906c6f266ba0b20fe553ebb0844cf08c8c8e933ea3902554206e5e212e'
exec(compile(basis.read_text(), str(basis), 'exec'))
