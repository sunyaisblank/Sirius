"""Read-only comparison of G's exact uniform ISPC source and copied receipt."""
from pathlib import Path
import hashlib
import sys

assert not sys.argv[1:], 'H reuses existing G source; generation is forbidden'
basis = Path.cwd() / 'out/native-cpu-simd-rows/derive_ispc.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '8875e5bb3978e7c3072a3d053cc3ab215368d07995411be22cb4676573363722'
exec(compile(basis.read_text(), str(basis), 'exec'))
