"""Reuse the original owned serial compiler with H-only object outputs."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-simd-rows/build_ispc.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'cc7545ec60842366540e185adaf931c09fac0ed2b3dd1b23168f2ef2aa32a0d9'
source = basis.read_text()
before = "destination = root / 'bin/linux-gcc/native-cpu-simd-rows/objects-uniform'"
after = "destination = root / 'bin/linux-gcc/native-cpu-simd-address32/objects-uniform'"
assert source.count(before) == 1
exec(compile(source.replace(before, after), str(basis), 'exec'))
