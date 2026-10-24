"""Reuse G's exact strict adapter recipe with independent I object output."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-simd-rows/prepare_adapter.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'de6a248b6aefa4e33d54f6000e80a6c9732a6e218e10ef4d7430e248fab24ced'
source = basis.read_text()
before = "output = root / 'bin/linux-gcc/native-cpu-simd-rows/objects/rows_adapter.o'"
after = "output = root / 'bin/linux-gcc/native-cpu-simd-packet-policy/objects/rows_adapter.o'\noutput.parent.mkdir(parents=True, exist_ok=True)"
assert source.count(before) == 1
exec(compile(source.replace(before, after), str(basis), 'exec'))
