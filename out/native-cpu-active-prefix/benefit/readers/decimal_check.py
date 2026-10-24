"""Reuse the sealed original authority; path binding changes only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-original-input-ownership/readers/decimal_check.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'd27a3e94775edbdba42b48569d0e09a794b720fc13bc94d92787f9c62b1f305b'
source = basis.read_text()
changes = [("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/native-cpu-active-prefix/benefit/runtime'")]
for before, after in changes:
    assert source.count(before) >= 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
