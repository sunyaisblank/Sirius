"""Bind the new scheduler provider to genuine, unchanged D object receipts."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6a414a90ea995c5cfb93d3c132c91397fe9203160d930f80125d2f8785ed85cb'
source = basis.read_text()
changes = [
    ("native = root / 'out/native-cpu-active-prefix'", "native = root / 'out/native-cpu-workgroup-scheduler'"),
    ("doc(native / 'batch-compile-result.json')", "doc(root / 'out/native-cpu-static-program/batch-compile-result.json')"),
    ('nativeCPU active-prefix providerB', 'nativeCPU workgroup scheduler providerE'),
]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
