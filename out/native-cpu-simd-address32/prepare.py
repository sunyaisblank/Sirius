"""Seal the sole original C8 against the genuine eight-object SIMD provider."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6a414a90ea995c5cfb93d3c132c91397fe9203160d930f80125d2f8785ed85cb'
source = basis.read_text()
for before, after in [
    ("native = root / 'out/native-cpu-active-prefix'", "native = root / 'out/native-cpu-simd-address32'"),
    ("doc(native / 'batch-compile-result.json')", "doc(native / 'adapter-compile-result.json')"),
    ('nativeCPU active-prefix providerB', 'nativeCPU SIMD rows providerH'),
    ('First actual native packet/return perstage observed passively;',
     'First actual rows wrapper, pure preflight and raw SIMD entry/return perstage observed passively, with exact concurrent parent/range joins;'),
]:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
