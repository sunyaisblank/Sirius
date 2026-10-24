"""Seal original C8 against genuine L SIMD, N adapter and N provider bytes."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6a414a90ea995c5cfb93d3c132c91397fe9203160d930f80125d2f8785ed85cb'
source = basis.read_text()
changes = [
    ("native = root / 'out/native-cpu-active-prefix'", "native = root / 'out/native-cpu-simd-medium-tails'"),
    ("doc(native / 'batch-compile-result.json')", "doc(root / 'out/native-cpu-simd-endpoint-regions/ispc-uniform-compile-result.json')"),
    ("expected.update(compiled['source_inputs'])",
     "expected.update(compiled['source_inputs'])\n"
     "expected.update(doc(native / 'adapter-compile-result.json')['source_inputs'])\n"
     "expected.update(doc(native / 'provider/source-contract.json')['source_inputs'])\n"
     "for name in ('ispc-uniform-compile-result.json', 'compile/ispc-uniform/owner.json'):\n"
     "    path = root / 'out/native-cpu-simd-endpoint-regions' / name\n"
     "    expected[str(path)] = pin(path)\n"
     "for name in ('adapter-compile-result.json', 'compile/adapter/owner.json'):\n"
     "    path = native / name\n"
     "    expected[str(path)] = pin(path)"),
    ('nativeCPU active-prefix providerB', 'nativeCPU medium-tail providerN with7 unchanged D/F/L objects and1 strict N adapter'),
    ('First actual native packet/return perstage observed passively;',
     'First actual rows wrapper, pure preflight and raw SIMD4..8 entry/return perstage observed passively, with exact concurrent parent/range joins and no partial-packet resampling;'),
]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
