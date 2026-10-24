"""Reuse original provider configuration, sealing the new adapter and SIMD objects."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_config.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5ce6ffb5e7ff923cf0a34696ec088132915e90b57205b9da4cbaad7471c04762'
source = basis.read_text()
for before, after in [
    ("(W / 'build/actions.json').read_text())['typed-three-stage-batch']", "(W / 'compile/actions.json').read_text())['ispc-uniform']"),
    ("(W / 'batch-compile-result.json').read_text()", "(W / 'ispc-uniform-compile-result.json').read_text()"),
    ("(W / 'build' / name).read_text()", "(W / 'compile' / name).read_text()"),
    ('consume3 sealed baseline-ISA native objects', 'consume8 sealed native objects under the unchanged CPU ISA guard'),
]:
    assert before in source
    source = source.replace(before, after)
before = "for item in json.loads((W / 'ispc-uniform-compile-result.json').read_text())['objects']:\n    inputs.add(item['object']['path'])"
after = "for item in json.loads((P / 'source-contract.json').read_text())['native_objects']:\n    inputs.add(item['path'])"
assert source.count(before) == 1
source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
