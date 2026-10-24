"""Seal the existing Mesa cache and all eight genuine reused L objects."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_config.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5ce6ffb5e7ff923cf0a34696ec088132915e90b57205b9da4cbaad7471c04762'
source = basis.read_text()
changes = [
    ("(W / 'build/actions.json').read_text())['typed-three-stage-batch']",
     "(R / 'out/native-cpu-simd-endpoint-regions/compile/actions.json').read_text())['ispc-uniform']"),
    ("inputs = set(base['frozen_inputs'])",
     "inputs = set(base['frozen_inputs'])\n"
     "inputs |= set(json.loads((P / 'source-contract.json').read_text())['source_inputs'])"),
    ("for item in json.loads((W / 'batch-compile-result.json').read_text())['objects']:\n    inputs.add(item['object']['path'])",
     "for item in json.loads((P / 'source-contract.json').read_text())['native_objects']:\n    inputs.add(item['path'])"),
    ("(W / 'build' / name).read_text()", "(W / 'compile' / name).read_text()"),
    ('consume3 sealed baseline-ISA native objects', 'reuse all8 sealed L native objects while adding M dead-SSA cleanup'),
]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
