"""Reuse the current Mesa cache/options; seal the eleven real native objects."""
from pathlib import Path

basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_config.py'
source = basis.read_text()
changes = [
    ("(W / 'build/actions.json').read_text())['typed-three-stage-batch']", "(W / 'compile/actions.json').read_text())['cpp']"),
    ("inputs = set(base['frozen_inputs'])", "inputs = set(base['frozen_inputs'])\ninputs |= set(json.loads((P / 'source-contract.json').read_text())['source_inputs'])"),
    ("for item in json.loads((W / 'batch-compile-result.json').read_text())['objects']:\n    inputs.add(item['object']['path'])", "for item in json.loads((P / 'source-contract.json').read_text())['native_objects']:\n    inputs.add(item['path'])"),
    ("(W / 'build' / name).read_text()", "(W / 'compile' / name).read_text()"),
    ('consume3 sealed baseline-ISA native objects', 'link8 sealed objects with three actual new scalar-zero objects and unchanged S scheduling')]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
