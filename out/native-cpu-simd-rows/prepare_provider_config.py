"""Reuse original provider configuration, sealing the new adapter and SIMD objects."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_config.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5ce6ffb5e7ff923cf0a34696ec088132915e90b57205b9da4cbaad7471c04762'
source = basis.read_text()
for before, after in [
    ("(W / 'build/actions.json').read_text())['typed-three-stage-batch']", "(W / 'compile/actions.json').read_text())['adapter']"),
    ("(W / 'batch-compile-result.json').read_text()", "(W / 'adapter-compile-result.json').read_text()"),
    ('consume3 sealed baseline-ISA native objects', 'consume8 sealed baseline-ISA native objects'),
]:
    assert before in source
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
