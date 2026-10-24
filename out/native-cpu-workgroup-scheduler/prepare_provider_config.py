"""Seal the original compiled D objects; no duplicate native compile."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare_provider_config.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '333cf237b0cf71129475b4c3e3c1c6eee9891c6e969e8f04891de062607ff81d'
source = basis.read_text()
# Apply the sealed B adapter to the original preparer, then bind genuine D receipts.
original = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_config.py'
assert hashlib.sha256(original.read_bytes()).hexdigest() == '5ce6ffb5e7ff923cf0a34696ec088132915e90b57205b9da4cbaad7471c04762'
source = original.read_text()
changes = [
    ("(W / 'build/actions.json').read_text())['typed-three-stage-batch']", "(R / 'out/native-cpu-static-program/build/actions.json').read_text())['compile']"),
    ("(W / 'batch-compile-result.json').read_text()", "(R / 'out/native-cpu-static-program/batch-compile-result.json').read_text()"),
]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(original), 'exec'))
