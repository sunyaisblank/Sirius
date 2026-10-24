"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/build_provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '3eb7e466d8a3d1cbfc9f17d9dcfe4b65557318bc4b24ec77bdc2af9b042cd1b5'
source = basis.read_text()
changes = [('native-cpu-active-prefix/provider', 'native-cpu-dispatch-admission/provider')]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
