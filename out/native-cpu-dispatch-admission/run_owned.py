"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/run_owned.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6ae4510d2a1672ffb6b1fee3e37bdbc62c3aee3e993cafb1c210224b55cf2dd2'
source = basis.read_text()
changes = []
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
