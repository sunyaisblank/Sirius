"""Reuse the sealed original authority; path binding changes only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-original-input-ownership/readers/run_abba.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '14177971258ceb4a269549f1c0f28055dec6d165f83e46f6047953e53be0d7eb'
source = basis.read_text()
changes = [("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/native-cpu-active-prefix/benefit/runtime'"), ("acceptance=W/'control-8-acceptance.json'", "acceptance=R/'out/native-cpu-active-prefix/control-8-acceptance.json'")]
for before, after in changes:
    assert source.count(before) >= 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
