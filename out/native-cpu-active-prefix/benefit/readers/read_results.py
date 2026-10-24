"""Reuse the sealed original authority; path binding changes only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-original-input-ownership/readers/read_results.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'e7ef5a184fd2f25f5d0d606baadcae0e8ec8fb2f4e14e88c41e24915b556034e'
source = basis.read_text()
changes = [("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/native-cpu-active-prefix/benefit/runtime'"), ("doc(W.parent/'compiler-result.json')", "doc(R/'out/native-cpu-row-ownership/compiler-result.json')"), ("folder=RUN/label;owner=doc(folder/'owner.json')", "folder=(R/'out/native-cpu-active-prefix/control-8' if label=='control-8' else RUN/label);owner=doc(folder/'owner.json')"), ("doc(W/'control-8-acceptance.json')", "doc(R/'out/native-cpu-active-prefix/control-8-acceptance.json')"), ("pin(RUN/'control-8/owner.json')", "pin(R/'out/native-cpu-active-prefix/control-8/owner.json')"), ("doc(W.parent/'previous-control-8-acceptance.json')", "doc(R/'out/native-cpu-row-ownership/previous-control-8-acceptance.json')")]
for before, after in changes:
    assert source.count(before) >= 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
