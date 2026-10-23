"""Reuse exact original acceptance authorities; bind the changed provider separately."""
from pathlib import Path
import hashlib, json, os
R = Path.cwd(); W = Path(__file__).resolve().parent
OLD = R / 'out/native-cooperative-transport-intergroup'
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def put(p, data):
    assert not p.exists(); p.parent.mkdir(parents=True, exist_ok=True); p.write_text(data)
def adapter(target, basis, replacements):
    lines = ['"""Reuse a sealed original authority with exact current artifact paths."""',
        'from pathlib import Path', 'import hashlib',
        '_basis = Path.cwd() / ' + repr(str(basis.relative_to(R))),
        'assert hashlib.sha256(_basis.read_bytes()).hexdigest() == ' + repr(pin(basis)['sha256']),
        '_source = _basis.read_text()']
    for old, new in replacements:
        lines += ['assert ' + repr(old) + ' in _source', '_source = _source.replace(' + repr(old) + ', ' + repr(new) + ')']
    lines += ['exec(compile(_source, str(_basis), "exec"))', '']
    put(target, '\n'.join(lines))
adapter(W / 'run_case.py', OLD / 'run_case.py', [])
target = W / 'runtime/run_owned.py'; target.parent.mkdir(exist_ok=True)
put(target, (OLD / 'runtime/run_owned.py').read_text())
target = W / 'runtime/owned_guard.py'
assert not target.exists(); target.symlink_to('../../native-cpu-row-ownership/cold-frame/owned_guard.py')
for name in ['compiler-result.json', 'matched-result.json']:
    target = W / name; assert not target.exists(); target.symlink_to('../native-cooperative-transport-intergroup/' + name)
for name in ['native_evidence.py', 'numeric_gates.py']:
    target = W / 'readers' / name; target.parent.mkdir(exist_ok=True)
    assert not target.exists(); target.symlink_to('../../native-original-input-ownership/readers/' + name)
put(W / 'readers/expected-readback-layout.json', (OLD / 'readers/expected-readback-layout.json').read_text())
for name in ['read_results.py', 'run_abba.py', 'decimal_check.py']:
    adapter(W / 'readers' / name, OLD / 'readers' / name,
        [('out/native-cooperative-transport-intergroup', 'out/lavapipe-owner-machine-codegen')])
adapter(W / 'readers/read_controls.py', OLD / 'readers/read_controls.py', [
    ('out/native-cooperative-transport-intergroup/compiler-result.json', 'out/lavapipe-owner-machine-codegen/compiler-result.json'),
    ('out/lavapipe-intergroup-row-owner/ordinary/endpoint', 'out/lavapipe-owner-machine-codegen/ordinary/transport')])
print(json.dumps({'original_reader_authorities_reused': True, 'acceptance_predicates_changed': False}))
