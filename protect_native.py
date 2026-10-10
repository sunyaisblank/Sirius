"""Preserve completed native cases with recoverable immutable input joins."""
from pathlib import Path
import hashlib, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(path): return json.loads(path.read_text())
def identity(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
def key(value): return value['bytes'], value['sha256']
bindings = read(WORK / 'execution-bindings.json')
head = bindings['candidate_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
comparison = read(WORK / 'native-comparison.json')
assert not comparison['pass'] and len(comparison['checks']) == 20
assert len(comparison['failed_checks']) == 10
for name in ('controls', 'matched'):
    sequence = read(WORK / ('native-' + name + '-sequence.json'))
    assert sequence['completed'] and 'failure' not in sequence
    assert all(v['ownership_verified'] and not v['cleanup_errors'] for v in sequence['owned_calls'])

DEST = ROOT / 'attestations/native-vulkan' / head[:7] / 'submission-fence-rejected'
assert not DEST.exists()
recoveries = {}
software = ROOT / read(WORK / 'software-protection.json')['archive']
software_manifest = read(software / 'manifest.json')
for record in software_manifest['payloads']:
    path = software / record['path']
    value = identity(path)
    assert value == {k: record[k] for k in ('bytes', 'sha256')}
    recoveries[key(value)] = {'kind': 'retained-archive-payload', 'path': str(path.relative_to(ROOT))}
for revision in (bindings['baseline_revision'], head):
    directory = ROOT / 'attestations/native-build' / revision[:7]
    for path in directory.rglob('*'):
        assert not path.is_symlink()
        if path.is_file():
            recoveries[key(identity(path))] = {'kind': 'retained-native-build-payload', 'path': str(path.relative_to(ROOT))}

trees = {}
for revision in (bindings['baseline_revision'], head):
    tree = {}
    for row in subprocess.check_output(['git', 'ls-tree', '-r', '-z', revision]).split(b'\0'):
        if row:
            info, path = row.split(b'\t', 1)
            mode, kind, oid = info.decode().split()
            assert kind == 'blob'
            tree[path.decode()] = oid
    trees[revision] = tree

payloads = {}
def preserve(path, value):
    assert path.is_file() and not path.is_symlink()
    assert identity(path) == value, str(path)
    if key(value) in recoveries: return recoveries[key(value)]
    target = 'payloads/' + value['sha256']
    output = DEST / target
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, output)
    assert identity(output) == value
    payloads[target] = {'path': target, **value}
    recovery = {'kind': 'local-archive-payload', 'path': str(output.relative_to(ROOT))}
    recoveries[key(value)] = recovery
    return recovery

ledger = []
for origin, value in sorted(bindings['files'].items()):
    path = ROOT / origin
    assert identity(path) == value, origin
    revision, source_path = head, origin
    for candidate in (bindings['baseline_revision'], head):
        prefix = 'bin/windows-msvc/native-' + candidate[:7] + '/recorded-source/'
        if origin.startswith(prefix): revision, source_path = candidate, origin[len(prefix):]
    if source_path in trees[revision]:
        blob = trees[revision][source_path]
        assert subprocess.check_output(['git', 'hash-object', '--no-filters', str(path)], text=True).strip() == blob
        recovery = {'kind': 'immutable-Git-blob', 'revision': revision, 'blob': blob, 'path': source_path}
    else:
        recovery = preserve(path, value)
    ledger.append({'origin': origin, **value, 'recovery': recovery})

diagnostics = []
for path in sorted(WORK.rglob('*')):
    assert not path.is_symlink()
    if path.is_file():
        value = identity(path)
        diagnostics.append({'origin': str(path.relative_to(ROOT)), **value, 'recovery': preserve(path, value)})
for record in (*ledger, *diagnostics):
    recovery = record['recovery']
    expected = {k: record[k] for k in ('bytes', 'sha256')}
    if recovery['kind'] == 'immutable-Git-blob':
        data = subprocess.check_output(['git', 'cat-file', 'blob', recovery['blob']])
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected
    else: assert identity(ROOT / recovery['path']) == expected
    assert identity(ROOT / record['origin']) == expected
assert {str(p.relative_to(DEST)) for p in DEST.rglob('*') if p.is_file()} == set(payloads)
manifest = {'source_revision': head, 'baseline_revision': bindings['baseline_revision'],
            'executed_input_ledger': ledger, 'diagnostics': diagnostics,
            'payloads': list(payloads.values()), 'payload_count': len(payloads),
            'payload_bytes': sum(v['bytes'] for v in payloads.values()),
            'scope': 'Two original native controls and one frozen A1/B1/B2/A2 comparison. Benefit gate rejected; numerical and bounded ownership controls passed. No full frame, science estate, release or cleanup clearance.'}
(DEST / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
receipt = {'archive': str(DEST.relative_to(ROOT)), 'manifest': identity(DEST / 'manifest.json'),
           'payload_count': len(payloads), 'payload_bytes': manifest['payload_bytes'],
           'executed_input_count': len(ledger), 'diagnostic_count': len(diagnostics),
           'all_origin_and_recovery_bytes_verified': True, 'cleanup_authorized': False}
(WORK / 'native-protection.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
