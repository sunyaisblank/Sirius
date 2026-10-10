"""Preserve finite original science, failed initial build and exact input versions."""
from pathlib import Path
import hashlib, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def document(path):
    return json.loads(path.read_text())
def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
source = document(WORK / 'source.json')
head = source['source_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert document(WORK / 'independent-software-review.json')['passed']
DEST = ROOT / 'attestations/review' / head[:7] / 'fp64-science-observation'
assert not DEST.exists()
selected, ledgers, trees = {}, {}, {}
def select(path, target):
    assert path.is_file() and not path.is_symlink()
    assert target not in selected or identity(selected[target]) == identity(path)
    selected[target] = path
def tree(revision):
    if revision not in trees:
        result = {}
        for row in subprocess.check_output(['git', 'ls-tree', '-r', '-z', revision]).split(b'\0'):
            if row:
                info, path = row.split(b'\t', 1)
                mode, kind, oid = info.decode().split()
                assert kind == 'blob'
                result[path.decode()] = oid
        trees[revision] = result
    return trees[revision]
alternates = [WORK / n for n in ('actions-076adbc.json', 'actions-corrected-configure.json',
    'validate-observations-sealed-before-science.py')]
for case in ('configure', 'build', 'corrected-configure', 'corrected-build', 'science'):
    folder = WORK / case
    owner = document(folder / 'owner.json')
    assert owner['passed'] is (case != 'build')
    assert owner['child_exit'] == (1 if case == 'build' else 0)
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors']
    assert not owner['remaining_owned_processes'] and all(v['absent'] for v in owner['observed_births'])
    revision = owner['source']['revision']
    assert owner['source']['status'] == ''
    inputs = document(folder / 'inputs-before.json')
    assert inputs == document(folder / 'inputs-after.json')
    ledger = []
    for relative, seal in sorted(inputs.items()):
        origin = ROOT / relative
        assert seal['resolved_path'] == str(origin.resolve())
        value = {k: seal[k] for k in ('bytes', 'sha256')}
        if relative in tree(revision):
            data = subprocess.check_output(['git', 'cat-file', 'blob', tree(revision)[relative]])
            assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == value
            recovery = {'kind': 'immutable-Git-blob', 'revision': revision, 'blob': tree(revision)[relative], 'path': relative}
        else:
            path = origin
            if identity(path) != value:
                matches = [p for p in alternates if identity(p) == value]
                assert len(matches) == 1, (case, relative)
                path = matches[0]
            assert identity(path) == value
            target = 'executed-inputs/' + value['sha256']
            select(path, target)
            recovery = {'kind': 'local-archive-payload', 'path': target}
        ledger.append({'origin': relative, **value, 'recovery': recovery})
    ledgers[case] = {'source_revision': revision, 'passed': owner['passed'], 'inputs': ledger}
    for path in folder.rglob('*'):
        if path.is_file():
            select(path, 'diagnostic/' + str(path.relative_to(WORK)))
for name in ('source.json', 'source-076adbc.json', 'actions.json', 'actions-076adbc.json',
    'actions-corrected-configure.json', 'baseline.json', 'baseline-arrays.json',
    'expected-stream.json', 'build-and-arrays.json', 'linux-readonly-payloads.json',
    'source-review.json', 'independent-software-review.json', 'preregister.md', 'preregister-readback.json',
    'provider-before.json', 'provider-live.json', 'provider-accepted.json', 'run_owned.py', 'owned_guard.py',
    'bind_linux.py', 'observer.patch', 'candidate.patch', 'candidate-076adbc.patch',
    'validate-observations-sealed-before-science.py', 'validate-observations-rejected-encoding.py',
    'validate-observations-local-accepted.py', 'validator-encoding-correction.json', 'protect_software.py'):
    select(WORK / name, 'diagnostic/' + name)
payloads = []
for target, path in sorted(selected.items()):
    seal = identity(path)
    output = DEST / target
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, output)
    assert identity(output) == seal
    payloads.append({'path': target, 'origin': str(path.relative_to(ROOT)), **seal})
assert {str(p.relative_to(DEST)) for p in DEST.rglob('*') if p.is_file()} == {v['path'] for v in payloads}
for item in payloads:
    assert identity(DEST / item['path']) == {k: item[k] for k in ('bytes', 'sha256')}
manifest = {'revision': head, 'payloads': payloads, 'payload_count': len(payloads),
    'payload_bytes': sum(v['bytes'] for v in payloads), 'executed_cases': ledgers,
    'scope': 'Failed 076 initial strict build preserved as failed; corrected 95 strict configure/build and one original software science pass. 160 unchanged complete readonly ELF payload joins. Recoverable original frozen input versions. No native/frame/full qualification claim.'}
(DEST / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
receipt = {'archive': str(DEST.relative_to(ROOT)), 'manifest': identity(DEST / 'manifest.json'),
    'payload_count': len(payloads), 'payload_bytes': manifest['payload_bytes'],
    'case_input_counts': {k: len(v['inputs']) for k, v in ledgers.items()}, 'fresh_byte_checks_pass': True}
(WORK / 'software-protection.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
