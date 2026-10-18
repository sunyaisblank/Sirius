import datetime, hashlib, json, os, pathlib, stat, subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT)
BASE = '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'
PARENT = '8554fd1969a1cd77e6393d34508550e19e7bfc7a'
BRANCH = 'refs/heads/development/native-fp32-packed-scratch-runtime-2026-10-18'
TAG = 'refs/tags/evidence/native-constant-controls-benefit-pass-5b2df58b'
BUNDLE = pathlib.Path('attestations/qualification/5b2df58b/native-constant-controls-benefit-pass')
RECEIPTS = pathlib.Path('out/native-constant-preservation')

def run(*args, data=None, env=None):
    return subprocess.run(args, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          check=True, env=env).stdout

def tree(rev):
    result = {}
    for row in run('git', 'ls-tree', '-rz', rev).split(b'\0'):
        if row:
            meta, name = row.split(b'\t')
            mode, kind, blob = meta.decode().split()
            assert kind == 'blob'
            result[name.decode()] = (mode, blob)
    return result

assert run('git', 'rev-parse', 'HEAD').decode().strip() == BASE
assert run('git', 'rev-parse', BRANCH).decode().strip() == PARENT
assert not run('git', 'status', '--porcelain')
canonical = tree(BASE)
assert len(canonical) == 538
physical_mode_differences = {}
for name, (mode, blob) in canonical.items():
    p = pathlib.Path(name)
    assert run('git', 'hash-object', '--', name).decode().strip() == blob
    physical_mode = '100755' if p.stat().st_mode & 0o111 else '100644'
    if physical_mode != mode:
        physical_mode_differences[name] = dict(Git=mode, physical=physical_mode)
assert physical_mode_differences == {'scripts/verify-operating-model.py': dict(Git='100755', physical='100644')}

index = (ROOT / RECEIPTS / 'preservation.index')
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
files = sorted(p for p in BUNDLE.rglob('*') if p.is_file())
records = {}
try:
    run('git', 'read-tree', PARENT, env=env)
    index_rows = []
    for p in files:
        assert not p.is_symlink()
        b = p.read_bytes()
        blob = run('git', 'hash-object', '-w', '--stdin', data=b).decode().strip()
        mode = '100755' if p.stat().st_mode & 0o111 else '100644'
        index_rows.append(f'{mode} {blob}\t{p}\0'.encode())
        records[str(p)] = dict(mode=mode, blob=blob, bytes=len(b), sha256=hashlib.sha256(b).hexdigest())
    run('git', 'update-index', '-z', '--index-info', data=b''.join(index_rows), env=env)
    new_tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', new_tree, '-p', PARENT,
                 data=b'Preserve native constant constructor original controls and 22-gate pass\n').decode().strip()
    actual = tree(commit)
    assert all(actual[n] == pair for n, pair in canonical.items())
    seen = set()
    for n, r in records.items():
        assert actual[n] == (r['mode'], r['blob'])
        if r['blob'] not in seen:
            b = run('git', 'cat-file', 'blob', r['blob'])
            assert len(b) == r['bytes'] and hashlib.sha256(b).hexdigest() == r['sha256']
            seen.add(r['blob'])
    run('git', 'update-ref', BRANCH, commit, PARENT)
    run('git', 'update-ref', TAG, commit, '0' * 40)
    receipt = dict(schema='sirius-native-constant-controls-benefit-Git-v1',
                   created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   commit=commit, parent=PARENT, branch=BRANCH, tag=TAG, files=records,
                   canonical_538_files_unchanged=True, stored_blob_bytes_verified=True, pushed=False)
    receipt['preexisting_checkout_physical_mode_difference'] = physical_mode_differences
    (RECEIPTS / 'root-actual-Git-preservation.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k:v for k,v in receipt.items() if k != 'files'}))
finally:
    index.unlink(missing_ok=True)
