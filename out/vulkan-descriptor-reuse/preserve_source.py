from pathlib import Path
import hashlib, json, os, subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
tag = 'evidence/issue-97-descriptor-reuse-8a27321'
revision = json.loads((W / 'audited-benefit/retention-gate-frozen.json').read_text())['candidate_revision']
query = lambda args, **kw: subprocess.check_output(args, **kw)
assert query(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert query(['git', 'status', '--porcelain'], text=True) == ''
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
index = W / 'source.index'
assert not index.exists()
environment = dict(os.environ, GIT_INDEX_FILE=str(index))
sources = sorted(W.rglob('*.py')) + [W / 'audited-observer.cpp']
sources += sorted(W.glob('actions*.json'))
for folder in ['benefit', 'audited-benefit']:
    sources += [W / folder / name for name in ['actions.json', 'retention-gate-frozen.json', 'expected-readback-layout.json']]
assert all(p.is_file() for p in sources)
records = []
try:
    subprocess.run(['git', 'read-tree', revision], env=environment, check=True)
    for p in sources:
        relative = str(p.relative_to(R))
        blob = query(['git', 'hash-object', '-w', str(p)], text=True).strip()
        subprocess.run(['git', 'update-index', '--add', '--cacheinfo', '100644,' + blob + ',' + relative], env=environment, check=True)
        records.append({'path': relative, 'blob': blob, 'bytes': p.stat().st_size,
                        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    tree = query(['git', 'write-tree'], env=environment, text=True).strip()
    commit = query(['git', 'commit-tree', tree, '-p', revision], input='Preserve reproducible #97 descriptor binding qualification commands\n', text=True).strip()
    subprocess.run(['git', 'update-ref', 'refs/tags/' + tag, commit], check=True)
finally:
    index.unlink(missing_ok=True)
subprocess.run(['git', 'push', 'origin', 'refs/tags/' + tag], check=True)
remote = query(['git', 'ls-remote', 'origin', 'refs/tags/' + tag], text=True).strip()
assert remote == commit + '\trefs/tags/' + tag
assert query(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert query(['git', 'status', '--porcelain'], text=True) == ''
(W / 'preserved-source.json').write_text(json.dumps({'tag': tag, 'commit': commit, 'parent': revision,
    'source_records': records, 'remote_readback': remote, 'working_tree_and_head_unchanged': True}, indent=2) + '\n')
print(json.dumps({'tag': tag, 'commit': commit, 'source_files': len(records), 'pushed': True}))
