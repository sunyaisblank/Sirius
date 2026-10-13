"""Preserve reusable issue98 integrated control tools without changing production HEAD."""
from pathlib import Path
import hashlib, json, os, subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
tag = 'evidence/issue-98-guarded-product-high-controls-4bcd987'
revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
query = lambda args, **kw: subprocess.check_output(args, **kw)
assert query(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert query(['git', 'status', '--porcelain'], text=True) == ''
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
index = W / 'source.index'
assert not index.exists()
environment = dict(os.environ, GIT_INDEX_FILE=str(index))
sources = sorted(W.glob('*.py')) + [W / item for item in ['audited-observer.cpp', 'actions.json', 'README.md']]
historical_preparers = ['81431fd39d42e31270a4b3fd35e16d0ee0d69705a5c49567b657591f1245220d',
                       'c19b80943485c298e99eaf6ff02f08e489b858dec7f01eed8e0dde93c3a65749']
sources += [W / 'input-versions' / value for value in historical_preparers]
assert all(p.is_file() for p in sources)
records = []
try:
    subprocess.run(['git', 'read-tree', revision], env=environment, check=True)
    for p in sources:
        suffix = ('source-versions/prepare-' + p.name[:12] + '.py'
                  if p.name in historical_preparers else p.name)
        relative = 'diagnostics/guarded-product-high-controls/' + suffix
        origin = str(p.relative_to(R))
        blob = query(['git', 'hash-object', '-w', str(p)], text=True).strip()
        subprocess.run(['git', 'update-index', '--add', '--cacheinfo', '100644,' + blob + ',' + relative], env=environment, check=True)
        records.append({'path': relative, 'origin': origin, 'blob': blob, 'bytes': p.stat().st_size,
                        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    tree = query(['git', 'write-tree'], env=environment, text=True).strip()
    commit = query(['git', 'commit-tree', tree, '-p', revision], input='Preserve issue #98 matched original integrated-control producer and verification\n', text=True).strip()
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
