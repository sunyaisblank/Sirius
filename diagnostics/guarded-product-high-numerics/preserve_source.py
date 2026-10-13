"""Preserve reusable issue98 numerical tools without changing production HEAD."""
from pathlib import Path
import hashlib, json, os, subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
tag = 'evidence/issue-98-guarded-product-high-numerics-4bcd987'
revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
query = lambda args, **kw: subprocess.check_output(args, **kw)
assert query(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert query(['git', 'status', '--porcelain'], text=True) == ''
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
index = W / 'source.index'
assert not index.exists()
environment = dict(os.environ, GIT_INDEX_FILE=str(index))
sources = sorted(W.glob('*.py')) + [W / name for name in ['residual_probe.slang', 'limb_arithmetic_probe.slang', 'actions.json', 'README.md']]
sources += [R / 'bin/linux-gcc/guarded-product-high-numerics/candidate/kernels/portable_binary32.h']
assert all(p.is_file() for p in sources)
records = []
try:
    subprocess.run(['git', 'read-tree', revision], env=environment, check=True)
    for p in sources:
        relative = 'diagnostics/guarded-product-high-numerics/' + ('candidate/' if p.name == 'portable_binary32.h' else '') + p.name
        origin = str(p.relative_to(R))
        blob = query(['git', 'hash-object', '-w', str(p)], text=True).strip()
        subprocess.run(['git', 'update-index', '--add', '--cacheinfo', '100644,' + blob + ',' + relative], env=environment, check=True)
        records.append({'path': relative, 'origin': origin, 'blob': blob, 'bytes': p.stat().st_size,
                        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    tree = query(['git', 'write-tree'], env=environment, text=True).strip()
    commit = query(['git', 'commit-tree', tree, '-p', revision], input='Preserve issue #98 actual arithmetic probes and original numerical validators\n', text=True).strip()
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
