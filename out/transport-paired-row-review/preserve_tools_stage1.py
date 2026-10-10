"""Preserve prototype/control source in a separate evidence tag, not adoption."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
TAG = 'evidence/issue-90-source-numerical-6845c81'
INDEX = TASK / 'source-git.index'
EXCLUDED = {'numerical-prepare-initial', 'numerical-prepare-reviewed-first',
            'native-temp', 'query-owner-temp', '__pycache__'}
baseline = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
assert baseline == '6845c8150ebf0df6992fe2609a4e9323c7f0612b'
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
assert not INDEX.exists()
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + TAG], cwd=ROOT).returncode == 1
sources = []
for path in sorted(TASK.rglob('*')):
    relative = path.relative_to(TASK)
    if not path.is_file() or EXCLUDED.intersection(relative.parts):
        continue
    if path.suffix in ('.py', '.ps1', '.cs') or path.name in ('retained_transport.slang', 'retained_scratch.slang'):
        sources.append(path)
environment = dict(os.environ, GIT_INDEX_FILE=str(INDEX))


def git(*args, **kwargs):
    return subprocess.check_output(['git', *args], cwd=ROOT, env=environment, **kwargs)


try:
    git('read-tree', baseline)
    git('add', '-f', '--', *(str(p.relative_to(ROOT)) for p in sources))
    tree = git('write-tree', text=True).strip()
    commit = git('commit-tree', tree, '-p', baseline, input='Preserve paired Transport prototype and native gate tools (#90)\n', text=True).strip()
    git('tag', TAG, commit)
    git('push', 'origin', 'refs/tags/' + TAG)
    assert subprocess.check_output(['git', 'rev-parse', TAG], cwd=ROOT, text=True).strip() == commit
    remote = subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + TAG], cwd=ROOT, text=True).split()[0]
    assert remote == commit
    records = [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size,
                'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
    receipt = {'tag': TAG, 'commit': commit, 'base': baseline, 'source_count': len(sources),
               'sources': records, 'pushed_verified': True,
               'scope': 'Source-only evidence preservation; development HEAD unchanged, prototype unadopted.'}
    (TASK / 'source-git-preservation.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'sources'}, indent=2))
finally:
    INDEX.unlink(missing_ok=True)
