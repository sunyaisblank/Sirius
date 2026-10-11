"""Push useful closing controller source and concise reviews without changing HEAD."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, subprocess
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def git(*args, **kw): return subprocess.check_output(['git', *args], **kw).decode().strip()
head = read(WORK / 'restoration/source.json')['source_revision']
assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}') and not git('status', '--porcelain')
review = read(WORK / 'cleanup-source-review.json'); assert review['accepted']
paths = [
    'protect_restoration.py', 'run_restoration_protection_owned.py', 'restoration-protection-actions.json',
    'restoration-protection-source-review.json', 'restoration-independent-review.json',
    'restoration-preservation-reviewed.json', 'fresh_native_consumers.ps1', 'run_fresh_census.py',
    'finish_cleanup.py', 'preserve_closing_tools.py', 'cleanup-source-review.json', 'cleanup-targets.json',
    'restoration/source.json', 'restoration/actions.json', 'restoration/selection-actions.json',
    'restoration/owned_guard.py', 'restoration/run_owned.py', 'restoration/run_selection_owned.py',
    'restoration/verify_restoration.py', 'restoration/bind_linux.py',
    'restoration/admission-registration-correction.json',
]
for name, value in review['controller_seals'].items():
    p = WORK / name
    assert {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} == value
tag = 'evidence/issue-95-rejected-trial-cleanup-' + head[:7]
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode != 0
assert not git('ls-remote', 'origin', 'refs/tags/' + tag)
index = WORK / 'closing-tools.index'; assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
subprocess.run(['git', 'read-tree', '--empty'], env=env, check=True)
records = []
for name in sorted(paths):
    p = WORK / name; assert p.is_file() and not p.is_symlink()
    data = p.read_bytes(); blob = git('hash-object', '-w', '--stdin', input=data)
    subprocess.run(['git', 'update-index', '--add', '--cacheinfo', '100644,' + blob + ',' + name], env=env, check=True)
    records.append({'path': name, 'origin': str(p.relative_to(ROOT)), 'blob': blob, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
tree = git('write-tree', env=env)
commit = git('commit-tree', tree, input=b'Preserve reviewed rejected #95 restoration and cleanup tools at 6be3132\n')
subprocess.run(['git', 'tag', tag, commit], check=True)
subprocess.run(['git', 'push', 'origin', 'refs/tags/' + tag], check=True)
assert git('ls-remote', 'origin', 'refs/tags/' + tag).split()[0] == commit
assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}') and not git('status', '--porcelain')
index.unlink()
value = {'tag': tag, 'commit': commit, 'source_revision': head, 'pushed': True, 'remote_commit_verified': True, 'files': records}
(WORK / 'closing-tools-preservation.json').write_text(json.dumps(value, indent=2) + '\n')
print(json.dumps({'tag': tag, 'commit': commit, 'files': len(records), 'pushed': True}))
