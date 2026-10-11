"""Push source-only diagnostic tools and concise finite findings without altering checkout."""
from pathlib import Path
import hashlib, json, os, subprocess

root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert head == 'b81e061fa081cdaed37eca3add950bd254bf296d'
assert not subprocess.check_output(['git', 'status', '--porcelain'])
tag = 'evidence/issue-94-rejected-upmultiply-cleanup-b81e061-v2'
parent = '8a1377192d00c5534f40288e08c9970b6b536f5c'
names = ['bind_restoration.py', 'protect_restoration.py', 'verify_restoration.py',
         'finish_cleanup.py', 'preserve_tools.py', 'preserve_tools_v3.py', 'cleanup-first-attempt.json', 'restoration/actions.json',
         'restoration/run_owned.py', 'restoration/owned_guard.py', 'restoration/source.json',
         'restoration-source-contract.json', 'restoration-source-independent-review.json',
         'restoration-collector-source-review.json', 'restoration-preservation-reviewed.json',
         'restoration-protection.json', 'restoration-source.json', 'cleanup-source-review.json',
         'trial-independent-disposition.json', 'trial-preservation-reviewed.json',
         'trial-protection.json', 'tools-preservation.json']
index = work / 'source-only-v3.index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
def git(*argv, input=None):
    return subprocess.check_output(['git', *argv], input=input, env=env)
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
files = {}
try:
    git('read-tree', parent)
    for name in names:
        path = work / name
        assert path.is_file() and not path.is_symlink()
        data = path.read_bytes()
        oid = git('hash-object', '-w', '--stdin', input=data).decode().strip()
        git('update-index', '--add', '--cacheinfo', '100644', oid, name)
        files[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    tree = git('write-tree').decode().strip()
    commit = git('commit-tree', tree, '-p', parent, input=b'Preserve rejected #94 restoration and cleanup sources\n\nSource-only evidence snapshot; never apply as a Sirius product patch.\n').decode().strip()
    subprocess.check_call(['git', 'tag', '-a', tag, commit, '-m', 'Source-only rejected #94 restoration and cleanup tools; not a product patch'])
    subprocess.check_call(['git', 'push', 'origin', 'refs/tags/' + tag])
    assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + tag + '^{}'], text=True).split()[0] == commit
    for name, value in files.items():
        data = subprocess.check_output(['git', 'show', tag + ':' + name])
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == value
    assert not subprocess.check_output(['git', 'status', '--porcelain'])
    receipt = {'source_revision': head, 'tag': tag, 'commit': commit, 'parent_source_only_snapshot': parent,
               'pushed': True, 'files': files, 'added_or_refreshed_files': len(files),
               'scope': 'Exact reviewed diagnostic/cleanup sources and concise finite findings. Parent retains trial tools. No product patch, raw binary backup or full qualification claim.'}
    (work / 'tools-preservation-v3.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'files'}))
finally:
    index.unlink(missing_ok=True)
