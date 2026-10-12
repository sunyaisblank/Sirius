"""Preserve machine-screen executable inputs without changing production HEAD."""
from pathlib import Path
import hashlib
import json
import os
import stat
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
revision = 'a05763a68fa371e927eda1b7f83cc97c98465323'
tag = 'evidence/issue-44-portable-indexed-machine-screen-a05763a'
prefix = 'attestations/source/issue-44-portable-indexed-machine-screen/'
index = work / 'machine-source-preservation.index'
assert not index.exists()
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
env = os.environ.copy()
env['GIT_INDEX_FILE'] = str(index)
def git(*args, data=None):
    return subprocess.check_output(['git', *args], input=data, cwd=root, env=env).decode().strip()
names = ['emit_backend.py', 'elf_layout.py', 'verify_machine.py', 'emission-settings.json',
         'actions-before-machine-emission.json', 'actions-before-machine-verification.json',
         'actions.json', 'preserve_machine_sources.py']
records = []
try:
    git('read-tree', revision)
    for name in names:
        source = work / name
        data = source.read_bytes()
        blob = git('hash-object', '-w', '--stdin', data=data)
        mode = '100755' if source.stat().st_mode & stat.S_IXUSR else '100644'
        alias = prefix + name
        git('update-index', '--add', '--cacheinfo', mode + ',' + blob + ',' + alias)
        records.append({'source': str(source.relative_to(root)), 'alias': alias,
                        'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                        'blob': blob, 'mode': mode})
    commit = git('commit-tree', git('write-tree'), '-p', revision,
                 data=b'Preserve indexed limb machine emission and complete JIT correspondence checks\n')
    git('tag', '-a', tag, commit, '-m', 'Immutable issue 44 indexed limb machine screen sources')
    subprocess.run(['git', 'push', 'origin', 'refs/tags/' + tag], check=True)
    remote = git('ls-remote', 'origin', 'refs/tags/' + tag, 'refs/tags/' + tag + '^{}')
    assert commit + '\trefs/tags/' + tag + '^{}' in remote
    for item in records:
        raw = subprocess.check_output(['git', 'show', commit + ':' + item['alias']])
        assert hashlib.sha256(raw).hexdigest() == item['sha256']
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
    assert not subprocess.check_output(['git', 'status', '--porcelain'])
    result = {'passed': True, 'tag': tag, 'commit': commit, 'parent_revision': revision,
              'aliases': records, 'source_aliases': len(records),
              'distinct_blobs': len({x['blob'] for x in records}),
              'remote_readback_exact': True, 'production_HEAD_unchanged': True}
    (work / 'machine-source-git-bindings.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'aliases'}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name + '.lock').unlink(missing_ok=True)
