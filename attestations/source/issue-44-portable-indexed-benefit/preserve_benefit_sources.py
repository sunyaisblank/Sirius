"""Preserve frozen original performance executable inputs without changing production HEAD."""
from pathlib import Path
import hashlib
import json
import os
import stat
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
revision = '6014f91d28ddeadbc7b916f213be7f89a4d1f16c'
tag = 'evidence/issue-44-portable-indexed-benefit-6014f91'
prefix = 'attestations/source/issue-44-portable-indexed-benefit/'
index = work / 'benefit-source-preservation.index'
assert not index.exists()
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
env = os.environ.copy()
env['GIT_INDEX_FILE'] = str(index)
def git(*args, data=None):
    return subprocess.check_output(['git', *args], input=data, cwd=root, env=env).decode().strip()
names = ['prepare_benefit.py', 'compare_benefit.py', 'run_limb_benefit_slot.py', 'run_limb_benefit_slot_v1.py', 'benefit/software_contracts.py', 'benefit/compare_software.py', 'benefit/run_sequence.py', 'benefit/expected-readback-layout.json', 'benefit/run_owned.py', 'benefit/owned_guard.py', 'benefit/run_limb_benefit_slot.py', 'benefit/baseline.json', 'benefit/build-and-arrays.json', 'benefit/actions.json', 'benefit/retention-gate-frozen.json', 'benefit/edits.json', 'actions-before-benefit-preparation.json', 'actions.json', 'preserve_benefit_sources.py']
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
                 data=b'Preserve indexed limb original frozen ABBA source and all20plus2 gates\n')
    git('tag', '-a', tag, commit, '-m', 'Immutable issue 44 indexed limb original performance source and frozen gates')
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
    (work / 'benefit-source-git-bindings.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'aliases'}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name + '.lock').unlink(missing_ok=True)
