"""Preserve executable diagnostic sources under an evidence tag without changing HEAD."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
revision = '95c4293f598c1dbf4cd13ece7b06ffff5383b81a'
tag = 'evidence/issue-44-portable-execution-attribution-95c4293'
prefix = 'attestations/source/issue-44-portable-execution-attribution/'
def git(*args, data=None, env=None):
    return subprocess.check_output(['git', *args], input=data, cwd=root, env=env, timeout=60)
assert git('rev-parse', 'HEAD').decode().strip() == revision
assert not git('status', '--porcelain')
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag],cwd=root).returncode == 1
owners = list(work.glob('*/owner.json'))
assert len(owners) == 5
for path in owners:
    value = json.loads(path.read_text())
    assert value['status'] == 'terminal' and value['passed'] and not value['remaining_owned_processes']
names = ['emit_backend.py', 'assemble_labels.py', 'assemble_llvm_labels.py',
         'assemble_machine_blocks.py', 'verify_emission.py', 'elf_layout.py',
         'verify_machine_blocks.py', 'owned_guard.py', 'run_owned.py', 'preserve_sources.py',
         'actions-before-labels.json', 'actions-before-LLVM-labels.json',
         'actions-before-machine-blocks.json', 'actions-before-machine-verification.json',
         'actions.json', 'emission-settings.json']
index = work / 'source-preservation.index'
receipt = work / 'source-git-bindings.json'
assert not index.exists() and not receipt.exists()
env = os.environ.copy()
env['GIT_INDEX_FILE'] = str(index)
bindings = []
try:
    git('read-tree', revision, env=env)
    for name in names:
        source = work / name
        content = source.read_bytes()
        blob = git('hash-object', '-w', '--stdin', data=content).decode().strip()
        target = prefix + name
        git('update-index', '--add', '--cacheinfo', '100644,'+blob+','+target, env=env)
        bindings.append({'source': str(source.relative_to(root)), 'git_path': target,
                         'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest(), 'blob': blob})
    tree = git('write-tree', env=env).decode().strip()
    message = ('Preserve exact offline Endpoint machine-block attribution sources\n\n'
               'Offline code/constant/unwind and symbolic relocation equality gates; '
               'complete emitted block and historical sample joins. No production '
               'implementation, performance gain or release qualification claim.\n')
    commit = git('commit-tree', tree, '-p', revision, data=message.encode()).decode().strip()
    for binding in bindings:
        content = git('show', commit+':'+binding['git_path'])
        assert len(content)==binding['bytes'] and hashlib.sha256(content).hexdigest()==binding['sha256']
    git('tag', '-a', tag, commit, '-m', 'Exact offline Endpoint block attribution executable evidence')
    result={'tag':tag,'commit':commit,'parent':revision,'bindings':bindings,
            'production_HEAD_unchanged':git('rev-parse','HEAD').decode().strip()==revision,
            'push_pending':True}
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'tag':tag,'commit':commit,'source_files':len(bindings)}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name+'.lock').unlink(missing_ok=True)
