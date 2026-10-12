"""Pushable cleanup source binding; keep production HEAD unchanged."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root=Path.cwd()
archive=Path(__file__).resolve().parent
revision='5b11d31b813b47a55117f5bccff5119683e9a7a0'
tag='evidence/issue-44-portable-indexed-limb-cleanup-5b11d31'
prefix='attestations/source/issue-44-portable-indexed-limb-cleanup/'
def git(*args,data=None,env=None):
    return subprocess.check_output(['git',*args],cwd=root,input=data,env=env,timeout=60)
assert git('rev-parse','HEAD').decode().strip()==revision and not git('status','--porcelain')
assert not git('diff','c8199bd','HEAD')
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag],cwd=root).returncode==1
assert json.loads((archive/'recorded-input-recovery.json').read_text())['passed']
raw=archive/'files/out/portable-indexed-limb-accumulation'
sources=[(archive/'scoped_cleanup.py','scoped_cleanup.py'),
    (archive/'preserve_cleanup_source.py','preserve_cleanup_source.py'),
    (raw/'archive_batch.py','archive_batch.py'),(raw/'recover_recorded_inputs.py','recover_recorded_inputs.py')]
index=archive/'cleanup-source-preservation.index'
receipt=archive/'cleanup-source-git-bindings.json'
assert not index.exists() and not receipt.exists()
env=os.environ.copy();env['GIT_INDEX_FILE']=str(index)
aliases=[]
try:
    git('read-tree',revision,env=env)
    for path,name in sources:
        content=path.read_bytes()
        blob=git('hash-object','-w','--stdin',data=content).decode().strip()
        target=prefix+name
        git('update-index','--add','--cacheinfo','100644,'+blob+','+target,env=env)
        aliases.append({'source':str(path.relative_to(root)),'alias':target,'blob':blob,
            'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()})
    tree=git('write-tree',env=env).decode().strip()
    commit=git('commit-tree',tree,'-p',revision,data=b'Preserve exact portable limb raw recovery and scoped cleanup source\n\nRetired trial only. Require exact retired-batch file recovery, every consumed input version, protected artifacts, synced source and fresh scoped consumer census before removal. No scientific/release qualification claim.\n').decode().strip()
    for binding in aliases:
        data=git('show',commit+':'+binding['alias'])
        assert len(data)==binding['bytes'] and hashlib.sha256(data).hexdigest()==binding['sha256']
    git('tag','-a',tag,commit,'-m','Portable limb retired-trial recovery and guarded cleanup source')
    receipt.write_text(json.dumps({'tag':tag,'commit':commit,'parent':revision,'aliases':aliases,
        'production_HEAD_unchanged':True,'push_pending':True},indent=2)+'\n')
    print(json.dumps({'tag':tag,'commit':commit,'source_aliases':len(aliases)}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name+'.lock').unlink(missing_ok=True)
