"""Preserve the frozen original performance specification before timed execution."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root=Path.cwd()
work=Path(__file__).resolve().parent
benefit=work/'benefit'
revision='661e157b53f9f9fbbba0347cfda32863a1684c75'
tag='evidence/issue-44-portable-limb-benefit-freeze-v2-661e157'
prefix='attestations/source/issue-44-portable-limb-benefit-freeze-v2/'
def git(*args,data=None,env=None):
    return subprocess.check_output(['git',*args],cwd=root,input=data,env=env,timeout=60)
assert git('rev-parse','HEAD').decode().strip()==revision and not git('status','--porcelain')
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag],cwd=root).returncode==1
assert not any((benefit/key).exists() for key in ('a1','b1','b2','a2'))
owner=json.loads((work/'prepare-limb-benefit-v2/owner.json').read_text())
assert owner['passed'] and owner['status']=='terminal' and owner['whole_inputs_unchanged']
sources=[(work/name,name) for name in ('prepare_limb_benefit_v2.py','run_limb_benefit_slot.py',
    'preserve_limb_benefit_sources_v2.py','actions.json')]
sources += [(p,'benefit/'+p.name) for p in sorted(benefit.iterdir()) if p.is_file()]
index=work/'limb-benefit-v2-source-preservation.index'
receipt=work/'limb-benefit-v2-source-git-bindings.json'
assert not index.exists() and not receipt.exists()
env=os.environ.copy();env['GIT_INDEX_FILE']=str(index)
bindings=[]
try:
    git('read-tree',revision,env=env)
    for source,name in sources:
        content=source.read_bytes()
        blob=git('hash-object','-w','--stdin',data=content).decode().strip()
        target=prefix+name
        git('update-index','--add','--cacheinfo','100644,'+blob+','+target,env=env)
        bindings.append({'source':str(source.relative_to(root)),'git_path':target,
            'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'blob':blob})
    tree=git('write-tree',env=env).decode().strip()
    commit=git('commit-tree',tree,'-p',revision,data=b'Preserve corrected original portable limb ABBA specification before execution\n\nUnchanged twenty numerical performance gates plus two inclusive API gates, original slots and scientific/source/stream outputs. Exact sealed wrapper/action ownership binding corrects the preserved prelaunch-only metadata mismatch. No timed slot launched and no benefit acceptance claimed.\n').decode().strip()
    for binding in bindings:
        content=git('show',commit+':'+binding['git_path'])
        assert len(content)==binding['bytes'] and hashlib.sha256(content).hexdigest()==binding['sha256']
    git('tag','-a',tag,commit,'-m','Portable limb original benefit gates frozen before execution')
    receipt.write_text(json.dumps({'tag':tag,'commit':commit,'parent':revision,'bindings':bindings,
        'production_HEAD_unchanged':True,'push_pending':True},indent=2)+'\n')
    print(json.dumps({'tag':tag,'commit':commit,'source_aliases':len(bindings)}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name+'.lock').unlink(missing_ok=True)
