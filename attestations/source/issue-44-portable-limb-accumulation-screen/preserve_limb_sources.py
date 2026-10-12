"""Preserve the source-equivalent limb compiler screen without advancing production."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root=Path.cwd()
work=Path(__file__).resolve().parent
revision='95c4293f598c1dbf4cd13ece7b06ffff5383b81a'
tag='evidence/issue-44-portable-limb-accumulation-screen-95c4293'
prefix='attestations/source/issue-44-portable-limb-accumulation-screen/'
def git(*args,data=None,env=None):
    return subprocess.check_output(['git',*args],cwd=root,input=data,env=env,timeout=60)
assert git('rev-parse','HEAD').decode().strip()==revision
assert not git('status','--porcelain')
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag],cwd=root).returncode==1
owners=list(work.glob('*/owner.json'))
assert len(owners)==11
for path in owners:
    owner=json.loads(path.read_text())
    assert owner['status']=='terminal' and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    if path.parent.name in ('compile-limb-screen','compile-short-input-probe'):
        assert owner['child_exit']==1 and not owner['passed'] and owner['stop_reason'] is None
    else:
        assert owner['child_exit']==0 and owner['passed']
sources=[(p,p.name) for p in sorted(work.glob('*.py'))]
sources += [(work/'short_input_jit_probe.cpp','short_input_jit_probe.cpp'),
            (work/'limb-screen-source.patch','limb-screen-source.patch')]
sources += [(p,p.name) for p in sorted(work.glob('actions*.json'))]
sources += [(work/name,name) for name in ('emission-settings.json','limb-screen-source.json','short-input-preset-recipe.json')]
stage=root/'bin/linux-gcc/portable-execution-attribution/limb-screen/kernels'
sources += [(p,'stage-src/kernels/'+p.name) for p in sorted(stage.iterdir()) if p.is_file()]
old=root/'attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding'
sources += [(old/'optimizer_observer.c','optimizer_observer.c'),(old/'observer-build.json','observer-build.json')]
index=work/'limb-source-preservation.index'
receipt=work/'limb-source-git-bindings.json'
assert not index.exists() and not receipt.exists()
env=os.environ.copy(); env['GIT_INDEX_FILE']=str(index)
bindings=[]
try:
    git('read-tree',revision,env=env)
    for source,name in sources:
        content=source.read_bytes()
        blob=git('hash-object','-w','--stdin',data=content).decode().strip()
        target=prefix+name
        git('update-index','--add','--cacheinfo','100644,'+blob+','+target,env=env)
        bindings.append({'source':str(source.relative_to(root)),'git_path':target,'bytes':len(content),
                         'sha256':hashlib.sha256(content).hexdigest(),'blob':blob})
    tree=git('write-tree',env=env).decode().strip()
    commit=git('commit-tree',tree,'-p',revision,data=b'Preserve reviewed portable limb accumulation compiler-screen sources\n\nSource and compiled carry/guard equivalence, source-bound machine counter removal and bounded compiler footprint. Numerical execution and whole-path benefit remain pending. Preserve both failed preparation versions.\n').decode().strip()
    for binding in bindings:
        content=git('show',commit+':'+binding['git_path'])
        assert len(content)==binding['bytes'] and hashlib.sha256(content).hexdigest()==binding['sha256']
    git('tag','-a',tag,commit,'-m','Portable limb accumulation source and compiler-feasibility evidence')
    receipt.write_text(json.dumps({'tag':tag,'commit':commit,'parent':revision,'bindings':bindings,
                                   'production_HEAD_unchanged':True,'push_pending':True},indent=2)+'\n')
    print(json.dumps({'tag':tag,'commit':commit,'source_aliases':len(bindings)}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name+'.lock').unlink(missing_ok=True)
