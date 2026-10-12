"""Preserve useful arithmetic/consumer witnesses without moving production HEAD."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root=Path.cwd()
work=Path(__file__).resolve().parent
revision='c8199bd25e4851b2f139f047b913c82352bace41'
tag='evidence/issue-44-portable-limb-closure-c8199bd'
prefix='attestations/source/issue-44-portable-limb-closure/'
def git(*args,data=None,env=None):
    return subprocess.check_output(['git',*args],cwd=root,input=data,env=env,timeout=60)
assert git('rev-parse','HEAD').decode().strip()==revision
assert not git('status','--porcelain')
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag],cwd=root).returncode==1
actions=json.loads((work/'actions.json').read_text())
assert len(actions)==51
for key in actions:
    owner=json.loads((work/key/'owner.json').read_text())
    assert owner['status']=='terminal' and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    if key in ('compile-limb-screen','compile-short-input-probe'):
        assert owner['child_exit']==1 and not owner['passed'] and owner['stop_reason'] is None
    else:
        assert owner['child_exit']==0 and owner['passed']
sources=[(p,p.name) for p in sorted(work.glob('*.py'))]
sources += [(p,p.name) for suffix in ('*.cpp','*.h','*.slang') for p in sorted(work.glob(suffix))]
sources += [(p,p.name) for p in sorted(work.glob('actions*.json'))]
sources += [(work/name,name) for name in ('emission-settings.json','limb-screen-source.json',
    'short-input-preset-recipe.json','limb-arithmetic-baseline-source.json','baseline-reuse.json',
    'limb-screen-source.patch')]
sources += [(p,'original-controls/candidate/'+p.name) for p in sorted((work/'original-controls/candidate').glob('*.cpp'))]
sources += [(p,'control-contracts/'+p.name) for p in sorted((work/'control-contracts').iterdir()) if p.is_file()]
sources += [(work/'benefit'/name,'benefit/'+name) for name in ('software_contracts.py','compare_software.py','run_sequence.py','run_owned.py','owned_guard.py','run_limb_benefit_slot.py','expected-readback-layout.json','actions.json','retention-gate-frozen.json','edits.json','baseline.json','build-and-arrays.json')]
for label in ('a1','b1','b2','a2'):
    owner=json.loads((work/'benefit'/label/'owner.json').read_text())
    assert owner['passed'] and owner['status']=='terminal' and owner['remaining_owned_processes']==[]
assert json.loads((work/'benefit/software-comparison.json').read_text())['pass'] is False
baseline=root/'bin/linux-gcc/portable-execution-attribution/limb-screen/baseline-kernels'
sources += [(p,'baseline-kernels/'+p.name) for p in sorted(baseline.iterdir()) if p.is_file()]
old=root/'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation'
sources += [(old/'generate_triples.py','original-generate-triples.py')]
assert len({name for _,name in sources})==len(sources)
index=work/'limb-closure-source-preservation.index'
receipt=work/'limb-closure-source-git-bindings.json'
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
    commit=git('commit-tree',tree,'-p',revision,data=b'Preserve retired portable limb numerical and original performance witnesses\n\nFour-mode exact Fraction numerical/refusal/recovery and complete-word comparison; complete actual application-consumer payload joins; original integrated control8. Sole original ABBA rejects adoption with18/20 numerical timing and0/2 API gates; source restored exactly. No retries, new CI, frame or full scientific/release qualification claim. Preserve both earlier failed preparation versions.\n').decode().strip()
    for binding in bindings:
        content=git('show',commit+':'+binding['git_path'])
        assert len(content)==binding['bytes'] and hashlib.sha256(content).hexdigest()==binding['sha256']
    git('tag','-a',tag,commit,'-m','Portable limb retirement and exact original performance evidence sources')
    receipt.write_text(json.dumps({'tag':tag,'commit':commit,'parent':revision,'bindings':bindings,
        'production_HEAD_unchanged':True,'push_pending':True},indent=2)+'\n')
    print(json.dumps({'tag':tag,'commit':commit,'source_aliases':len(bindings)}))
finally:
    index.unlink(missing_ok=True)
    index.with_name(index.name+'.lock').unlink(missing_ok=True)
