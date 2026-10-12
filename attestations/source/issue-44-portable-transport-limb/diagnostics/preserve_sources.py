"""Preserve executable diagnostics and frozen specifications without changing HEAD."""
from pathlib import Path
import hashlib, json, os, stat, subprocess
r=Path.cwd();w=Path(__file__).resolve().parent
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert revision=='d910f3c3b186bb5ed28a8102faf84b24698d679d'
assert subprocess.check_output(['git','status','--porcelain'],text=True)==''
namespace='attestations/source/issue-44-portable-transport-limb'
tag='evidence/issue-44-portable-transport-limb-d910f3c'
index=w/'preservation.index';assert not index.exists()
env=os.environ.copy();env['GIT_INDEX_FILE']=str(index)
def git(*args,input=None):
    return subprocess.check_output(['git',*args],input=input,cwd=r,env=env).decode().strip()
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag]).returncode==1
sources=[]
for name in ['run_owned.py','owned_guard.py','build_transport_products.py','verify_transport_products.py',
             'run_consumer_compile.py','verify_consumers.py','run_control.py','verify_control.py',
             'prepare_benefit.py','run_limb_benefit_slot.py','compare_benefit.py','preserve_sources.py',
             'actions-before-products.json','actions-before-product-verification.json','actions-before-consumers.json',
             'actions-before-consumer-verification.json','actions-before-integrated-control.json',
             'actions-before-integrated-verification.json','actions-before-benefit-preparation.json',
             'actions-before-benefit-comparison.json','actions.json']:
    sources.append((w/name,'diagnostics/'+name))
for folder in ['original-controls/candidate','control-contracts','benefit']:
    for source in sorted((w/folder).iterdir()):
        if source.is_file() and (source.suffix=='.py' or source.name in ['actions.json','retention-gate-frozen.json','expected-readback-layout.json','baseline.json','build-and-arrays.json']):
            sources.append((source,folder+'/'+source.name))
for source in [r/'src/sirius/kernels/retained_triple.slang',r/'scripts/build-retained-kernels.py']:
    sources.append((source,str(source.relative_to(r))))
assert all(source.is_file() for source,_ in sources)
records=[]
try:
    git('read-tree',revision)
    for source,alias in sources:
        data=source.read_bytes();blob=git('hash-object','-w','--stdin',input=data)
        mode='100755' if source.stat().st_mode&stat.S_IXUSR else '100644'
        path=namespace+'/'+alias
        git('update-index','--add','--cacheinfo',mode+','+blob+','+path)
        records.append({'source':str(source.relative_to(r)),'alias':path,'bytes':len(data),
                        'sha256':hashlib.sha256(data).hexdigest(),'blob':blob,'mode':mode})
    tree=git('write-tree')
    commit=git('commit-tree',tree,'-p',revision,input=b'Preserve Transport-only limb diagnostics and frozen original comparison\n')
    git('tag','-a',tag,commit,'-m','Immutable issue 44 Transport-only limb sources and frozen original gates')
    subprocess.run(['git','push','origin','refs/tags/'+tag],cwd=r,check=True)
    remote=git('ls-remote','origin','refs/tags/'+tag,'refs/tags/'+tag+'^{}')
    assert commit+'\trefs/tags/'+tag+'^{}' in remote
    for record in records:
        data=subprocess.check_output(['git','show',commit+':'+record['alias']],cwd=r)
        assert hashlib.sha256(data).hexdigest()==record['sha256']
    result={'passed':True,'tag':tag,'peeled_commit':commit,'parent_revision':revision,
            'aliases':len(records),'distinct_blobs':len({x['blob'] for x in records}),'records':records,
            'remote_readback_exact':True,'root_HEAD_unchanged':git('rev-parse','HEAD')==revision}
    assert result['root_HEAD_unchanged']
    (w/'source-git-bindings.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key!='records'}))
finally:
    if index.exists():index.unlink()
