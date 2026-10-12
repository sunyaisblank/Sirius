"""Preserve this staged compiler screen without changing production HEAD."""
from pathlib import Path
import hashlib,json,os,stat,subprocess
r=Path.cwd();w=Path(__file__).resolve().parent;b=r/'bin/linux-gcc/portable-indexed-limb-accumulation'
revision='a05763a68fa371e927eda1b7f83cc97c98465323'
tag='evidence/issue-44-portable-indexed-limb-screen-a05763a'
prefix='attestations/source/issue-44-portable-indexed-limb-screen/'
index=w/'source-preservation.index';assert not index.exists()
env=os.environ.copy();env['GIT_INDEX_FILE']=str(index)
def git(*args,data=None):return subprocess.check_output(['git',*args],input=data,cwd=r,env=env).decode().strip()
assert git('rev-parse','HEAD')==revision and git('status','--porcelain')==''
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag]).returncode==1
sources=[(p,'kernels/'+p.name) for p in sorted((b/'kernels').iterdir()) if p.is_file()]
for name in ['run_owned.py','owned_guard.py','compile_screen.py','verify_screen.py','run_capture.py',
             'capture_short_input_gdb.py','short_input_jit_probe.cpp','short-input-preset-recipe.json',
             'candidate-source.patch','candidate-source.json','actions.json','actions-before-compiler.json',
             'actions-before-SPIRV-verification.json','actions-before-capture.json','preserve_sources.py']:
    sources.append((w/name,'diagnostics/'+name))
sources.append((r/'scripts/build-retained-kernels.py','scripts/build-retained-kernels.py'))
records=[]
try:
    git('read-tree',revision)
    for source,alias in sources:
        data=source.read_bytes();blob=git('hash-object','-w','--stdin',data=data)
        mode='100755' if source.stat().st_mode&stat.S_IXUSR else '100644';path=prefix+alias
        git('update-index','--add','--cacheinfo',mode+','+blob+','+path)
        records.append({'source':str(source.relative_to(r)),'alias':path,'bytes':len(data),
                        'sha256':hashlib.sha256(data).hexdigest(),'blob':blob,'mode':mode})
    commit=git('commit-tree',git('write-tree'),'-p',revision,data=b'Preserve private-array limb compiler screen and fixed observer inputs\n')
    git('tag','-a',tag,commit,'-m','Immutable issue 44 private-array compiler screen sources')
    subprocess.run(['git','push','origin','refs/tags/'+tag],check=True)
    remote=git('ls-remote','origin','refs/tags/'+tag,'refs/tags/'+tag+'^{}')
    assert commit+'\trefs/tags/'+tag+'^{}' in remote
    for item in records:
        raw=subprocess.check_output(['git','show',commit+':'+item['alias']])
        assert hashlib.sha256(raw).hexdigest()==item['sha256']
    assert git('rev-parse','HEAD')==revision and git('status','--porcelain')==''
    result={'passed':True,'tag':tag,'commit':commit,'parent_revision':revision,
            'aliases':records,'source_aliases':len(records),'distinct_blobs':len({x['blob'] for x in records}),
            'remote_readback_exact':True,'production_HEAD_unchanged':True}
    (w/'source-git-bindings.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key!='aliases'}))
finally:
    index.unlink(missing_ok=True);index.with_name(index.name+'.lock').unlink(missing_ok=True)
