"""Preserve useful executable source and exact local witness bytes before recovery."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
SOURCE=json.loads((WORK/'source.json').read_text());HEAD=SOURCE['source_revision']
ARCHIVE=ROOT/'attestations/native-vulkan'/HEAD[:7]/'stage-seven-word-witness'
STAGE=ROOT/SOURCE['stage_root']
assert len(sys.argv)==2 and sys.argv[1] in ('source','archive','refresh')
def git(*args,**kwargs):return subprocess.check_output(['git',*args],cwd=ROOT,timeout=60,**kwargs)
def seal(p):
    data=p.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def save(p,r):
    assert not p.exists();p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(r,indent=2)+'\n')
def inventory(base):
    assert base.is_dir() and not base.is_symlink() and base.resolve()==base
    paths=[]
    for p in sorted(base.rglob('*')):
        assert not p.is_symlink(),p
        if p.is_file():paths.append(p)
        else:assert p.is_dir()
    return paths
assert ROOT.resolve()==ROOT and WORK==ROOT/'out/native-stage-seven-value-witness'
assert not git('status','--porcelain')
assert json.loads((WORK/'actual-results-review.json').read_text())['accepted']
for case in ('device-timestamps','mixed-transport','factored-endpoints'):
    for file,key in [('native-control-bridge.json','accepted'),('retirement-run.json','accepted'),('native-retirement.json','pass')]:
        assert json.loads((WORK/case/file).read_text(encoding='utf-8-sig'))[key]
if sys.argv[1]=='source':
    assert git('rev-parse','HEAD').decode().strip()==HEAD
    assert git('rev-parse','@{u}').decode().strip()==HEAD
    tag='evidence/issue-44-native-stage7-word-witness-'+HEAD[:7]
    candidates=[p for p in WORK.iterdir() if p.is_file() and p.suffix in ('.py','.ps1')]
    candidates += [p for d in ('kernels','versions') for p in (WORK/d).iterdir() if p.is_file() and p.suffix in ('.py','.slang')]
    candidates += [WORK/'actions.json',WORK/'value-schema.json']
    assert not git('tag','--list',tag).strip()
    index=WORK/'source-preservation.index';assert not index.exists()
    env=dict(os.environ,GIT_INDEX_FILE=str(index))
    try:
        git('read-tree',HEAD,env=env)
        records={}
        for p in sorted(candidates):
            name=str(p.relative_to(WORK));target=str(ARCHIVE.relative_to(ROOT)/'source'/name)
            data=p.read_bytes();blob=git('hash-object','-w','--stdin',input=data).decode().strip()
            git('update-index','--add','--cacheinfo','100644,'+blob+','+target,env=env)
            records[name]={'git_path':target,**seal(p)}
        tree=git('write-tree',env=env).decode().strip()
        message=('Establish finite native stage-seven/Endpoint word witness (#44)\n\n'
                 'All three original bounded native controls pass. The162 curved high pairs have exact normalized keys and5022 complete five-word results, with original32 digests/51 ordered markers. Preserve executable observer/reader/ownership source. Invalid intermediates remain distinct from missing records and refusal;68 observed invalid records occur in successful standalone rows. No universal word theorem, safe reuse, complete geometry owner, ordinary benefit, frame or full qualification is established.\n')
        commit=git('commit-tree',tree,'-p',HEAD,input=message.encode()).decode().strip()
        git('update-ref','refs/tags/'+tag,commit)
        save(WORK/'source-preservation-local.json',{'source_revision':HEAD,'source_tag':tag,'source_commit':commit,'source_versions':records,'accepted':True,'pushed':False})
        git('push','origin','refs/tags/'+tag)
        assert git('ls-remote','origin','refs/tags/'+tag).split()[0].decode()==commit
        for name,r in records.items():assert git('show',commit+':'+r['git_path'])==(WORK/name).read_bytes()
        save(WORK/'source-preservation-pushed.json',{'source_revision':HEAD,'source_tag':tag,'source_commit':commit,'source_versions':records,'accepted':True,'pushed':True})
        print(json.dumps({'source_tag':tag,'source_commit':commit,'useful_source_files':len(records)}))
    finally:
        index.unlink(missing_ok=True)
elif sys.argv[1]=='archive':
    assert git('rev-parse','HEAD').decode().strip()==HEAD and not ARCHIVE.exists()
    assert json.loads((WORK/'source-preservation-pushed.json').read_text())['pushed']
    products=WORK/'diagnostic-products';assert not products.exists();products.mkdir()
    for name in ('bin/linux-gcc/src/sirius/app/sirius','bin/linux-gcc/tests/backend/sirius_backend_tests','bin/linux-gcc/src/sirius/app/sirius_app_tests','bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h','bin/linux-gcc/generated/sirius/alignment_receipt.json','bin/linux-gcc/src/sirius/app/resources/model/alignment_receipt.json','bin/linux-gcc/CMakeCache.txt'):
        p=ROOT/name;q=products/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);assert seal(p)==seal(q)
    records=[]
    for base,part in [(WORK,'files'),(STAGE,'native-stage')]:
        for p in inventory(base):
            expected=seal(p);q=ARCHIVE/part/p.relative_to(base);q.parent.mkdir(parents=True,exist_ok=True)
            assert not q.exists();shutil.copy2(p,q);assert seal(p)==expected==seal(q)
            records.append({'origin':str(p.relative_to(ROOT)),'recovery':str(q.relative_to(ROOT)),**expected})
    save(ARCHIVE/'initial-inventory.json',{'accepted':True,'source_revision':HEAD,'records':records,'scope':'Exact finished witness WORK, source/configuration, official export, complete native stage and actual diagnostic Linux products. Raw binaries/evidence remain local; pushed tag preserves useful executable source.'})
    print(json.dumps({'archive':str(ARCHIVE.relative_to(ROOT)),'files':len(records),'bytes':sum(r['bytes'] for r in records)}))
else:
    assert (ARCHIVE/'initial-inventory.json').exists()
    records=[]
    for base,part in [(WORK,'files'),(STAGE,'native-stage')]:
        for p in inventory(base):
            q=ARCHIVE/part/p.relative_to(base);expected=seal(p)
            if q.exists():assert seal(q)==expected,(p,q)
            else:q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
            assert seal(q)==expected;records.append({'origin':str(p.relative_to(ROOT)),'recovery':str(q.relative_to(ROOT)),**expected})
    n=len(list(ARCHIVE.glob('refresh-*.json')))
    save(ARCHIVE/('refresh-'+str(n)+'.json'),{'accepted':True,'source_revision':HEAD,'records':records})
    print(json.dumps({'refreshed_files':len(records),'bytes':sum(r['bytes'] for r in records)}))
