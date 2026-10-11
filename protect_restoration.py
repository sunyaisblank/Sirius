"""Preserve restored product/input epochs; retain the finished protection owner."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def expected(v):return {k:v[k] for k in ('bytes','sha256')}
def key(v):return v['bytes'],v['sha256']
source=read(WORK/'restoration/source.json');head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{u}'])
review=read(WORK/'restoration-independent-review.json');assert review['accepted'] and review['source_revision']==head
assert subprocess.check_output(['git','rev-parse','HEAD^{tree}'])==subprocess.check_output(['git','rev-parse',source['baseline_revision']+'^{tree}'])
for case in ('configure','build','products','bindings','admission','selection'):
    o=read(WORK/'restoration'/case/'owner.json')
    assert o['passed'] and o['source']=={'revision':head,'status':''} and not o['stop_reason'] and not o['cleanup_errors'] and o['remaining_owned_processes']==[]
trial=read(WORK/'trial-protection.json');prior=ROOT/trial['archive']
assert seal(prior/'manifest.json')==trial['manifest']
old=read(prior/'manifest.json');recoveries={}
for v in [*old['executed_input_ledger'],*old['diagnostics'],*old['volatile_product_and_native_stage_files']]:
    recoveries[key(v)]=v['recovery']
DEST=ROOT/'attestations/review'/head[:7]/'native-narrow-fma-restoration'
assert not DEST.exists()
payloads={}
def resolve(p,v):
    k=key(v)
    if k not in recoveries:
        assert p.is_file() and seal(p)==expected(v),str(p)
        name='payloads/'+v['sha256'];q=DEST/name;q.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,q);assert seal(q)==expected(v)
        payloads[name]={'path':name,**expected(v)}
        recoveries[k]={'kind':'local-archive-payload','path':str(q.relative_to(ROOT))}
    r=recoveries[k]
    if r['kind']=='immutable-Git-blob':
        b=subprocess.check_output(['git','cat-file','blob',r['blob']]);assert {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}==expected(v)
    else:assert seal(ROOT/r['path'])==expected(v),r
    return r

ledger=[];seen=set()
for folder in (WORK/'restoration',WORK/'protection'):
    for p in [*folder.glob('*/inputs-before.json'),*folder.glob('*/inputs-after.json'),*folder.glob('inputs-before.json'),*folder.glob('inputs-after.json')]:
        for origin,v in read(p).items():
            if (origin,*key(v)) in seen:continue
            seen.add((origin,*key(v)));q=ROOT/origin
            current=q.is_file() and seal(q)==expected(v)
            if folder.name=='restoration':assert current,origin
            ledger.append({'origin':origin,**expected(v),'ledger':str(p.relative_to(ROOT)),
                'origin_current_at_preservation':current,'recovery':resolve(q,v)})
records=[]
paths=[WORK/'restoration',WORK/'protection',WORK/'trial-protection.json',
    WORK/'trial-preservation-reviewed.json',WORK/'restoration-independent-review.json',
    WORK/'protect_restoration.py',WORK/'run_restoration_protection_owned.py',
    WORK/'restoration-protection-actions.json',WORK/'restoration-protection-source-review.json',
    ROOT/'bin/linux-gcc/src/sirius/backend/retained',ROOT/'bin/linux-gcc/generated',
    ROOT/'bin/linux-gcc/CMakeCache.txt',ROOT/'bin/linux-gcc/build.ninja']
for directory in paths:
    assert directory.exists() and not directory.is_symlink()
    for p in sorted(directory.rglob('*')) if directory.is_dir() else [directory]:
        assert not p.is_symlink()
        if p.is_file():
            v=seal(p);records.append({'origin':str(p.relative_to(ROOT)),**v,'recovery':resolve(p,v)})
manifest={'source_revision':head,'baseline_revision':source['baseline_revision'],
    'trial_manifest':{'path':str((prior/'manifest.json').relative_to(ROOT)),**trial['manifest']},
    'review':review,'executed_input_ledger':ledger,'records':records,'payloads':list(payloads.values()),
    'payload_count':len(payloads),'payload_bytes':sum(v['bytes'] for v in payloads.values()),
    'scope':'Exact restored source/build/product and host-only 4+1 admission evidence. Initial stale filter and missing-case correction remain distinct. Completed primary protection owner preserved; current restoration-protection owner excluded until closing publication. Original two obsolete SPVs are retained pending final cleanup. No native/frame/full-science/release or cleanup authorization.'}
DEST.mkdir(parents=True,exist_ok=True);(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'archive':str(DEST.relative_to(ROOT)),'manifest':seal(DEST/'manifest.json'),'payload_count':len(payloads),
    'payload_bytes':manifest['payload_bytes'],'executed_origin_version_pairs':len(ledger),'records':len(records),
    'all_recovery_bytes_verified':True,'cleanup_authorized':False}
(WORK/'restoration-protection.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
