"""Preserve the finished finite #94 trial before overwriting build outputs."""
from pathlib import Path
import hashlib,json,shutil,subprocess
R=Path.cwd(); W=Path(__file__).resolve().parent
def read(path):return json.loads(Path(path).read_text())
def identity(path):
    data=Path(path).read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
head=read(W/'retention-gate-frozen.json')['candidate_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
assert read(W/'primitive-result.json')['pass']
assert read(W/'software-sequence.json')['completed']
comparison=read(W/'software-comparison.json')
assert len(comparison['checks'])==20 and comparison['readback_identity_module_and_completion_prerequisites_pass']
dest=R/'attestations/software-vulkan'/head[:7]/'guarded-upmultiply-trial'
assert not dest.exists()
tree={}
for row in subprocess.check_output(['git','ls-tree','-r','-z',head]).split(b'\0'):
    if row:
        info,path=row.split(b'\t',1); mode,kind,oid=info.decode().split();assert kind=='blob';tree[path.decode()]=oid
versions={}; selected={}
for path in W.rglob('*'):
    if path.is_file() and not path.is_symlink():
        value=identity(path); versions[(value['bytes'],value['sha256'])]=path
        selected['diagnostic/'+str(path.relative_to(W))]=path
ledgers={}; git_values={}
cases=['configure','build','prepare-primitive','admission','primitive-baseline','primitive-candidate',
    'transport-control','endpoint-control','projection-control','wide-control','a1','b1','b2','a2']
for case in cases:
    folder=W/case; owner=read(folder/'owner.json')
    assert owner['passed'] and owner['child_exit']==0 and owner['status']=='terminal'
    assert owner['source']=={'revision':head,'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and owner['remaining_owned_processes']==[]
    assert owner['observed_births'] and all(item['absent'] for item in owner['observed_births'])
    before=read(folder/'inputs-before.json');assert before==read(folder/'inputs-after.json')
    ledger=[]
    for origin,seal in sorted(before.items()):
        value={key:seal[key] for key in ['bytes','sha256']}
        key=(value['bytes'],value['sha256'])
        if origin in tree:
            if origin not in git_values:
                data=subprocess.check_output(['git','cat-file','blob',tree[origin]])
                git_values[origin]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
            assert git_values[origin]==value
            recovery={'kind':'immutable-Git-blob','revision':head,'blob':tree[origin],'path':origin}
        else:
            live=Path(origin) if Path(origin).is_absolute() else R/origin
            if live.is_file() and identity(live)==value:
                assert str(live.resolve(strict=True))==seal['resolved_path']
                source=live
            else:
                assert key in versions,('No exact preserved historical input version',case,origin,value)
                source=versions[key]
            target='executed-inputs/'+value['sha256']
            if target in selected:assert identity(selected[target])==value
            selected[target]=source
            recovery={'kind':'local-archive-payload','path':target}
        ledger.append({'origin':origin,**value,'resolved_path_at_execution':seal['resolved_path'],'recovery':recovery})
    ledgers[case]={'source_revision':head,'inputs':ledger,'passed':True}
payloads=[]
for target,source in sorted(selected.items()):
    value=identity(source); output=dest/target;output.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,output);assert identity(output)==value
    payloads.append({'path':target,'origin':str(source.relative_to(R)) if source.is_relative_to(R) else str(source),**value})
assert {str(path.relative_to(dest)) for path in dest.rglob('*') if path.is_file()}=={item['path'] for item in payloads}
for case in cases:
    for item in ledgers[case]['inputs']:
        if item['recovery']['kind']=='local-archive-payload':assert identity(dest/item['recovery']['path'])=={key:item[key] for key in ['bytes','sha256']}
manifest={'source_revision':head,'baseline_revision':comparison['baseline_revision'],
    'retention_pass':comparison['pass'],'payloads':payloads,'payload_count':len(payloads),
    'payload_bytes':sum(item['bytes'] for item in payloads),'executed_cases':ledgers,
    'scope':'Recoverable finite #94 configure/build/preparation, exact-function words, original affected controls and one software ABBA. Earlier diagnostic versions retain original byte/epoch bounds. Provider/ICD bytes and loaded mappings preserved; no all-external-library/process/consumer closure, original frame/full-science/native/release qualification claim.'}
(dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'archive':str(dest.relative_to(R)),'manifest':identity(dest/'manifest.json'),
    'payload_count':len(payloads),'payload_bytes':manifest['payload_bytes'],
    'case_input_counts':{case:len(item['inputs']) for case,item in ledgers.items()},'fresh_byte_checks_pass':True}
(W/'trial-protection.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
