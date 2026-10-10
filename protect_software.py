"""Retain completed software cases and recoverable whole-input joins."""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def identity(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
source=read(WORK/'source.json');head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
DEST=ROOT/'attestations/review'/head[:7]/'submission-fence-restoration';assert not DEST.exists()
inputs={};selected={}
def select(p,target):
    assert p.is_file() and not p.is_symlink()
    assert target not in selected or selected[target]==p;selected[target]=p
for name in ('configure','build','lifetime','pair','timestamps','shared-tracer'):
    folder=WORK/('candidate-'+name);owner=read(folder/'owner.json')
    assert owner['passed'] and not owner['stop_reason'] and not owner['cleanup_errors'] and owner['remaining_owned_processes']==[] and all(e['absent'] for e in owner['observed_births'])
    before=read(folder/'inputs-before.json');assert before==read(folder/'inputs-after.json')
    for path,record in before.items():
        assert record['resolved_path']==str((ROOT/path).resolve())
        value={k:record[k] for k in ('bytes','sha256')}
        assert path not in inputs or inputs[path]==value;inputs[path]=value
    for p in folder.rglob('*'):
        if p.is_file():select(p,'diagnostic/'+str(p.relative_to(WORK)))
for path,value in read(WORK/'linux-readonly-payloads.json')['actual_consumers'].items():
    value={k:value[k] for k in ('bytes','sha256')};assert path not in inputs or inputs[path]==value;inputs[path]=value
tree={}
for row in subprocess.check_output(['git','ls-tree','-r','-z',head]).split(b'\0'):
    if row:
        info,path=row.split(b'\t',1);mode,kind,oid=info.decode().split();assert kind=='blob';tree[path.decode()]=oid
tracked=sorted(set(inputs)&set(tree))
actual=subprocess.check_output(['git','hash-object','--no-filters','--stdin-paths'],input=''.join(p+'\n' for p in tracked),text=True).splitlines()
assert len(actual)==len(tracked) and all(oid==tree[p] for p,oid in zip(tracked,actual))
ledger=[]
for path,value in sorted(inputs.items()):
    p=ROOT/path;assert identity(p)==value,path
    if path in tree: recovery={'kind':'immutable-Git-blob','revision':head,'blob':tree[path],'path':path}
    else:
        target='executed-inputs/'+path;select(p,target);recovery={'kind':'local-archive-payload','path':target}
    ledger.append({'origin':path,**value,'recovery':recovery})
for name in ('source.json','software-controls.json','build-and-arrays.json','linux-readonly-payloads.json','validate_results.py','protect_software.py'):
    select(WORK/name,'diagnostic/'+name)
payloads=[]
for target,p in sorted(selected.items()):
    value=identity(p);q=DEST/target;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);assert identity(q)==value
    payloads.append({'path':target,'origin':str(p.relative_to(ROOT)),**value})
assert {str(p.relative_to(DEST)) for p in DEST.rglob('*') if p.is_file()}=={v['path'] for v in payloads}
for v in payloads:assert identity(DEST/v['path'])=={k:v[k] for k in ('bytes','sha256')}
for path,value in inputs.items():assert identity(ROOT/path)==value
manifest={'revision':head,'payloads':payloads,'payload_count':len(payloads),'payload_bytes':sum(v['bytes'] for v in payloads),'executed_input_ledger':ledger,'scope':'Exact accepted tree restoration, six strict configure/build/original focused software executions and complete160 unchanged kernel joins. Git source and ignored input recovery; no full frame/science/runtime/release claim.'}
(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'archive':str(DEST.relative_to(ROOT)),'manifest':identity(DEST/'manifest.json'),'payload_count':len(payloads),'payload_bytes':manifest['payload_bytes'],'recoverable_inputs':len(ledger),'Git_blob_inputs':len(tracked),'fresh_byte_checks_pass':True}
(WORK/'software-protection.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
