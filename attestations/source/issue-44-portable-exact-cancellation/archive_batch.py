"""Seal complete owned raw evidence; refuse live owned work."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys
r=Path.cwd();w=Path(__file__).resolve().parent;b=r/'bin/linux-gcc/portable-exact-cancellation';a=r/'attestations/software-vulkan/1d59f17/portable-exact-cancellation'
assert sys.argv[1:]==['--archive'];assert not a.exists()
identity=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
receipts=[]
for p in w.rglob('owner.json'):
 v=json.loads(p.read_text());assert v['status']=='terminal' and v['remaining_owned_processes']==[] and all(x['absent'] for x in v['observed_births']),p;receipts.append(str(p.relative_to(r)))
assert receipts
for name,v in json.loads((w/'protected-before.json').read_text()).items():assert identity(r/name)==v,name
disposition=json.loads((w/'disposition.json').read_text());assert disposition['disposition']=='REJECT_FOR_ADOPTION';assert disposition['all_owned_execution_terminal']
a.mkdir(parents=True)
manifest={};modes={}
for root in (w,b):
 for p in sorted(root.rglob('*')):
  if not p.is_file():continue
  assert not p.is_symlink(),p
  relative=p.relative_to(r);target=a/'files'/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);v=identity(p);assert identity(target)==v;manifest[str(relative)]=v;modes[str(relative)]=oct(p.stat().st_mode&0o777);assert (target.stat().st_mode&0o777)==(p.stat().st_mode&0o777)
record={'source_revision':'1d59f1742e43c55019aabba9dce4073ad8ee6e4e','roots':[str(w.relative_to(r)),str(b.relative_to(r))],'files':manifest,'preserved_file_modes':modes,'file_count':len(manifest),'total_bytes':sum(x['bytes'] for x in manifest.values()),'terminal_owner_receipts':receipts,'canonical_protected_files':len(json.loads((w/'protected-before.json').read_text())),'meaning':'Whole raw owned evidence exact recovery; not qualification or source adoption; no remote raw backup claim'}
(a/'whole-byte-manifest.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'archive':str(a.relative_to(r)),'files':len(manifest),'bytes':record['total_bytes'],'terminal_receipts':len(receipts)}))
