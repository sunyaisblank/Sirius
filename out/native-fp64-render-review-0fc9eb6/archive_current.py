"""Preserve one finished native batch and exact input recovery mappings."""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
DEST=ROOT/'attestations/native-vulkan/0fc9eb6/original-fp64-science-and-frame'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def identity(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def key(v):return v['bytes'],v['sha256']
bindings=read(WORK/'execution-bindings.json');revision=bindings['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==revision
assert not subprocess.check_output(['git','status','--porcelain'])
assert not DEST.exists();DEST.mkdir(parents=True)
pool={};tracked={}
for line in subprocess.check_output(['git','ls-tree','-r',revision],text=True).splitlines():
 meta,path=line.split('\t');mode,kind,oid=meta.split();assert kind=='blob'
 data=subprocess.check_output(['git','cat-file','blob',oid])
 v={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
 assert identity(ROOT/path)==v,path
 retained={'kind':'git-blob','revision':revision,'path':path,'object_id':oid,**v}
 tracked[path]=retained;pool.setdefault(key(v),retained)
for path,v in bindings['files'].items():
 if path.startswith('attestations/native-build/0fc9eb6/'):
  assert identity(ROOT/path)==v,path
  pool.setdefault(key(v),{'kind':'protected-official-export','path':path,**v})
header='attestations/native-vulkan/0fc9eb6/paired-transport-restoration/source-build-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
pool.setdefault(key(identity(ROOT/header)),{'kind':'protected-restoration-input','path':header,**identity(ROOT/header)})
scratch={}
for p in sorted(WORK.rglob('*')):
 if not p.is_file():continue
 assert not p.is_symlink(),p
 relative=p.relative_to(WORK);q=DEST/'diagnostic'/relative;q.parent.mkdir(parents=True,exist_ok=True)
 shutil.copyfile(p,q);v=identity(p);assert identity(q)==v
 scratch[str(p.relative_to(ROOT))]={'kind':'bundled-diagnostic','path':str(q.relative_to(ROOT)),**v}
stage={}
for p in sorted((ROOT/bindings['stage']).rglob('*')):
 if not p.is_file():continue
 assert not p.is_symlink(),p
 v=identity(p);retained=pool.get(key(v))
 if retained is None:
  q=DEST/'unmapped-stage'/p.relative_to(ROOT/bindings['stage']);q.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(p,q);assert identity(q)==v
  retained={'kind':'bundled-stage','path':str(q.relative_to(ROOT)),**v}
 stage[str(p.relative_to(ROOT))]={'original':v,'retained':retained}
input_ledger={}
sealed=read(WORK/'post-frame-whole-input-seals.json')
for group,items in sealed['files'].items():
 for item in items:
  path=item['path'];p=Path(path) if path.startswith('/') else ROOT/path
  v={k:item[k] for k in ['bytes','sha256']};assert identity(p)==v,path
  retained=scratch.get(path) or tracked.get(path)
  if retained is None and path in stage:retained=stage[path]['retained']
  if retained is None:retained=pool.get(key(v))
  if retained is None and path.startswith('/mnt/c/'):
   retained={'kind':'installed-provider-retained','path':path,**v}
  assert retained is not None,path
  assert key(retained)==key(v),path
  input_ledger.setdefault(path,{'original':v,'retained':retained,'groups':[]})['groups'].append(group)
assert len(bindings['files'])==734 and len(stage)==94
assert set(bindings['files'])<=set(input_ledger)
ledger={'source_revision':revision,'frozen_count':734,'post_frame_predecessor_count':sealed['predecessor_count'],'frozen_and_predecessor_inputs':input_ledger,'all_stage_files':stage,'all_scratch_files':scratch,'scope':'Exact whole-byte recovery mapping; installed provider references remain installed. No native CTest reconstruction or full runtime/release qualification.'}
(DEST/'input-retention-ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
manifest={'source_revision':revision,'kind':'Original current FP64 science pass and incomplete frame','files':{}}
for p in sorted(DEST.rglob('*')):
 if p.is_file():manifest['files'][str(p.relative_to(DEST))]=identity(p)
(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'archive':str(DEST.relative_to(ROOT)),'payload_count':len(manifest['files']),'payload_bytes':sum(v['bytes'] for v in manifest['files'].values()),'manifest':identity(DEST/'manifest.json'),'stage_count':len(stage),'unmapped_stage_count':sum(v['retained']['kind']=='bundled-stage' for v in stage.values()),'input_count':len(input_ledger)}))
