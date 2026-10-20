"""Repack identical owned payloads once, preserving every logical path and mode."""
from pathlib import Path
import hashlib,json,stat,zipfile
R=Path.cwd();B=Path(__file__).resolve().parent
m=json.loads((B/'MANIFEST.json').read_text());old=m['archives'];seen={};backings={}
for name,row in m['files'].items():
 key=(row['bytes'],row['sha256'],row['mode'])
 if key not in seen:seen[key]=name
 backings[name]=seen[key]
target=B/'owned-payloads.zip';assert not target.exists()
with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for name in seen.values():
  row=m['files'][name];p=R/name;b=p.read_bytes()
  assert len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256']
  i=zipfile.ZipInfo(name,date_time=(2026,10,21,0,0,0));i.create_system=3;i.compress_type=zipfile.ZIP_DEFLATED;i.external_attr=(stat.S_IFREG|row['mode'])<<16
  z.writestr(i,b)
assert target.stat().st_size<95_000_000
with zipfile.ZipFile(target) as z:
 assert z.testzip() is None and len(z.namelist())==len(seen)
 for name,member in backings.items():
  b=z.read(member);row=m['files'][name]
  assert len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256']
  assert stat.S_IMODE(z.getinfo(member).external_attr>>16)==row['mode']
  row['backing_member']=member
with target.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
m['archives']={target.name:{'bytes':target.stat().st_size,'sha256':digest,'members':len(seen)}}
m['schema']='sirius-viewer-input-reinitialization-preservation-v2'
m['deduplication']='Every original logical file and mode is retained; backing_member selects its exact regular ZIP member. Identical current source/math/resource bytes are stored once, without replacing any recorded execution input.'
(B/'MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n')
(B/'container-repack.json').write_text(json.dumps({'passed':True,'logical_files':len(m['files']),'unique_members':len(seen),'superseded_local_container_encodings_not_retained':old,'scope':'Local preservation-container repack only. Every original payload byte/path/mode is conserved; no provider/execution/qualification data changed.'},indent=2)+'\n')
for name in old:(B/name).unlink()
print(json.dumps({'logical_files':len(m['files']),'unique_members':len(seen),'archive_bytes':target.stat().st_size,'sha256':digest}))
