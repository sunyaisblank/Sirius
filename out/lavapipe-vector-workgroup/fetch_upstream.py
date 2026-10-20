from pathlib import Path
import hashlib,json,tarfile,time,urllib.request
W=Path(__file__).resolve().parent
URL='https://archive.mesa3d.org/mesa-25.2.8.tar.xz'
EXPECTED='097842f3e49d996868b38688db87b006f7d4541e93ce86d2f341d8b3e7be7c93'
p=W/'mesa-25.2.8.tar.xz';assert not p.exists()
start=time.monotonic();size=0;digest=hashlib.sha256()
with urllib.request.urlopen(URL,timeout=60) as response,p.open('xb') as out:
 while block:=response.read(1024*1024):
  size+=len(block);assert size<=100_000_000 and time.monotonic()-start<180
  digest.update(block);out.write(block)
assert digest.hexdigest()==EXPECTED
with tarfile.open(p) as t:
 members=t.getmembers()
 assert all(m.name=='mesa-25.2.8' or m.name.startswith('mesa-25.2.8/') for m in members)
 for m in members:
  assert not Path(m.name).is_absolute() and '..' not in Path(m.name).parts
 t.extractall(W/'upstream',filter='data')
s=W/'upstream/mesa-25.2.8';changed={}
for p in (W/'candidate').rglob('*'):
 if p.is_file():
  n=p.relative_to(W/'candidate');upstream=s/n;original=W/'original'/n
  assert upstream.read_bytes()==original.read_bytes(),n
  changed[str(n)]={'original_sha256':hashlib.sha256(upstream.read_bytes()).hexdigest(),'candidate_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
  upstream.write_bytes(p.read_bytes())
(W/'upstream-admission.json').write_text(json.dumps({'passed':True,'url':URL,'official_checksum_source':'https://docs.mesa3d.org/relnotes/25.2.8.html','archive_bytes':size,'archive_sha256':digest.hexdigest(),'tar_members':len(members),'changed_files':changed,'elapsed_seconds':time.monotonic()-start,'scope':'Official archive checksum and exact original source joins before private patch; no provider build or runtime qualification.'},indent=2)+'\n')
print(json.dumps({'passed':True,'archive_bytes':size,'tar_members':len(members),'changed_files':len(changed),'elapsed_seconds':time.monotonic()-start}))
