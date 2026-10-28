"""Retire only the archived, rejected, never-adopted compact DP products."""
from pathlib import Path
import hashlib,json,os,subprocess,zipfile

R=Path.cwd();B=R/'bin/linux-gcc/native-cpu-layer-projection'
A=R/'attestations/qualification/b2cbd4db/native-cpu-private-context-interface-gates'
J=R/'out/native-cpu-layer-projection/retirement.json'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def save(d):J.write_text(json.dumps(d,indent=2)+'\n')


assert not J.exists() and B.resolve()==R/'bin/linux-gcc/native-cpu-layer-projection'
m=json.loads((A/'MANIFEST.json').read_text());recover={}
for bundle in m['bundles'].values():
 p=R/bundle['path'];assert pin(p)=={k:bundle[k] for k in ['bytes','sha256']}
 subprocess.run(['git','ls-files','--error-unmatch',str(p.relative_to(R))],check=True,stdout=subprocess.DEVNULL)
 with zipfile.ZipFile(p) as z:
  for name,expected in bundle['members'].items():
   if name.startswith(str(B.relative_to(R))+'/'):
    target=R/name;assert target.is_file() and not target.is_symlink() and target.parent==B
    assert target.stat().st_nlink==1
    assert pin(target)==expected
    with z.open(name) as s:assert hashlib.file_digest(s,'sha256').hexdigest()==expected['sha256']
    recover[name]=dict(identity=expected,bundle=bundle['path'],member=name,allocated_bytes=target.stat().st_blocks*512)
assert set(recover)=={str(p.relative_to(R)) for p in B.iterdir()} and len(recover)==6
for name in ['plugin','lower','cleanup','native']:
 o=json.loads((R/'out/native-cpu-layer-projection'/name/'owner.json').read_text())
 assert o['passed'] and o['status']=='terminal' and not o['remaining_owned_processes'] and not o['cleanup_errors']
consumers=[];unreadable=0
for proc in Path('/proc').glob('[0-9]*'):
 try:
  if proc.stat().st_uid!=os.getuid():continue
  for fd in (proc/'fd').iterdir():
   try:
    if os.readlink(fd).startswith(str(B)+'/'):consumers.append(dict(pid=int(proc.name),fd=fd.name))
   except (FileNotFoundError,ProcessLookupError):pass
  for line in (proc/'maps').read_text().splitlines():
   if str(B)+'/' in line:consumers.append(dict(pid=int(proc.name),mapping=True))
 except (FileNotFoundError,ProcessLookupError):pass
 except PermissionError:unreadable+=1
assert not consumers and not unreadable,(consumers,unreadable)
before=os.statvfs(R);free_before=before.f_bavail*before.f_frsize
record=dict(status='started',source_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 manifest=dict(path=str((A/'MANIFEST.json').relative_to(R)),**pin(A/'MANIFEST.json')),
 verified_recovery=recover,live_same_uid_consumers=consumers,unreadable_same_uid_processes=unreadable,
 original_semantic_source_and_receipts_retained=True,active_math_and_SDK_caches_preserved=True,unlinked=[])
save(record)
for name,expected in recover.items():
 target=R/name;assert pin(target)==expected['identity'];target.unlink();record['unlinked'].append(name);save(record)
B.rmdir();after=os.statvfs(R)
record.update(status='complete',allocated_bytes_retired=sum(x['allocated_bytes'] for x in recover.values()),
 observed_free_bytes_change=after.f_bavail*after.f_frsize-free_before,removed_directory=str(B.relative_to(R)))
save(record);print(json.dumps({k:record[k] for k in ['status','allocated_bytes_retired','observed_free_bytes_change','removed_directory']}))
