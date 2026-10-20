"""Retire only the exact pushed local viewer-input diagnostic roots."""
from pathlib import Path
import hashlib,json,os,shutil,stat,subprocess,zipfile
R=Path.cwd();B=Path(__file__).resolve().parent
assert B==R/'attestations/qualification/6d2af336/viewer-input-reinitialization'
def git(*args):return subprocess.check_output(['git',*args])
def doc(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def ident(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
head=git('rev-parse','HEAD').decode().strip();branch=git('branch','--show-current').decode().strip()
assert branch=='development/viewer-input-evidence-2026-10-21'
assert not git('status','--porcelain')
assert git('ls-remote','--heads','origin',branch).decode().split()==[head,'refs/heads/'+branch]
for p in B.rglob('*'):
 if p.is_file():assert p.read_bytes()==git('show',head+':'+str(p.relative_to(R)))
assert doc(B/'ci-review.json')['passed']
m=doc(B/'MANIFEST.json');assert m['schema']=='sirius-viewer-input-reinitialization-preservation-v2'
name,archive=next(iter(m['archives'].items()));assert ident(B/name)=={k:archive[k] for k in ['bytes','sha256']}
roots=[R/n for n in m['roots']]
assert all(p.resolve(strict=True)==p and not p.is_symlink() for p in roots)
assert {str(p.relative_to(R)) for root in roots for p in root.rglob('*') if p.is_file()}==set(m['files'])
with zipfile.ZipFile(B/name) as z:
 assert z.testzip() is None and len(z.namelist())==archive['members']
 for name,row in m['files'].items():
  p=R/name;assert ident(p)=={k:row[k] for k in ['bytes','sha256']}
  assert not p.is_symlink() and stat.S_IMODE(p.stat().st_mode)==row['mode']
  b=z.read(row['backing_member']);assert len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256']
  assert stat.S_IMODE(z.getinfo(row['backing_member']).external_attr>>16)==row['mode']
births=set()
def walk(v):
 if isinstance(v,list):
  for x in v:walk(x)
 elif isinstance(v,dict):
  for key in ['saved_start_ticks','start_ticks']:
   if type(v.get('pid')) is int and isinstance(v.get(key),str) and v[key].isdigit():births.add((v['pid'],v[key]))
  for x in v.values():walk(x)
for n in m['files']:
 if n.endswith('.json'):
  try:walk(doc(R/n))
  except (ValueError,UnicodeError):pass
observed=[]
for pid,saved in sorted(births):
 try:current=Path('/proc')/str(pid)/'stat';ticks=current.read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):ticks=None
 assert ticks!=saved,(pid,saved,ticks)
 observed.append({'pid':pid,'saved_birth':saved,'current_birth':ticks,'original_absent':True})
consumers=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit() or int(p.name)==os.getpid():continue
 candidates=[]
 for key in ['cwd','exe']:
  try:candidates.append(os.readlink(p/key))
  except OSError:pass
 try:
  for f in (p/'fd').iterdir():
   try:candidates.append(os.readlink(f))
   except OSError:pass
 except OSError:pass
 for value in candidates:
  if any(value==str(root) or value.startswith(str(root)+'/') for root in roots):consumers.append({'pid':p.name,'path':value})
assert not consumers,consumers
for root in roots:shutil.rmtree(root)
assert all(not p.exists() for p in roots)
result={'schema':'sirius-viewer-input-retirement-v1','passed':True,'preserved_remote_head_before_cleanup':head,'pushed_Git_blob_and_ZIP_logical_member_byte_mode_checks_passed':True,'removed_roots':m['roots'],'removed_files':len(m['files']),'unique_archive_members':archive['members'],'recorded_linux_births_absent':observed,'accessible_owned_root_consumers':consumers,'canonical_presets_and_toolchains_retained':True,'full_qualification_claimed':False}
(B/'retirement.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':True,'removed_roots':len(roots),'logical_files':len(m['files']),'recorded_births':len(observed),'preserved_remote_head':head}))
