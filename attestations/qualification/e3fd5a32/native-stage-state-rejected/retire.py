"""Verify pushed exact backings and retired owners before owned-root removal."""
from pathlib import Path
import hashlib,json,os,shutil,stat,subprocess,zipfile
R=Path.cwd();B=Path(__file__).resolve().parent
assert B==R/'attestations/qualification/e3fd5a32/native-stage-state-rejected'
def doc(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def identity(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
def git(*args):return subprocess.check_output(['git',*args])
m=doc(B/'MANIFEST.json');head=git('rev-parse','HEAD').decode().strip()
branch=git('branch','--show-current').decode().strip()
assert branch=='development/native-stage-state-evidence-2026-10-21'
assert not git('status','--porcelain')
remote=git('ls-remote','--heads','origin',branch).decode().split()
assert remote==[head,'refs/heads/'+branch]
for p in B.iterdir():
 if p.is_file():
  assert p.read_bytes()==git('show',head+':'+str(p.relative_to(R)))
roots=[R/n for n in m['roots']]
assert all(p.resolve(strict=True)==p and p.is_dir() and not p.is_symlink() for p in roots)
assert {str(p.relative_to(R)) for root in roots for p in root.rglob('*') if p.is_file()}==set(m['files'])
for name,row in m['files'].items():
 p=R/name
 assert not p.is_symlink() and identity(p)=={k:row[k] for k in ('bytes','sha256')},name
 assert stat.S_IMODE(p.stat().st_mode)==row['mode'],name
zips={}
for name,row in m['archives'].items():
 p=B/name;assert identity(p)=={k:row[k] for k in ('bytes','sha256')}
 with zipfile.ZipFile(p) as z:
  assert z.testzip() is None and len(z.namelist())==row['members']
  for i in z.infolist():
   assert i.filename not in zips and i.filename in m['files']
   b=z.read(i);v=m['files'][i.filename]
   assert len(b)==v['bytes'] and hashlib.sha256(b).hexdigest()==v['sha256']
   assert stat.S_IMODE(i.external_attr>>16)==v['mode']
   zips[i.filename]=name
assert set(zips)==set(m['files'])
assert doc(R/'out/native-stage-state-ci/precleanup-native-consumer-review.json')['passed'] is True
births=set()
def walk(v):
 if isinstance(v,list):
  for x in v:walk(x)
 elif isinstance(v,dict):
  for pk,bk in [('linux_controller_pid','linux_controller_birth'),('linux_child_pid','linux_child_birth'),('pid','saved_start_ticks'),('pid','saved_birth')]:
   pid,b=v.get(pk),v.get(bk)
   if type(pid) is int and isinstance(b,str) and b.isdigit():births.add((pid,b))
  for x in v.values():walk(x)
for n in m['files']:
 if n.endswith('.json'):
  try:walk(doc(R/n))
  except (UnicodeError,ValueError):pass
observed=[]
for pid,saved in sorted(births):
 try:current=Path('/proc')/str(pid)/'stat';ticks=current.read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):ticks=None
 assert ticks!=saved,(pid,saved,ticks)
 observed.append({'pid':pid,'saved_birth':saved,'current_birth':ticks,'original_absent':True})
accessible=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit() or int(p.name)==os.getpid():continue
 candidates=[]
 for key in ['cwd','exe']:
  try:candidates.append(os.readlink(p/key))
  except (OSError,PermissionError):pass
 try:
  for f in (p/'fd').iterdir():
   try:candidates.append(os.readlink(f))
   except OSError:pass
 except OSError:pass
 for name in candidates:
  if any(name==str(root) or name.startswith(str(root)+'/') for root in roots):accessible.append({'pid':int(p.name),'path':name})
assert not accessible,accessible
# The source/output is outside every removed root and is already pushed above.
for root in roots:shutil.rmtree(root)
assert all(not p.exists() for p in roots)
result={'schema':'sirius-native-stage-state-retirement-v1','passed':True,'preserved_remote_head_before_cleanup':head,'all_Git_blob_archive_member_current_byte_mode_checks_passed':True,'removed_roots':m['roots'],'removed_files':len(m['files']),'linux_original_births_absent':observed,'accessible_linux_owned_root_consumers':accessible,'native_births_and_accessible_owned_path_consumers_absent':True,'native_review_scope':'17 exact saved births and accessible Win32_Process paths at the retained precleanup observation. Hidden protected processes are outside this claim.','canonical_presets_and_toolchains_retained':True,'production_shader_adopted':False,'full_qualification_claimed':False}
(B/'retirement.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':True,'removed_roots':len(roots),'removed_files':len(m['files']),'linux_births':len(observed),'preserved_remote_head':head}))
