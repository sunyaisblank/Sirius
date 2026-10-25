"""Retire only three recovered integer identity-fixture products."""
from pathlib import Path
import hashlib,json,os,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-local-identities';C=R/'out/native-cpu-local-identities-cohort';D=R/'out/native-cpu-local-identities-compatibility';J=lambda p:json.loads(p.read_text())
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{upstream}'])
source=J(C/'control-8/owner.json')['source']['revision'];bundles=list((R/'attestations/qualification'/source[:8]).glob('native-cpu-local-identities-*/MANIFEST.json'));assert len(bundles)==1
manifest_path=bundles[0];manifest=J(manifest_path);assert manifest['finite_component_passed'] and manifest['original_C8_passed']
for p in [manifest_path,*(manifest_path.parent/a['file'] for a in manifest['archives'])]:
 assert subprocess.check_output(['git','show','@{upstream}:'+str(p.relative_to(R))])==p.read_bytes()
for root in (W,C,D):
 for p in root.rglob('owner.json'):
  o=J(p);assert o['status']=='terminal' and not o['remaining_owned_processes'] and not o['cleanup_errors']
  assert o['source_unchanged'] and o['whole_inputs_unchanged']
  for row in o['observed_births']:
   try:t=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
   except (FileNotFoundError,ProcessLookupError):t=None
   assert t!=str(row['saved_start_ticks'])
frozen={(R/name).resolve() for root in (W,C,D,R/'out/native-cpu-owned-lower-increments') for p in root.rglob('actions.json') if p!=W/'fixture/actions.json' for a in J(p).values() for name in a['frozen_inputs']}
tracked=set(subprocess.check_output(['git','ls-files'],text=True).splitlines());members=[];archives=[]
products=[R/x['product']['path'] for x in J(W/'fixture/build-result.json')['commands']]
groups=[(manifest_path,products)]
for p,paths in groups:
 assert subprocess.check_output(['git','show','@{upstream}:'+str(p.relative_to(R))])==p.read_bytes()
 a,=[x for x in J(p)['archives'] if x['file']=='actual-products.zip'];archive=p.parent/a['file']
 assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==a['identity']['sha256']
 assert hashlib.sha256(subprocess.check_output(['git','show','@{upstream}:'+str(archive.relative_to(R))])).hexdigest()==a['identity']['sha256']
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None and set(z.namelist())==set(a['members'])
  for key,v in a['members'].items():
   data=z.read(key);assert len(data)==v['bytes'] and hashlib.sha256(data).hexdigest()==v['sha256']
   assert stat.S_IMODE(z.getinfo(key).external_attr>>16)==v['mode']
  for path in paths:
   relative=str(path.relative_to(R));key='files/'+relative;v=a['members'][key];s=path.stat()
   assert path.is_file() and not path.is_symlink() and relative not in tracked and path.resolve() not in frozen
   assert s.st_nlink==1 and s.st_size==v['bytes'] and stat.S_IMODE(s.st_mode)==v['mode'] and s.st_mtime_ns==v['mtime_ns']
   assert path.read_bytes()==z.read(key)
   members.append(dict(path=relative,archive=str(archive.relative_to(R)),member=key,sha256=v['sha256'],bytes=s.st_size,allocated_bytes=s.st_blocks*512,removed=False))
 archives.append(dict(path=str(archive.relative_to(R)),sha256=a['identity']['sha256'],verified_members=len(a['members'])))
assert len(members)==3;targets={str(R/x['path']) for x in members};checked=0;privileged=[]
code="""from pathlib import Path
import json,os,sys
p=Path(sys.argv[1]);targets=set(json.loads(sys.argv[3]));matches=[];count=0;maps=''
try:
 if p.stat().st_uid==int(sys.argv[2]):
  maps=(p/'maps').read_text()
  for fd in (p/'fd').iterdir():
   try:target=os.readlink(fd)
   except FileNotFoundError:continue
   count+=1
   if target in targets:matches.append(fd.name)
except FileNotFoundError:pass
print(json.dumps(dict(matches=matches,count=count,mapped=[x for x in targets if x in maps])))
"""
for p in Path('/proc').iterdir():
 if not p.name.isdecimal():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  try:
   maps=(p/'maps').read_text();assert not any(target in maps for target in targets)
   for fd in (p/'fd').iterdir():
    try:target=os.readlink(fd)
    except FileNotFoundError:continue
    assert target not in targets,(p,fd)
  except PermissionError:
   result=json.loads(subprocess.check_output(['sudo','-n','/usr/bin/python3','-B','-c',code,str(p),str(os.getuid()),json.dumps(sorted(targets))],text=True))
   assert not result['matches'] and not result['mapped'];privileged.append(dict(pid=int(p.name),descriptors_checked=result['count']))
  checked+=1
 except (FileNotFoundError,ProcessLookupError):continue
build=(R/'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja').read_text();assert not any(target in build for target in targets)
report=dict(status='started',source=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),members=members,recovery_archives=archives,
 same_user_map_and_fd_tables_checked=checked,privileged_read_only_metadata_checks=privileged,gross_allocated_bytes_retired=sum(x['allocated_bytes'] for x in members),
 scope='Three untracked recovered diagnostic products only: renamed original scalar-reference driver, integer-only ISPC fixture object and linked ABI check. Full pushed archive SHA/CRC/member/mode/time verified. No current AN/H non-fixture frozen-input, Ninja, same-user mapped or open-FD consumer. The sole historical fixture action is terminal, all births absent, and preserved in pushed recovery. All current raw/provider/full38 products, unchanged AM math/provider and unique source/concise evidence/resources/reusable caches retained. Gross retired allocation, not net drive gain or scientific/release qualification.')
record=manifest_path.parent/'RETIRED-DIAGNOSTIC-PRODUCTS.json';assert not record.exists()
def save():record.write_text(json.dumps(report,indent=2)+'\n')
save()
for row in members:
 p=R/row['path'];assert hashlib.file_digest(p.open('rb'),'sha256').hexdigest()==row['sha256'];p.unlink();row['removed']=True;save()
(R/'bin/linux-gcc/native-cpu-local-identities/fixture').rmdir()
report['status']='complete';save();print(json.dumps(dict(retired_product_names=3,gross_allocated_bytes=report['gross_allocated_bytes_retired'])))
