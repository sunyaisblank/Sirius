"""Retire recovered AJ/AK provider pairs after the corrected AL cohort is saved."""
from pathlib import Path
import hashlib,json,os,stat,subprocess,zipfile
R=Path.cwd();C=R/'out/native-cpu-reachable-tag-writers-cohort';W=R/'out/native-cpu-reachable-tag-writers';D=R/'out/native-cpu-reachable-tag-writers-compatibility';J=lambda p:json.loads(p.read_text())
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{upstream}'])
source=J(C/'control-8/owner.json')['source']['revision'];bundles=list((R/'attestations/qualification'/source[:8]).glob('native-cpu-reachable-tag-writers-*/MANIFEST.json'));assert len(bundles)==1
manifest_path=bundles[0];manifest=J(manifest_path);assert J(W/'fixture/result.json')['old_stale_tag_arithmetic_counterexamples']==12
assert manifest['all22_original_gates_pass'] and manifest['all8_original_controls_pass']
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
frozen={(R/name).resolve() for root in (W,C,D,R/'out/native-cpu-owned-lower-increments') for p in root.rglob('actions.json') for a in J(p).values() for name in a['frozen_inputs']}
tracked=set(subprocess.check_output(['git','ls-files'],text=True).splitlines());members=[];archives=[]
for slug,name in [('native-cpu-normalization-tags','attestations/qualification/0ce77927/native-cpu-normalization-tags-frame-timeout/MANIFEST.json'),('native-cpu-predicated-refinement','attestations/qualification/5bc67064/native-cpu-predicated-refinement-finite-compatibility/MANIFEST.json')]:
 p=R/name;assert subprocess.check_output(['git','show','@{upstream}:'+name])==p.read_bytes()
 a,=[x for x in J(p)['archives'] if x['file']=='actual-products.zip'];archive=p.parent/a['file']
 assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==a['identity']['sha256']
 assert hashlib.sha256(subprocess.check_output(['git','show','@{upstream}:'+str(archive.relative_to(R))])).hexdigest()==a['identity']['sha256']
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None and set(z.namelist())==set(a['members'])
  for key,v in a['members'].items():
   data=z.read(key);assert len(data)==v['bytes'] and hashlib.sha256(data).hexdigest()==v['sha256']
   assert stat.S_IMODE(z.getinfo(key).external_attr>>16)==v['mode']
  for leaf in ['libvulkan_lvp.so','lvp_icd.json']:
   path=R/'bin/linux-gcc'/slug/'provider'/leaf;relative=str(path.relative_to(R));key='files/'+relative;v=a['members'][key];s=path.stat()
   assert path.is_file() and not path.is_symlink() and relative not in tracked and path.resolve() not in frozen
   assert s.st_size==v['bytes'] and stat.S_IMODE(s.st_mode)==v['mode'] and s.st_mtime_ns==v['mtime_ns']
   assert path.read_bytes()==z.read(key)
   members.append(dict(path=relative,archive=str(archive.relative_to(R)),member=key,sha256=v['sha256'],bytes=s.st_size,allocated_bytes=s.st_blocks*512,removed=False))
 archives.append(dict(path=str(archive.relative_to(R)),sha256=a['identity']['sha256'],verified_members=len(a['members'])))
assert len(members)==4;targets={str(R/x['path']) for x in members};checked=0;privileged=[]
for p in Path('/proc').iterdir():
 if not p.name.isdecimal():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  try:maps=(p/'maps').read_text()
  except PermissionError:
   maps=subprocess.check_output(['sudo','-n','/usr/bin/cat',str(p/'maps')],text=True);privileged.append(int(p.name))
  assert not any(target in maps for target in targets),p;checked+=1
 except (FileNotFoundError,ProcessLookupError):continue
build=(R/'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja').read_text();assert not any(target in build for target in targets)
report=dict(status='started',source=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),members=members,recovery_archives=archives,
 same_user_map_tables_checked=checked,privileged_read_only_map_checks=privileged,gross_allocated_bytes_retired=sum(x['allocated_bytes'] for x in members),
 scope='Four untracked superseded AJ/AK provider names whose coherent-tag claims are disproved. All old archive payloads independently verified against pushed recovery; no current corrected AL/owned-host frozen-input, Ninja or same-user mapping consumer. Old source, raw/CPP objects, finite records and authentic recovery retained. Corrected AL products and all scientific/public resources/caches retained. Gross retired allocation, not net drive gain or release qualification.')
record=manifest_path.parent/'RETIRED-INVALID-PROVIDERS.json';assert not record.exists()
def save():record.write_text(json.dumps(report,indent=2)+'\n')
save()
for row in members:
 p=R/row['path'];assert hashlib.file_digest(p.open('rb'),'sha256').hexdigest()==row['sha256'];p.unlink();row['removed']=True;save()
for slug in ['native-cpu-normalization-tags','native-cpu-predicated-refinement']:(R/'bin/linux-gcc'/slug/'provider').rmdir()
report['status']='complete';save();print(json.dumps(dict(retired_provider_names=4,gross_allocated_bytes=report['gross_allocated_bytes_retired'])))
