"""Retire only recovered, untracked AK ledgers outside current run closures."""
from pathlib import Path
import hashlib,json,os,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-predicated-refinement';C=R/'out/native-cpu-predicated-refinement-cohort'
B=R/'attestations/qualification/5bc67064/native-cpu-predicated-refinement-finite-compatibility'
J=lambda p:json.loads(p.read_text())
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{upstream}'])
manifest=J(B/'MANIFEST.json')
assert manifest['finite_compatibility_passed'] and not manifest['accepted'] and not manifest['frame_executed']
assert len(manifest['owners'])==8 and all(o['passed'] and o['status']=='terminal' for o in manifest['owners'].values())
for name in manifest['owners']:
 o=J(R/name);assert not o['remaining_owned_processes'] and not o['cleanup_errors']
 for row in o['observed_births']:
  try:ticks=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
  except (FileNotFoundError,ProcessLookupError):ticks=None
  assert ticks!=str(row['saved_start_ticks'])
for p in [B/'MANIFEST.json',*(B/a['file'] for a in manifest['archives'])]:
 assert subprocess.check_output(['git','show','@{upstream}:'+str(p.relative_to(R))])==p.read_bytes()
metadata,=[a for a in manifest['archives'] if a['file']=='source-build-and-controls.zip'];archive=B/metadata['file']
assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==metadata['identity']['sha256']
frozen={(R/name).resolve() for p in [C/'actions.json',C/'affected-control/actions.json'] for a in J(p).values() for name in a['frozen_inputs']}
tracked=set(subprocess.check_output(['git','ls-files'],text=True).splitlines());rows=[]
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None and set(z.namelist())==set(metadata['members'])
 for phase in [W/'cpp/build',W/'cpp/run',W/'provider/configure',W/'provider/build']:
  for leaf in ['inputs-before.json','inputs-after.json','samples.jsonl']:
   p=phase/leaf;name=str(p.relative_to(R));key='files/'+name;v=metadata['members'][key];s=p.stat()
   assert p.is_file() and not p.is_symlink() and name not in tracked and p.resolve() not in frozen
   assert s.st_size==v['bytes'] and stat.S_IMODE(s.st_mode)==v['mode'] and s.st_mtime_ns==v['mtime_ns']
   assert p.read_bytes()==z.read(key) and hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256']
   rows.append(dict(path=name,sha256=v['sha256'],bytes=s.st_size,allocated_bytes=s.st_blocks*512,removed=False))
assert len(rows)==12
targets={str((R/r['path']).resolve()) for r in rows}
checked=0;privileged=[]
metadata_check="""from pathlib import Path
import json,os,sys
p=Path(sys.argv[1]);targets=set(json.loads(sys.argv[3]));matches=[];count=0
try:
 if p.stat().st_uid==int(sys.argv[2]):
  for fd in (p/'fd').iterdir():
   try:target=os.readlink(fd)
   except FileNotFoundError:continue
   count+=1
   if target in targets:matches.append(fd.name)
except FileNotFoundError:pass
print(json.dumps(dict(matches=matches,descriptors_checked=count)))
"""
for p in Path('/proc').iterdir():
 if not p.name.isdecimal():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  try:
   for fd in (p/'fd').iterdir():
    try:target=os.readlink(fd)
    except FileNotFoundError:continue
    assert target not in targets,(p,fd,target)
  except PermissionError:
   value=json.loads(subprocess.check_output(['sudo','-n','/usr/bin/python3','-B','-c',metadata_check,str(p),str(os.getuid()),json.dumps(sorted(targets))],text=True))
   assert not value['matches'],(p,value)
   privileged.append(dict(pid=int(p.name),descriptors_checked=value['descriptors_checked']))
  checked+=1
 except (FileNotFoundError,ProcessLookupError):continue
report=dict(status='started',source=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 recovered_archive=str(archive.relative_to(R)),archive_sha256=metadata['identity']['sha256'],members=rows,
 current_same_user_fd_tables_checked=checked,privileged_read_only_fd_link_checks=privileged,gross_allocated_bytes_retired=sum(r['allocated_bytes'] for r in rows),
 scope='Twelve expanded compiler/configure/build ledgers only, whole SHA/CRC/member/mode/time verified against pushed source archive, no live owned birth/open same-user FD or current C8/C5 frozen-input consumer. Actual products, original/AJ sources, current C8/C5/full38 records, all eight owner outcomes, useful caches and scoped reconstruction recipe remain. Gross retired allocation, not net drive gain or qualification.')
record=B/'RETIRED-EXPANDED-LEDGERS.json';assert not record.exists()
def save():record.write_text(json.dumps(report,indent=2)+'\n')
save()
for r in rows:
 p=R/r['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];p.unlink();r['removed']=True;save()
report['status']='complete';save()
print(json.dumps(dict(retired_files=len(rows),gross_allocated_bytes=report['gross_allocated_bytes_retired'])))
