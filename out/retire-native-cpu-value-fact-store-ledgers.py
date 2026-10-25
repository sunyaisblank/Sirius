"""Retire recovered compiler/build ledgers outside the retained cohort closure."""
from pathlib import Path
import hashlib,json,os,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-value-fact-stores';C=R/'out/native-cpu-value-fact-stores-cohort';S=R/'out/native-cpu-value-facts';D=R/'out/native-cpu-value-fact-stores-compatibility';J=lambda p:json.loads(p.read_text())
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{upstream}'])
source=J(C/'control-8/owner.json')['source']['revision'];paths=list((R/'attestations/qualification'/source[:8]).glob('native-cpu-value-fact-stores-*/MANIFEST.json'));assert len(paths)==1
B=paths[0].parent;manifest=J(paths[0]);assert manifest['all22_original_gates_pass'] and manifest['all8_original_controls_pass']
for name in manifest['owners']:
 o=J(R/name);assert o['status']=='terminal' and not o['remaining_owned_processes'] and not o['cleanup_errors']
 for row in o['observed_births']:
  try:t=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
  except (FileNotFoundError,ProcessLookupError):t=None
  assert t!=str(row['saved_start_ticks'])
for p in [paths[0],*(B/a['file'] for a in manifest['archives'])]:
 assert subprocess.check_output(['git','show','@{upstream}:'+str(p.relative_to(R))])==p.read_bytes()
metadata,=[a for a in manifest['archives'] if a['file']=='source-build-and-trial.zip'];archive=B/metadata['file']
assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==metadata['identity']['sha256']
frozen={(R/name).resolve() for root in (C,D,R/'out/native-cpu-owned-lower-increments') for p in root.rglob('actions.json') for a in J(p).values() for name in a['frozen_inputs']}
tracked=set(subprocess.check_output(['git','ls-files'],text=True).splitlines());rows=[]
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None and set(z.namelist())==set(metadata['members'])
 for root,phase in [(S,'cpp/build'),(S,'cpp/run'),(W,'cpp/build'),(W,'cpp/run'),(W,'fixture/build'),(W,'fixture/run'),(W,'compile/ispc-dense'),(W,'provider/configure'),(W,'provider/build')]:
  for leaf in ['inputs-before.json','inputs-after.json','samples.jsonl']:
   p=root/phase/leaf;name=str(p.relative_to(R));v=metadata['members']['files/'+name];s=p.stat()
   assert p.is_file() and not p.is_symlink() and name not in tracked and p.resolve() not in frozen
   assert s.st_size==v['bytes'] and stat.S_IMODE(s.st_mode)==v['mode'] and s.st_mtime_ns==v['mtime_ns']
   assert p.read_bytes()==z.read('files/'+name) and hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256']
   rows.append(dict(path=name,sha256=v['sha256'],bytes=s.st_size,allocated_bytes=s.st_blocks*512,removed=False))
assert len(rows)==27;targets={str(R/r['path']) for r in rows};checked=0;privileged=[]
code="""from pathlib import Path
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
print(json.dumps(dict(matches=matches,count=count)))
"""
for p in Path('/proc').iterdir():
 if not p.name.isdecimal():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  try:
   for fd in (p/'fd').iterdir():
    try:target=os.readlink(fd)
    except FileNotFoundError:continue
    assert target not in targets,(p,fd)
  except PermissionError:
   result=json.loads(subprocess.check_output(['sudo','-n','/usr/bin/python3','-B','-c',code,str(p),str(os.getuid()),json.dumps(sorted(targets))],text=True))
   assert not result['matches'];privileged.append(dict(pid=int(p.name),descriptors_checked=result['count']))
  checked+=1
 except (FileNotFoundError,ProcessLookupError):continue
report=dict(status='started',source=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),recovered_archive=str(archive.relative_to(R)),archive_sha256=metadata['identity']['sha256'],members=rows,
 same_user_fd_tables_checked=checked,privileged_read_only_metadata_checks=privileged,gross_allocated_bytes_retired=sum(x['allocated_bytes'] for x in rows),
 scope='Twenty-seven expanded supporting/current CPP/fixture/raw/configure/build ledgers only, SHA/CRC/member/mode/time checked against pushed recovery; all owned births absent and no same-user open FD/current C8/ABBA/all8/frame frozen-input consumer. All original successful owners and concise results preserved. Actual products, sources, owner outcomes and current numerical/cohort inputs remain. Gross retirement, not net drive gain or scientific/release qualification.')
record=B/'RETIRED-EXPANDED-LEDGERS.json';assert not record.exists()
def save():record.write_text(json.dumps(report,indent=2)+'\n')
save()
for row in rows:
 p=R/row['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'];p.unlink();row['removed']=True;save()
report['status']='complete';save();print(json.dumps(dict(retired_files=27,gross_allocated_bytes=report['gross_allocated_bytes_retired'])))
