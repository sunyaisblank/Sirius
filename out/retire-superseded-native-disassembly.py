"""Retire exactly three archived, untracked, unused generated text dumps."""
from pathlib import Path
import hashlib,json,stat,subprocess,zipfile
R=Path.cwd();J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();assert revision==subprocess.check_output(['git','rev-parse','@{upstream}'],text=True).strip()
frame=J(R/'out/native-cpu-direct-compound-cohort/cold-frame/frame/owner.json')
assert frame['status']=='terminal' and frame['source_unchanged'] and frame['whole_inputs_unchanged']
assert not frame['remaining_owned_processes'] and not frame['cleanup_errors']
frozen=set()
for root in ('out/native-cpu-direct-compound','out/native-cpu-direct-compound-cohort','out/native-cpu-direct-compound-compatibility','out/native-cpu-owned-lower-increments','out/native-cpu-frame-route','out/native-cpu-dense-values'):
 for p in (R/root).rglob('actions.json'):
  for action in J(p).values():frozen.update((R/name).resolve() for name in action.get('frozen_inputs',[]))
tracked=set(subprocess.check_output(['git','ls-files'],text=True).splitlines())
rows=[]
for stage,directory in (('native-cpu-static-program','ead88dde/native-cpu-static-program-finite-rejection'),('native-cpu-workgroup-scheduler','cf1ba045/native-cpu-workgroup-scheduler-control3-timeout'),('native-cpu-dispatch-admission','629e3794/native-cpu-dispatch-admission-frame-timeout')):
 name='out/'+stage+'/provider/linked-native.disassembly.txt';p=R/name
 assert p.is_file() and not p.is_symlink() and name not in tracked and p.resolve() not in frozen
 manifest=R/'attestations/qualification'/directory/'MANIFEST.json';m=J(manifest)
 assert manifest.read_bytes()==subprocess.check_output(['git','show','@{upstream}:'+str(manifest.relative_to(R))])
 member=m['members'][name];archive=manifest.parent/member['archive'];identity=m['archives'][member['archive']]
 assert pin(archive)==identity
 assert hashlib.sha256(subprocess.check_output(['git','show','@{upstream}:'+str(archive.relative_to(R))])).hexdigest()==identity['sha256']
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  matches=[k for k in z.namelist() if k.removeprefix('files/')==name];assert len(matches)==1
  key=matches[0];data=z.read(key)
  assert len(data)==member['bytes'] and hashlib.sha256(data).hexdigest()==member['sha256']
  assert stat.S_IMODE(z.getinfo(key).external_attr>>16)==member['mode']
  assert p.read_bytes()==data
 s=p.stat();assert stat.S_IMODE(s.st_mode)==member['mode']
 rows.append({'path':name,'archive':str(archive.relative_to(R)),'archive_identity':identity,'member':key,'identity':pin(p),'mode':stat.S_IMODE(s.st_mode),'mtime_ns':s.st_mtime_ns,'allocated_bytes':s.st_blocks*512,'removed':False})
build=(R/'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja').read_text();assert all(str(R/row['path']) not in build for row in rows)
record=R/'attestations/qualification/9275d0aa/native-cpu-direct-compound-frame-timeout/SUPERSEDED-DISASSEMBLY.json';assert not record.exists()
report={'source':revision,'status':'started','members':rows,'gross_allocated_bytes_retired':sum(v['allocated_bytes'] for v in rows),'scope':'Exactly three regenerable old text dumps. SHA/CRC/member bytes/mode verified against pushed original archives, outside current declared frozen inputs, Git tracked files and Mesa build commands. Historical before/after receipts remain historical and unchanged. Actual source/code/objects/qualification archives and resource caches retained. Old manifests omit nanosecond mtimes; current nonsemantic text-file mtime is recorded here. Gross retired allocation, not net drive gain.'}
def save():record.write_text(json.dumps(report,indent=2)+'\n')
save()
for row in rows:
 p=R/row['path'];assert pin(p)==row['identity'] and p.stat().st_mtime_ns==row['mtime_ns'];p.unlink();row['removed']=True;save()
report['status']='complete';save();print(json.dumps({'removed_text_dumps':len(rows),'gross_allocated_bytes_retired':report['gross_allocated_bytes_retired']}))
