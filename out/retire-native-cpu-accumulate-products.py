"""Retire six superseded R products while keeping its active driver/baseline."""
from pathlib import Path
import hashlib
import json
import os
import stat
import subprocess
import zipfile

R=Path.cwd()
A=R/'attestations/qualification/603a155e/native-cpu-accumulate-inline-frame-timeout'
record_path=A/'SUPERSEDED-PRODUCTS.json'
assert not record_path.exists()
assert not subprocess.check_output(['git','status','--porcelain']).strip()
manifest=json.loads((A/'MANIFEST.json').read_text())
bundle=next(row for row in manifest['archives'] if row['file']=='actual-products.zip')
archive=A/bundle['file']
relative=str(archive.relative_to(R))
upstream='origin/development/native-cpu-slang-bridge-2026-10-23'
assert subprocess.check_output(['git','rev-parse',upstream+':'+relative]).strip()==subprocess.check_output(['git','hash-object',relative]).strip()
assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==bundle['identity']['sha256']
selected=[]
with zipfile.ZipFile(archive) as zipped:
    assert zipped.testzip() is None and set(zipped.namelist())==set(bundle['members'])
    for name,saved in bundle['members'].items():
        data=zipped.read(name)
        assert len(data)==saved['bytes'] and hashlib.sha256(data).hexdigest()==saved['sha256']
        path=Path(saved['source'])
        if not (path.is_relative_to(R/'bin/linux-gcc/native-cpu-accumulate-inline') or
                path==R/'bin/linux-gcc/native-cpu-accumulate-compatibility/compatibility'):continue
        assert not path.is_symlink() and path.resolve()==path
        s=path.stat()
        assert s.st_nlink==1 and s.st_size==saved['bytes'] and s.st_mtime_ns==saved['mtime_ns']
        assert stat.S_IMODE(s.st_mode)==saved['mode']
        assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==saved['sha256']
        selected.append({'path':str(path.relative_to(R)),'member':name,'allocated_bytes':s.st_blocks*512,
                         **{k:saved[k] for k in ('bytes','sha256','mode','mtime_ns')}})
assert len(selected)==6
tracked={os.fsdecode(name) for name in subprocess.check_output(['git','ls-files','-z']).split(b'\0') if name}
assert not tracked.intersection(row['path'] for row in selected)
names={str(R/row['path']) for row in selected}
current=R/'out/native-cpu-scalar-zero-cohort'
actions=[R/'out/native-cpu-scalar-zero/compile/actions.json',R/'out/native-cpu-scalar-zero/provider/actions.json',
    R/'out/native-cpu-scalar-zero-compatibility/actions.json',current/'actions.json',
    current/'benefit/runtime/actions.json',current/'controls/runtime/actions.json']
for path in actions:
    for action in json.loads(path.read_text()).values():
        assert not names.intersection(str((R/name).resolve()) for name in action['frozen_inputs'])
for path in (R/'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja',
             R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/targets/lavapipe/meson.build'):
    text=path.read_text();assert not any(name in text for name in names)
for process in Path('/proc').glob('[0-9]*'):
    try:text=(process/'maps').read_text()
    except (PermissionError,FileNotFoundError,ProcessLookupError):continue
    assert not any(name in text for name in names),process
owners=list(current.rglob('owner.json'))
assert len(owners)==8  # C8, four ABBA slots and original controls1/2/3.
for path in owners:
    owner=json.loads(path.read_text())
    assert owner['status']=='terminal' and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert json.loads((current/'controls/runtime/control-3/owner.json').read_text())['stop_reason']=='diagnostic_timeout'
before=os.statvfs(R)
record={'status':'started','recovery_archive':relative,'recovery_identity':bundle['identity'],
    'files':selected,'retired_files':[],'allocated_bytes':sum(row['allocated_bytes'] for row in selected),
    'logical_bytes':sum(row['bytes'] for row in selected),'free_bytes_before':before.f_bavail*before.f_frsize,
    'scope':'Only six exact archived superseded R object/provider/diagnostic products. The R driver and renamed Q baseline, current zlib, source/generators, S/U objects, active build cache, original resources and unrelated work remain. All active U closures/build commands exclude retired paths; historical product paths recover from recorded ZIP member/mode/mtime.'}


def save():record_path.write_text(json.dumps(record,indent=2)+'\n')


save()
for row in selected:
    path=R/row['path'];assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==row['sha256']
    path.unlink();record['retired_files'].append(row['path']);save()
for folder in (R/'bin/linux-gcc/native-cpu-accumulate-inline/objects',R/'bin/linux-gcc/native-cpu-accumulate-inline/provider'):
    folder.rmdir()
after=os.statvfs(R)
record.update(status='complete',free_bytes_after=after.f_bavail*after.f_frsize,
    free_bytes_observed_delta=(after.f_bavail*after.f_frsize)-(before.f_bavail*before.f_frsize),
    free_space_attribution='Observed filesystem delta includes concurrent unrelated writes; gross removed allocated bytes measured separately.')
save()
assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==bundle['identity']['sha256']
print(json.dumps({k:record[k] for k in ('status','logical_bytes','allocated_bytes','free_bytes_observed_delta')}))
