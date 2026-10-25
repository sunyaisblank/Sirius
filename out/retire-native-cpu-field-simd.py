"""Retire only archived rejected-T files after the active Mesa cache moved to U."""
from pathlib import Path
import hashlib
import json
import os
import stat
import subprocess
import zipfile

R=Path.cwd()
A=R/'attestations/qualification/ae2c2ae6/native-cpu-field-simd-benefit-timeout'
assert not (A/'CLEANUP.json').exists()
assert not subprocess.check_output(['git','status','--porcelain']).strip()
subprocess.run(['git','merge-base','--is-ancestor','d92215eb','origin/development/native-cpu-slang-bridge-2026-10-23'],check=True)
manifest=json.loads((A/'MANIFEST.json').read_text())
archives={}
for bundle in manifest['archives']:
    path=A/bundle['file']
    assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==bundle['identity']['sha256']
    with zipfile.ZipFile(path) as zipped:
        assert zipped.testzip() is None and set(zipped.namelist())==set(bundle['members'])
        for name,row in bundle['members'].items():
            data=zipped.read(name)
            assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256']
            archives[row['source']]={**row,'archive':bundle['file'],'member':name}
roots=[R/'bin/linux-gcc/native-cpu-field-simd',R/'bin/linux-gcc/native-cpu-field-compatibility',
       R/'out/native-cpu-field-simd',R/'out/native-cpu-field-cohort',R/'out/native-cpu-field-compatibility']
tracked={os.fsdecode(name) for name in subprocess.check_output(['git','ls-files','-z']).split(b'\0') if name}
files=[]
for root in roots:
    assert root.resolve()==root and not root.is_symlink()
    for path in root.rglob('*'):
        if not path.is_file() or str(path.relative_to(R)) in tracked:continue
        assert not path.is_symlink() and path.resolve().is_relative_to(root)
        saved=archives[str(path)];s=path.stat()
        assert s.st_nlink==1 and s.st_size==saved['bytes'] and s.st_mtime_ns==saved['mtime_ns']
        assert stat.S_IMODE(s.st_mode)==saved['mode']
        assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==saved['sha256']
        files.append({'path':str(path.relative_to(R)),'allocated_bytes':s.st_blocks*512,
                      **{k:saved[k] for k in ('bytes','sha256','mode','mtime_ns','archive','member')}})
assert len(files)==144
# Only these failed prototype files disappear; every input of the active U
# compile/component/provider actions, including the useful fixed ELF reader,
# must remain outside the retirement set.
retired={row['path'] for row in files}
for path in (R/'out/native-cpu-scalar-zero/compile/actions.json',
             R/'out/native-cpu-scalar-zero/provider/actions.json',
             R/'out/native-cpu-scalar-zero-compatibility/actions.json'):
    for action in json.loads(path.read_text()).values():
        assert not retired.intersection(str((R/name).relative_to(R)) for name in action['frozen_inputs'] if (R/name).is_relative_to(R))
for path in (R/'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja',
             R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/targets/lavapipe/meson.build'):
    text=path.read_text()
    assert not any(str(R/name) in text for name in retired)
provider=R/'out/native-cpu-scalar-zero/provider'
assert json.loads((provider/'build/owner.json').read_text())['passed']
assert json.loads((provider/'linked-native-proof.json').read_text())['passed']
readable_maps=0
for process in Path('/proc').glob('[0-9]*'):
    try:
        text=(process/'maps').read_text();readable_maps+=1
    except (PermissionError,FileNotFoundError,ProcessLookupError):continue
    assert not any(str(root) in text for root in roots[:2]),process
before=os.statvfs(R)
record={'status':'started','recovery_commit':'d92215eb','files':files,'retired_files':[],
    'logical_bytes':sum(row['bytes'] for row in files),'allocated_bytes':sum(row['allocated_bytes'] for row in files),
    'free_bytes_before':before.f_bavail*before.f_frsize,'readable_process_maps_checked':readable_maps,
    'scope':'Only exact ZIP-verified untracked files of rejected field-SIMD T, after committed/pushed recovery and active Mesa/U input/command migration. All tracked source, useful fixed ELF reader, U inputs, baseline/driver, original resources, active build cache and unrelated work preserved. Archives are immutable; originals recover by recorded member/path/mode/mtime.'}


def save():
    (A/'CLEANUP.json').write_text(json.dumps(record,indent=2)+'\n')


save()
for row in files:
    path=R/row['path']
    assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==row['sha256']
    path.unlink();record['retired_files'].append(row['path']);save()
directories={root for root in roots}
for root in roots:directories.update(path for path in root.rglob('*') if path.is_dir())
removed=[]
for path in sorted(directories,key=lambda path:len(path.parts),reverse=True):
    if not any(path.iterdir()):path.rmdir();removed.append(str(path.relative_to(R)))
after=os.statvfs(R)
record.update(status='complete',removed_empty_directories=removed,free_bytes_after=after.f_bavail*after.f_frsize,
    free_bytes_observed_delta=(after.f_bavail*after.f_frsize)-(before.f_bavail*before.f_frsize),
    net_free_space_attribution='Observed filesystem delta includes concurrent unrelated writes; gross removed allocated file bytes are separately measured. Archive creation predates this measurement.')
save()
for bundle in manifest['archives']:
    path=A/bundle['file'];assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==bundle['identity']['sha256']
print(json.dumps({k:record[k] for k in ('status','logical_bytes','allocated_bytes','free_bytes_observed_delta')}))
