"""Retire Y's exact archived provider pair after Z replaces its active build."""
from pathlib import Path
import hashlib
import json
import os
import stat
import subprocess
import zipfile

R = Path.cwd()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert revision == subprocess.check_output(['git', 'rev-parse', '@{upstream}'], text=True).strip()
assert json.loads((R / 'out/native-cpu-direct-scratch-cohort/control-8-acceptance.json').read_text())['passed']
frozen = { (R / name).resolve() for root in ('out/native-cpu-direct-scratch', 'out/native-cpu-direct-scratch-cohort', 'out/native-cpu-direct-scratch-compatibility')
    for p in (R / root).rglob('actions.json') for action in json.loads(p.read_text()).values() for name in action['frozen_inputs'] }
manifest = R / 'attestations/qualification/bee63496/native-cpu-task-cardinality-frame-timeout/MANIFEST.json'
metadata, = [v for v in json.loads(manifest.read_text())['archives'] if v['file'] == 'actual-products.zip']
archive = manifest.parent / metadata['file']
assert hashlib.file_digest(archive.open('rb'), 'sha256').hexdigest() == metadata['identity']['sha256']
assert hashlib.sha256(subprocess.check_output(['git', 'show', '@{upstream}:' + str(archive.relative_to(R))])).hexdigest() == metadata['identity']['sha256']
tracked = set(subprocess.check_output(['git', 'ls-files'], text=True).splitlines())
members = []
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None and set(z.namelist()) == set(metadata['members'])
    for key, v in metadata['members'].items():
        data = z.read(key)
        assert len(data) == v['bytes'] and hashlib.sha256(data).hexdigest() == v['sha256']
        assert stat.S_IMODE(z.getinfo(key).external_attr >> 16) == v['mode']
    for suffix in ('libvulkan_lvp.so', 'lvp_icd.json'):
        p = R / 'bin/linux-gcc/native-cpu-task-cardinality/provider' / suffix
        key = 'files/' + str(p.relative_to(R)); v = metadata['members'][key]; s = p.stat()
        assert p.is_file() and not p.is_symlink() and p.resolve() not in frozen
        assert str(p.relative_to(R)) not in tracked
        assert s.st_size == v['bytes'] and stat.S_IMODE(s.st_mode) == v['mode'] and s.st_mtime_ns == v['mtime_ns']
        assert p.read_bytes() == z.read(key)
        members.append({'path': str(p.relative_to(R)), 'member': key, 'sha256': v['sha256'],
            'bytes': s.st_size, 'allocated_bytes': s.st_blocks * 512, 'mode': v['mode'],
            'mtime_ns': v['mtime_ns'], 'removed': False})
paths = {str(R / row['path']) for row in members}
mapped = 0; privileged = []
for process in Path('/proc').iterdir():
    if not process.name.isdecimal():
        continue
    try:
        if process.stat().st_uid != os.getuid():
            continue
        try:
            content = (process / 'maps').read_text()
        except PermissionError:
            content = subprocess.check_output(['sudo', '-n', '/usr/bin/cat', str(process / 'maps')], text=True)
            privileged.append(int(process.name))
    except FileNotFoundError:
        continue
    assert all(path not in content for path in paths), process.name
    mapped += 1
build = (R / 'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja').read_text()
assert all(path not in build for path in paths)
record = manifest.parent / 'SUPERSEDED-PROVIDER.json'
assert not record.exists()
report = {'source': revision, 'status': 'started', 'members': members,
    'archive': str(archive.relative_to(R)), 'archive_sha256': metadata['identity']['sha256'],
    'verified_product_members': len(metadata['members']), 'same_user_maps_checked': mapped,
    'privileged_read_only_maps_checks': privileged,
    'gross_allocated_bytes_retired': sum(v['allocated_bytes'] for v in members),
    'scope': 'Two untracked Y provider names only, independently SHA/CRC/mode verified against both pushed archive products, excluded from all current Z frozen inputs/build commands/process maps. Y source/evidence and reused X raw3/diagnostic plus active Z/native objects/caches/resources retained. Gross retired allocation, not net drive gain.'}
def save():
    record.write_text(json.dumps(report, indent=2) + '\n')
save()
for row in members:
    p = R / row['path']; assert hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == row['sha256']
    p.unlink(); row['removed'] = True; save()
(R / 'bin/linux-gcc/native-cpu-task-cardinality/provider').rmdir()
assert all(not (R / row['path']).exists() for row in members)
report['status'] = 'complete'; save()
print(json.dumps({'removed_names': 2, 'gross_allocated_bytes_retired': report['gross_allocated_bytes_retired']}))
