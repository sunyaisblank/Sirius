"""Retire five additional provider pairs with complete recovery metadata."""
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
assert subprocess.check_output(['git', 'rev-parse', '@{upstream}'], text=True).strip() == revision
audit_path = R / 'out/native-cpu-superseded-provider-audit.json'
audit = json.loads(audit_path.read_text())
eligible = [row for row in audit['rows'] if not row['frozen'] and row['recovery']
            and row['path'] != 'bin/linux-gcc/native-cpu-simd-packet-policy/provider/libvulkan_lvp.so']
assert len(eligible) == 5
assert json.loads((R / 'out/native-cpu-product-inline-cohort/control-8-acceptance.json').read_text())['passed']
frozen = set()
for root in ('out/native-cpu-product-inline', 'out/native-cpu-product-inline-cohort', 'out/native-cpu-product-inline-compatibility'):
    for p in (R / root).rglob('actions.json'):
        for action in json.loads(p.read_text()).values():
            frozen.update((R / name).resolve() for name in action['frozen_inputs'])
tracked = set(subprocess.check_output(['git', 'ls-files'], text=True).splitlines())
members = []
archives = []
for candidate in eligible:
    recovery, = candidate['recovery']
    manifest = R / recovery['manifest']
    metadata = json.loads(manifest.read_text())['archives'][recovery['archive']]
    archive = manifest.parent / recovery['archive']
    assert archive.stat().st_size == metadata['bytes']
    assert hashlib.file_digest(archive.open('rb'), 'sha256').hexdigest() == metadata['sha256']
    assert hashlib.sha256(subprocess.check_output(['git', 'show', '@{upstream}:' + str(archive.relative_to(R))])).hexdigest() == metadata['sha256']
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None and set(z.namelist()) == set(metadata['members'])
        for key, v in metadata['members'].items():
            data = z.read(key)
            assert len(data) == v['bytes'] and hashlib.sha256(data).hexdigest() == v['sha256']
            assert stat.S_IMODE(z.getinfo(key).external_attr >> 16) == v['mode']
        for suffix in ('libvulkan_lvp.so', 'lvp_icd.json'):
            key = str(Path(candidate['path']).parent / suffix)
            v = metadata['members'][key]; p = R / key; s = p.stat()
            assert p.is_file() and not p.is_symlink() and p.resolve() not in frozen and key not in tracked
            assert s.st_size == v['bytes'] and stat.S_IMODE(s.st_mode) == v['mode'] and s.st_mtime_ns == v['mtime_ns']
            assert p.read_bytes() == z.read(key)
            members.append({'path': key, 'archive': str(archive.relative_to(R)), 'member': key,
                'bytes': s.st_size, 'allocated_bytes': s.st_blocks * 512, 'sha256': v['sha256'],
                'mode': v['mode'], 'mtime_ns': v['mtime_ns'], 'removed': False})
    archives.append({'archive': str(archive.relative_to(R)), 'sha256': metadata['sha256'],
                     'verified_product_members': len(metadata['members'])})
assert len(members) == 10
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
record = R / 'attestations/qualification/f69d74df/native-cpu-lazy-nir-frame-timeout/OLDER-SUPERSEDED-PROVIDERS.json'
assert not record.exists()
report = {'schema': 'sirius-older-superseded-provider-retirement-v1', 'source': revision,
    'status': 'started', 'members': members, 'archives': archives,
    'gross_allocated_bytes_retired': sum(row['allocated_bytes'] for row in members),
    'gross_logical_bytes_retired': sum(row['bytes'] for row in members),
    'same_user_process_maps_checked': mapped, 'privileged_read_only_maps_checks': privileged,
    'scope': 'Ten exact archived untracked provider names only, outside all current W frozen closures/build commands and observed process maps. Five full product ZIPs independently verified including SHA/CRC/modes and pushed Git bytes. I provider pair retained because its old archive metadata omits mode/mtime. Native objects, compiler/source/build caches, unique unarchived providers, current W and still-frozen NIR provider retained. Gross retired allocation is not net filesystem gain.'}
def save():
    record.write_text(json.dumps(report, indent=2) + '\n')
save()
for row in members:
    p = R / row['path']
    assert hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == row['sha256']
    p.unlink(); row['removed'] = True; save()
for p in {(R / row['path']).parent for row in members}:
    p.rmdir()
assert all(not (R / row['path']).exists() for row in members)
report['status'] = 'complete'; save()
print(json.dumps({'removed_files': 10, 'gross_allocated_bytes_retired': report['gross_allocated_bytes_retired']}))
