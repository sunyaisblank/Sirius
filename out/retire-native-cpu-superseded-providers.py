"""Retire six independently archived provider files after W replaces them."""
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
new = R / 'out/native-cpu-product-inline'
assert json.loads((new / 'provider/linked-native-proof.json').read_text())['passed']
assert json.loads((R / 'out/native-cpu-product-inline-cohort/control-8-acceptance.json').read_text())['passed']
frozen = set()
for p in (new, R / 'out/native-cpu-product-inline-cohort', R / 'out/native-cpu-product-inline-compatibility'):
    for f in p.rglob('actions.json'):
        for action in json.loads(f.read_text()).values():
            frozen.update((R / name).resolve() for name in action['frozen_inputs'])
tracked = set(subprocess.check_output(['git', 'ls-files'], text=True).splitlines())
members = []
archives = []
for name, recovery in (
    ('zero-term', '115b86f3/native-cpu-zero-term-frame-timeout'),
    ('scalar-zero', 'fb50bb2a/native-cpu-scalar-zero-control3-timeout'),
    ('lazy-nir', 'f69d74df/native-cpu-lazy-nir-frame-timeout')):
    manifest = R / 'attestations/qualification' / recovery / 'MANIFEST.json'
    saved = json.loads(manifest.read_text())
    for entry in saved['archives']:
        if entry['file'] != 'actual-products.zip':
            continue
        archive = manifest.parent / entry['file']
        assert hashlib.file_digest(archive.open('rb'), 'sha256').hexdigest() == entry['identity']['sha256']
        upstream = subprocess.check_output(['git', 'show', '@{upstream}:' + str(archive.relative_to(R))])
        assert hashlib.sha256(upstream).hexdigest() == entry['identity']['sha256']
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None and set(z.namelist()) == set(entry['members'])
            for key, v in entry['members'].items():
                data = z.read(key)
                assert len(data) == v['bytes'] and hashlib.sha256(data).hexdigest() == v['sha256']
                assert stat.S_IMODE(z.getinfo(key).external_attr >> 16) == v['mode']
            for suffix in ('libvulkan_lvp.so', 'lvp_icd.json'):
                p = R / ('bin/linux-gcc/native-cpu-' + name + '/provider') / suffix
                key = 'files/' + str(p.relative_to(R)); v = entry['members'][key]
                assert p.is_file() and not p.is_symlink() and p.resolve() not in frozen
                assert str(p.relative_to(R)) not in tracked
                s = p.stat()
                assert s.st_size == v['bytes'] and stat.S_IMODE(s.st_mode) == v['mode'] and s.st_mtime_ns == v['mtime_ns']
                assert p.read_bytes() == z.read(key)
                members.append({'path': str(p.relative_to(R)), 'archive': str(archive.relative_to(R)),
                    'member': key, 'bytes': s.st_size, 'allocated_bytes': s.st_blocks * 512,
                    'sha256': v['sha256'], 'mode': v['mode'], 'mtime_ns': v['mtime_ns'], 'removed': False})
        archives.append({'archive': str(archive.relative_to(R)), 'sha256': entry['identity']['sha256'],
                         'independently_verified_product_members': len(entry['members'])})
assert len(members) == 6
paths = {str(R / row['path']) for row in members}
mapped = []
privileged_maps = []
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
            privileged_maps.append(int(process.name))
    except FileNotFoundError:
        continue
    for path in paths:
        assert path not in content, (process.name, path)
    mapped.append(int(process.name))
build = (R / 'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja').read_text()
assert all(path not in build for path in paths)
record = R / 'attestations/qualification/f69d74df/native-cpu-lazy-nir-frame-timeout/SUPERSEDED-PROVIDERS.json'
assert not record.exists()
report = {'schema': 'sirius-superseded-provider-retirement-v1', 'source': revision,
    'status': 'started', 'members': members, 'recovery': archives,
    'gross_allocated_bytes_retired': sum(v['allocated_bytes'] for v in members),
    'gross_logical_bytes_retired': sum(v['bytes'] for v in members),
    'checked_same_user_maps': len(mapped),
    'privileged_read_only_maps_checks': privileged_maps,
    'scope': 'Six exact archived untracked S/U/V provider files only. W provider/current frozen inputs/native objects/source generators/active build cache/original resources unchanged. Product archives were independently CRC/SHA/mode verified and checked against pushed Git bytes before unlink. Gross retired allocation is not net project or filesystem space gain.'}
def save():
    record.write_text(json.dumps(report, indent=2) + '\n')
save()
for row in members:
    p = R / row['path']
    assert hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == row['sha256']
    p.unlink(); row['removed'] = True; save()
for directory in { (R / row['path']).parent for row in members }:
    directory.rmdir()
assert all(not (R / row['path']).exists() for row in members)
report['status'] = 'complete'; save()
print(json.dumps({'removed_files': 6, 'gross_allocated_bytes_retired': report['gross_allocated_bytes_retired']}))
