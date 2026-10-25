"""Retire archived raw W trial output outside every current Z/contingent AA input closure."""
from pathlib import Path
import hashlib
import json
import stat
import subprocess
import zipfile

R = Path.cwd(); C = R / 'out/native-cpu-product-inline-cohort'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert revision == subprocess.check_output(['git', 'rev-parse', '@{upstream}'], text=True).strip()
assert json.loads((R / 'out/native-cpu-direct-scratch-cohort/control-8-acceptance.json').read_text())['passed']
frozen = {(R / name).resolve() for root in ('out/native-cpu-direct-scratch', 'out/native-cpu-direct-scratch-cohort', 'out/native-cpu-direct-scratch-compatibility', 'out/native-cpu-static-registers', 'out/native-cpu-static-registers-cohort', 'out/native-cpu-static-registers-compatibility')
          for p in (R / root).rglob('actions.json') for action in json.loads(p.read_text()).values() for name in action['frozen_inputs']}
tracked = set(subprocess.check_output(['git', 'ls-files'], text=True).splitlines())
manifest = R / 'attestations/qualification/4595f5a6/native-cpu-product-inline-frame-timeout/MANIFEST.json'
saved, = [v for v in json.loads(manifest.read_text())['archives'] if v['file'] == 'source-build-and-trial.zip']
archive = manifest.parent / saved['file']
assert hashlib.file_digest(archive.open('rb'), 'sha256').hexdigest() == saved['identity']['sha256']
assert hashlib.sha256(subprocess.check_output(['git', 'show', '@{upstream}:' + str(archive.relative_to(R))])).hexdigest() == saved['identity']['sha256']
members = []
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None and set(z.namelist()) == set(saved['members'])
    for key, v in saved['members'].items():
        data = z.read(key)
        assert len(data) == v['bytes'] and hashlib.sha256(data).hexdigest() == v['sha256']
        assert stat.S_IMODE(z.getinfo(key).external_attr >> 16) == v['mode']
        p = Path(v['source'])
        if not p.is_relative_to(C) or str(p.relative_to(R)) in tracked or p.resolve() in frozen:
            continue
        assert p.is_file() and not p.is_symlink() and p.resolve().is_relative_to(C)
        s = p.stat(); assert s.st_size == v['bytes'] and s.st_mtime_ns == v['mtime_ns'] and stat.S_IMODE(s.st_mode) == v['mode']
        assert p.read_bytes() == data
        members.append({'path': str(p.relative_to(R)), 'member': key, 'bytes': s.st_size,
            'allocated_bytes': s.st_blocks * 512, 'sha256': v['sha256'], 'mode': v['mode'], 'mtime_ns': v['mtime_ns'], 'removed': False})
assert members
# Read-only inode consumer check, including protected service processes.
result = subprocess.run(['sudo', '-n', '/usr/bin/fuser', '-s', *[str(R / v['path']) for v in members]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
assert result.returncode == 1 and not result.stdout and not result.stderr, (result.returncode, result.stdout, result.stderr)
build = (R / 'bin/linux-gcc/lavapipe-acqrel-ordering/build.ninja').read_text()
assert all(str(R / v['path']) not in build for v in members)
record = manifest.parent / 'RAW-COHORT-RETIREMENT.json'; assert not record.exists()
report = {'source': revision, 'status': 'started', 'members': members,
    'archive': str(archive.relative_to(R)), 'archive_sha256': saved['identity']['sha256'],
    'verified_archive_members': len(saved['members']), 'open_inode_check': 'privileged read-only fuser: no consumers',
    'gross_allocated_bytes_retired': sum(v['allocated_bytes'] for v in members),
    'gross_logical_bytes_retired': sum(v['bytes'] for v in members),
    'scope': 'Only exact archived untracked raw W cohort files outside every current Z/contingent AA frozen closure and build command. Whole source ZIP independently CRC/SHA/mode verified and matched against pushed Git bytes. Tracked source/findings, W source/native objects/full38 standalone component and all current Z/contingent AA inputs/caches/resources preserved. Recovery by explicit archive member/mode/mtime. Gross retired allocation is not net drive gain.'}
def save():
    record.write_text(json.dumps(report, indent=2) + '\n')
save()
for v in members:
    p = R / v['path']; assert hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == v['sha256']
    p.unlink(); v['removed'] = True; save()
for p in sorted((p for p in C.rglob('*') if p.is_dir()), key=lambda p: len(p.parts), reverse=True):
    if not any(p.iterdir()):
        p.rmdir()
assert all(not (R / v['path']).exists() for v in members)
report['status'] = 'complete'; save()
print(json.dumps({'retired_raw_files': len(members), 'gross_allocated_bytes_retired': report['gross_allocated_bytes_retired']}))
