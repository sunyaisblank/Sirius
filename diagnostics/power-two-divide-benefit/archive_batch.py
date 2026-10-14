"""Retain exact completed original matched benefit evidence; leave live roots for cleanup review."""
from pathlib import Path
import hashlib
import json
import shutil
import stat
import subprocess
import sys

sys.dont_write_bytecode = True
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/power-two-divide-benefit'
archive = root / 'attestations/software-vulkan/4bcd987/power-two-divide-benefit'
revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
assert sys.argv[1:] == ['--archive'] and not archive.exists()
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not subprocess.check_output(['git', 'status', '--porcelain'])


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


receipts = sorted(work.rglob('owner.json'))
assert len(receipts) == 4
for path in receipts:
    owner = json.loads(path.read_text())
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(item['absent'] for item in owner['observed_births'])
    assert owner['source'] == {'revision': revision, 'status': ''}
    assert owner['whole_inputs_unchanged'] and owner['source_unchanged']
for name, expected in json.loads((work / 'protected-before.json').read_text()).items():
    assert identity(root / name) == {key: expected[key] for key in ['bytes', 'sha256']}, name
disposition = json.loads((work / 'disposition.json').read_text())
assert disposition['decision'] in ['REJECT_POWER_TWO_DIVIDE_COST', 'PROCEED_TO_PRODUCTION_IMPLEMENTATION']
assert not disposition['production_adopted'] and disposition['original_benefit_trial_completed']
sources = json.loads((work / 'preserved-source.json').read_text())
assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + sources['tag']], text=True).strip() == sources['remote_readback']
assert not (work / 'source.index').exists()
archive.mkdir(parents=True)
files, modes = {}, {}
directories, directory_identities = {}, {}
for owned in [work, stage]:
    assert owned.is_dir() and not owned.is_symlink()
    for path in [owned, *sorted(owned.rglob('*'))]:
        mode = path.lstat().st_mode
        assert not stat.S_ISLNK(mode) and (stat.S_ISREG(mode) or stat.S_ISDIR(mode)), path
        if stat.S_ISDIR(mode):
            name = str(path.relative_to(root))
            target = archive / 'files' / name
            target.mkdir(parents=True, exist_ok=True)
            target.chmod(mode & 0o777)
            directories[name] = oct(mode & 0o777)
            saved = path.stat()
            directory_identities[name] = {'device': saved.st_dev, 'inode': saved.st_ino}
            continue
        name = str(path.relative_to(root))
        target = archive / 'files' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        expected = identity(path)
        shutil.copy2(path, target)
        assert identity(path) == expected and identity(target) == expected, name
        files[name] = expected
        modes[name] = oct(mode & 0o777)
        assert oct(target.stat().st_mode & 0o777) == modes[name]
assert {str(path.relative_to(root)) for owned in [work, stage]
    for path in owned.rglob('*') if path.is_file()} == set(files)
current_directories = {str(path.relative_to(root)): path for owned in [work, stage]
    for path in [owned, *[p for p in owned.rglob('*') if p.is_dir()]]}
assert set(current_directories) == set(directories)
for name, path in current_directories.items():
    saved = path.stat()
    assert {'device': saved.st_dev, 'inode': saved.st_ino} == directory_identities[name]
    assert oct(saved.st_mode & 0o777) == directories[name]
manifest = {'source_revision': revision, 'evidence_tag': sources['tag'],
    'evidence_commit': sources['commit'],
    'roots': [str(path.relative_to(root)) for path in [work, stage]],
    'files': files, 'preserved_file_modes': modes, 'file_count': len(files),
    'preserved_directories': directories, 'live_directory_identities': directory_identities,
    'total_bytes': sum(record['bytes'] for record in files.values()),
    'terminal_owner_receipts': [str(path.relative_to(root)) for path in receipts],
    'protected_payloads_unchanged': 114,
    'scope': 'Exact local raw single original20+2matched trial recovery. Passed or rejected fixed benefit verdict retained; no full qualification, production adoption or remote binary backup claim.'}
(archive / 'whole-byte-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({'archive': str(archive.relative_to(root)), 'files': len(files),
    'bytes': manifest['total_bytes'], 'terminal_owners': len(receipts)}))
