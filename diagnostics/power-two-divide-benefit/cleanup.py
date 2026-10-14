"""Retire only two fully preserved matched-benefit roots after a fresh scoped audit."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/power-two-divide-benefit'
archive = root / 'attestations/software-vulkan/4bcd987/power-two-divide-benefit'
revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
assert work == root / 'out/power-two-divide-benefit'
load = lambda path: json.loads(path.read_text())


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


manifest_path = archive / 'whole-byte-manifest.json'
manifest = load(manifest_path)
expected = manifest['files']
actual = {str(path.relative_to(archive / 'files')): path
    for path in (archive / 'files').rglob('*') if path.is_file()}
assert set(actual) == set(expected)
for name, path in actual.items():
    assert not path.is_symlink() and identity(path) == expected[name]
    assert oct(path.stat().st_mode & 0o777) == manifest['preserved_file_modes'][name]
assert {str(path.relative_to(archive)) for path in archive.rglob('*') if path.is_file()} == {
    'whole-byte-manifest.json', *('files/' + name for name in expected)}
archived_directories = {str(path.relative_to(archive / 'files')): oct(path.stat().st_mode & 0o777)
    for name in manifest['roots']
    for path in [archive / 'files' / name, *[p for p in (archive / 'files' / name).rglob('*') if p.is_dir()]]}
assert archived_directories == manifest['preserved_directories']


def check_live_files():
    directories = {}
    for owned in [work, stage]:
        assert owned.is_dir() and not owned.is_symlink()
        for path in owned.rglob('*'):
            assert not path.is_symlink() and (path.is_file() or path.is_dir()), path
        for path in [owned, *[p for p in owned.rglob('*') if p.is_dir()]]:
            name = str(path.relative_to(root))
            saved = path.stat()
            directories[name] = oct(saved.st_mode & 0o777)
            assert {'device': saved.st_dev, 'inode': saved.st_ino} == manifest['live_directory_identities'][name]
    assert directories == manifest['preserved_directories']
    live = {str(path.relative_to(root)): path for owned in [work, stage]
        for path in owned.rglob('*') if path.is_file()}
    assert set(live) == set(expected)
    for name, path in live.items():
        assert not path.is_symlink() and identity(path) == expected[name]
        assert oct(path.stat().st_mode & 0o777) == manifest['preserved_file_modes'][name]
    return live


check_live_files()
sources = load(work / 'preserved-source.json')
assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + sources['tag']], text=True).strip() == sources['remote_readback']
for record in sources['source_records']:
    saved = subprocess.check_output(['git', 'show', sources['commit'] + ':' + record['path']])
    assert saved == (root / record['origin']).read_bytes()


def proc(pid):
    words = Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()
    return {'pid': pid, 'ppid': int(words[1]), 'start_ticks': words[19]}


boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
births = set()
receipts = sorted(work.rglob('owner.json'))
assert len(receipts) == 4
for path in receipts:
    owner = load(path)
    assert owner['status'] == 'terminal' and owner['passed'] and not owner['remaining_owned_processes']
    assert not owner['cleanup_errors']
    for item in owner['observed_births']:
        assert item['absent']
        births.add((owner['boot_id'], int(item['pid']), str(item['saved_start_ticks'])))
    births.add((owner['boot_id'], int(owner['controller']['pid']), str(owner['controller']['start_ticks'])))
for token in births:
    if token[0] != boot:
        continue
    try:
        assert proc(token[1])['start_ticks'] != token[2], token
    except (FileNotFoundError, ProcessLookupError):
        pass
excluded = set()
pid = os.getpid()
while pid > 0 and pid not in excluded:
    excluded.add(pid)
    try:
        pid = proc(pid)['ppid']
    except (FileNotFoundError, ProcessLookupError):
        break
roots = [str(work), str(stage)]
consumers, permission_limited = [], set()
for entry in Path('/proc').glob('[0-9]*'):
    pid = int(entry.name)
    if pid in excluded:
        continue
    try:
        before = proc(pid)
    except PermissionError:
        permission_limited.add(pid)
        continue
    except (FileNotFoundError, ProcessLookupError):
        continue
    strings = []
    for field in ['cwd', 'cmdline', 'maps']:
        try:
            if field == 'cwd':
                strings.append(os.readlink(entry / field))
            elif field == 'cmdline':
                strings += [os.fsdecode(value) for value in (entry / field).read_bytes().split(b'\0') if value]
            else:
                strings.append((entry / field).read_text())
        except PermissionError:
            permission_limited.add(pid)
        except (FileNotFoundError, ProcessLookupError):
            pass
    try:
        descriptors = list((entry / 'fd').iterdir())
    except PermissionError:
        permission_limited.add(pid)
        descriptors = []
    except (FileNotFoundError, ProcessLookupError):
        descriptors = []
    for fd in descriptors:
        try:
            strings.append(os.readlink(fd))
        except PermissionError:
            permission_limited.add(pid)
        except (FileNotFoundError, ProcessLookupError):
            pass
    if any(owned in value for owned in roots for value in strings):
        try:
            assert proc(pid)['start_ticks'] == before['start_ticks']
            consumers.append(before)
        except (FileNotFoundError, ProcessLookupError):
            pass
assert not consumers, consumers
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert subprocess.check_output(['git', 'worktree', 'list', '--porcelain'], text=True).count('worktree ') == 1
for name, pin in load(work / 'protected-before.json').items():
    assert identity(root / name) == {key: pin[key] for key in ['bytes', 'sha256']}
live = check_live_files()
report = {'roots': roots, 'removed_files': len(live),
    'removed_bytes': sum(path.stat().st_size for path in live.values()),
    'sealed_manifest_sha256': identity(manifest_path)['sha256'],
    'saved_births_absent': len(births), 'accessible_consumers': consumers,
    'permission_limited_processes': len(permission_limited),
    'source_revision': revision, 'evidence_tag': sources['tag'], 'evidence_commit': sources['commit'],
    'scope': 'Exact source/archive recovery and fresh saved-birth plus accessible cwd/argv/maps/fd checks. Observer/ancestors excluded. No global hidden-process or scheduled-consumer absence claim.'}
destination = root / 'attestations/cleanup/2026-10-14/power-two-divide-benefit/cleanup.json'
assert not destination.exists()
destination.parent.mkdir(parents=True, exist_ok=True)
shutil.rmtree(work)
shutil.rmtree(stage)
destination.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
