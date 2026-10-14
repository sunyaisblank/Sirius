"""Preserve completed issue99 prerequisites without changing production HEAD."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import stat
import subprocess
import sys

sys.dont_write_bytecode = True
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
tag = 'evidence/issue-99-power-two-divide-prerequisites-4bcd987'
archive = root / 'attestations/software-vulkan/4bcd987/power-two-divide-prerequisites'
names = ['screen', 'arithmetic', 'controls']
roots = [root / ('out/power-two-divide-' + name) for name in names]
roots += [root / ('bin/linux-gcc/power-two-divide-' + name) for name in names]
roots += [work]
query = lambda args, **kw: subprocess.check_output(args, **kw)


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def load(path):
    return json.loads(path.read_text())


assert sys.argv[1:] == ['--preserve'] and not archive.exists()
assert query(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not query(['git', 'status', '--porcelain'])
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + tag]).returncode == 1
controls = roots[2]
assert load(controls / 'disposition.json')['decision'] == 'PROCEED_TO_ORIGINAL_BENEFIT_GATES'
assert load(controls / 'control-acceptance.json')['passed']
assert load(controls / 'endpoint-route-acceptance.json')['passed']
assert load(roots[1] / 'arithmetic-result.json')['passed']
receipts = sorted(p for owned in roots for p in owned.rglob('owner.json'))
assert len(receipts) == 32
for path in receipts:
    owner = load(path)
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0, path
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors'], path
    assert all(item['absent'] for item in owner['observed_births']), path
    assert owner['source'] == {'revision': revision, 'status': ''}, path
    assert owner['whole_inputs_unchanged'] and owner['source_unchanged'], path
for name, expected in load(controls / 'protected-before.json').items():
    assert identity(root / name) == expected, name

index = work / 'source.index'
assert not index.exists()
environment = dict(os.environ, GIT_INDEX_FILE=str(index))
sources = []
for owned in [*roots[:3], work]:
    sources += [p for p in sorted(owned.iterdir()) if p.is_file()
                and p.suffix in ['.py', '.json', '.patch', '.cpp', '.slang']
                and not p.name.startswith('issue99-')]
sources += [roots[0] / 'candidate/src/sirius/kernels/portable_binary32.h']
sources += [controls / 'contract/software_contracts.py']
records = []
try:
    subprocess.run(['git', 'read-tree', revision], env=environment, check=True)
    for path in sources:
        relative = 'diagnostics/' + str(path.relative_to(root / 'out'))
        blob = query(['git', 'hash-object', '-w', str(path)], text=True).strip()
        subprocess.run(['git', 'update-index', '--add', '--cacheinfo',
                        '100644,' + blob + ',' + relative], env=environment, check=True)
        records.append({'path': relative, 'origin': str(path.relative_to(root)),
                        'blob': blob, **identity(path)})
    tree = query(['git', 'write-tree'], env=environment, text=True).strip()
    commit = query(['git', 'commit-tree', tree, '-p', revision], text=True,
                   input='Preserve issue #99 exact-power division compiler arithmetic and original controls\n').strip()
    subprocess.run(['git', 'update-ref', 'refs/tags/' + tag, commit], check=True)
finally:
    index.unlink(missing_ok=True)
subprocess.run(['git', 'push', 'origin', 'refs/tags/' + tag], check=True)
remote = query(['git', 'ls-remote', 'origin', 'refs/tags/' + tag], text=True).strip()
assert remote == commit + '\trefs/tags/' + tag
assert query(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not query(['git', 'status', '--porcelain'])
(work / 'preserved-source.json').write_text(json.dumps({
    'tag': tag, 'commit': commit, 'parent': revision, 'source_records': records,
    'remote_readback': remote, 'working_tree_and_head_unchanged': True}, indent=2) + '\n')

archive.mkdir(parents=True)
files, modes, directories, directory_identities = {}, {}, {}, {}
for owned in roots:
    assert owned.is_dir() and not owned.is_symlink()
    for path in [owned, *sorted(owned.rglob('*'))]:
        mode = path.lstat().st_mode
        assert not stat.S_ISLNK(mode) and (stat.S_ISREG(mode) or stat.S_ISDIR(mode)), path
        name = str(path.relative_to(root))
        target = archive / 'files' / name
        if stat.S_ISDIR(mode):
            target.mkdir(parents=True, exist_ok=True)
            target.chmod(mode & 0o777)
            directories[name] = oct(mode & 0o777)
            saved = path.stat()
            directory_identities[name] = {'device': saved.st_dev, 'inode': saved.st_ino}
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        expected = identity(path)
        shutil.copy2(path, target)
        assert identity(path) == expected and identity(target) == expected, name
        files[name] = expected
        modes[name] = oct(mode & 0o777)
        assert oct(target.stat().st_mode & 0o777) == modes[name]
assert {str(p.relative_to(root)) for owned in roots for p in owned.rglob('*') if p.is_file()} == set(files)
manifest = {
    'source_revision': revision, 'evidence_tag': tag, 'evidence_commit': commit,
    'roots': [str(p.relative_to(root)) for p in roots], 'files': files,
    'preserved_file_modes': modes, 'preserved_directories': directories,
    'live_directory_identities': directory_identities, 'file_count': len(files),
    'total_bytes': sum(pin['bytes'] for pin in files.values()),
    'terminal_owner_receipts': [str(p.relative_to(root)) for p in receipts],
    'scope': 'Exact local raw compiler/arithmetic/eight-control prerequisites; source and concise results pushed separately. No benefit, production adoption or full qualification claim.'}
(archive / 'whole-byte-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({'archive': str(archive.relative_to(root)), 'files': len(files),
                  'bytes': manifest['total_bytes'], 'terminal_owners': len(receipts),
                  'tag': tag, 'commit': commit, 'pushed': True}))
