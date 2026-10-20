import hashlib
import json
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
FOUNDATION = '1b02b28a3a29e30f0463fa9ac1578b92f9d6e30e'
PRIVATE = 'e1c618f83861c478fe18c5a22f5d9aa27f6280fd'
PREFIX = 'attestations/qualification/2c8ef086/native-coefficient-sharing-source-compiler-driver'


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


assert not git('status', '--porcelain')
refs = git('ls-remote', 'origin',
           'refs/heads/development/native-coefficient-sharing-diagnostic-2026-10-20',
           'refs/tags/evidence/native-coefficient-sharing-compiler-driver-2c8ef086').decode()
remote = {row.split()[1]: row.split()[0] for row in refs.splitlines()}
assert remote == {'refs/heads/development/native-coefficient-sharing-diagnostic-2026-10-20': PRIVATE,
                  'refs/tags/evidence/native-coefficient-sharing-compiler-driver-2c8ef086': FOUNDATION}
manifest = json.loads(git('show', FOUNDATION + ':' + PREFIX + '/MANIFEST.json'))
roots = ['out/native-coefficient-sharing-' + name for name in
         ['source', 'compiler', 'emitted', 'producer-source', 'producer-inputs', 'windows-source', 'driver']]
# All current compiler/query/runtime owners are terminal and already preserved.
# No next candidate or live consumer uses these completed private inputs.
tracked = set(git('ls-files', '-z').decode().split('\0'))
inventory = {}
for name in roots:
    source = ROOT / name
    assert source.is_dir() and not source.is_symlink() and source.resolve().is_relative_to(ROOT)
    rows = {}
    for path in source.rglob('*'):
        assert not path.is_symlink()
        if not path.is_file():
            continue
        rel = str(path.relative_to(ROOT))
        assert rel not in tracked
        key = 'files/' + rel
        if key in manifest['payloads']:
            record = manifest['payloads'][key]
            commit, stored = FOUNDATION, PREFIX + '/' + key
        else:
            commit, stored = PRIVATE, rel
            # Git records only the executable mode; this source set is 0644.
            record = {**pin(git('show', commit + ':' + stored)), 'mode': '0o644'}
        actual = {**pin(path.read_bytes()), 'mode': oct(stat.S_IMODE(path.stat().st_mode))}
        assert actual == {key: record[key] for key in ('bytes', 'sha256', 'mode')}, rel
        assert pin(git('show', commit + ':' + stored)) == {key: record[key] for key in ('bytes', 'sha256')}, rel
        assert git('ls-tree', commit, '--', stored).decode().split()[0] == '100644'
        rows[rel] = {**actual, 'Git_commit': commit, 'Git_path': stored}
    inventory[name] = rows
# No old runtime/compiler/provider child remains; the new compiler and driver
# query are terminal. No future new-candidate execution depends on these roots.
for name in roots:
    shutil.rmtree(ROOT / name)
    assert not (ROOT / name).exists()
report = {'schema': 'sirius-native-coefficient-sharing-published-input-retirement-v1', 'passed': True,
          'remote_refs': remote, 'retired_roots': inventory,
          'retired_files': sum(len(rows) for rows in inventory.values()),
          'retired_bytes': sum(row['bytes'] for rows in inventory.values() for row in rows.values()),
          'scope': 'Only exact Git-backed completed rejected candidate source/compiler/producer input roots. Canonical source/builds and host providers remain. No new candidate has been started.'}
(WORK / 'published-input-retirement.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: value for key, value in report.items() if key != 'retired_roots'}))
