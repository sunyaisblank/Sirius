import hashlib
import json
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
COMMIT = '1b02b28a3a29e30f0463fa9ac1578b92f9d6e30e'
PREFIX = 'attestations/qualification/2c8ef086/native-coefficient-sharing-source-compiler-driver'


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'mode': oct(stat.S_IMODE(path.stat().st_mode))}


receipt = json.loads((WORK / 'milestone-Git.json').read_text())
assert receipt['commit'] == COMMIT and receipt['pushed'] and receipt['fresh_remote_verified']
manifest = json.loads(git('show', COMMIT + ':' + PREFIX + '/MANIFEST.json'))
roots = ['out/native-rounded-error-closeout', 'out/native-rounded-error-preservation']
binding = json.loads((ROOT / 'out/native-coefficient-sharing-runtime-timestamp/execution-bindings.json').read_text())
assert not any(str((ROOT / name).resolve()) in str((ROOT / value).resolve())
               for name in roots for value in binding['files'])
inventories = {}
for name in roots:
    path = ROOT / name
    assert path.is_dir() and not path.is_symlink()
    expected = {row['source']: row for row in manifest['payloads'].values() if row['source'].startswith(name + '/')}
    actual = {}
    for item in path.rglob('*'):
        assert not item.is_symlink()
        if item.is_file():
            actual[str(item.relative_to(ROOT))] = pin(item)
    assert set(actual) == set(expected)
    for member, row in actual.items():
        assert row == {key: expected[member][key] for key in ('bytes', 'sha256', 'mode')}
        stored = git('show', COMMIT + ':' + PREFIX + '/files/' + member)
        assert len(stored) == row['bytes'] and hashlib.sha256(stored).hexdigest() == row['sha256']
    inventories[name] = actual
for name in roots:
    shutil.rmtree(ROOT / name)
report = {'schema': 'sirius-rounded-error-final-owned-scratch-retirement-v1', 'passed': True,
          'preserved_commit': COMMIT, 'retired_roots': inventories,
          'retired_files': sum(len(rows) for rows in inventories.values()),
          'retired_bytes': sum(row['bytes'] for rows in inventories.values() for row in rows.values()),
          'current_candidate_bindings_do_not_consume_these_roots': True,
          'scope': 'Completed rejected candidate preservation/cleanup helpers only; exact immutable bytes/modes verified, new live native control untouched.'}
(WORK / 'prior-scratch-retirement.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: value for key, value in report.items() if key != 'retired_roots'}))
