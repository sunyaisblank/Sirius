import hashlib
import json
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
COMMIT = 'b9a5703a18f06070fe809a9bfd14a414c4f2f906'
PREFIX = 'attestations/qualification/2c8ef086/native-coefficient-sharing-controls-benefit-outcome'


def run(*args):
    return subprocess.check_output(args, cwd=ROOT)


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'mode': oct(stat.S_IMODE(path.stat().st_mode))}


manifest = json.loads(run('git', 'show', COMMIT + ':' + PREFIX + '/MANIFEST.json'))
receipt = json.loads((ROOT / 'out/native-coefficient-sharing-preservation/root-Git.json').read_text())
assert receipt['commit'] == COMMIT and receipt['pushed'] and receipt['fresh_remote_verified']
review = json.loads((WORK / 'conservation.json').read_text())
assert review['passed'] is True
assert not run('git', 'status', '--porcelain')
tracked = run('git', 'ls-files', '-z').decode().split('\0')
expected = {}
absent = []
for name, row in manifest['owned_roots'].items():
    if name == 'out/native-coefficient-sharing-preservation':
        continue  # Receipt and preservation scripts stay until closeout is committed.
    if not row['exists']:
        assert not (ROOT / name).exists()
        absent.append(name)
        continue
    members = {entry['source']: entry for entry in manifest['payloads'].values()
               if entry['source'].startswith(name + '/')}
    assert len(members) == row['files'], name
    expected[name] = members
stage = 'bin/windows-msvc/native-coefficient-sharing-matched-current'
expected[stage] = manifest['archive_backed_stage']
inventories = {}
for name, rows in expected.items():
    path = ROOT / name
    assert path.is_dir() and not path.is_symlink() and path.resolve().is_relative_to(ROOT)
    assert not any(item.startswith(name + '/') for item in tracked)
    found = {}
    for item in path.rglob('*'):
        assert not item.is_symlink() and (item.is_file() or item.is_dir())
        if item.is_file():
            found[str(item.relative_to(ROOT))] = pin(item)
    assert set(found) == set(rows), name
    for member, actual in found.items():
        assert all(actual[key] == rows[member][key] for key in ['bytes', 'sha256', 'mode']), member
    inventories[name] = found
# All five original Windows-owned instances and terminal audits are closed and
# preserved. The new compiler seals do not depend on these retired roots.
# All five current scientific owners and their terminal audits are closed.
# No next candidate or live consumer uses any of these completed roots.
for name in expected:
    shutil.rmtree(ROOT / name)
    assert not (ROOT / name).exists()
record = {'schema': 'sirius-native-coefficient-sharing-outcome-retirement-v1', 'passed': True,
          'preserved_commit': COMMIT, 'preserved_manifest': PREFIX + '/MANIFEST.json',
          'separate_ROOT_conservation': pin(WORK / 'conservation.json'),
          'retired_roots': inventories, 'already_absent_roots': absent,
          'retired_files': sum(len(rows) for rows in inventories.values()),
          'retired_bytes': sum(row['bytes'] for rows in inventories.values() for row in rows.values()),
          'all_live_files_and_modes_matched_preserved_inventory': True,
          'limits': 'Only completed outcome runtime, reader, CI, audit, render and stage roots. '
                    'Private source/compiler inputs, preservation receipts, new coefficient-sharing '
                    'candidate, canonical builds and host providers remain outside this removal.'}
(WORK / 'retirement.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({key: value for key, value in record.items() if key != 'retired_roots'}))
