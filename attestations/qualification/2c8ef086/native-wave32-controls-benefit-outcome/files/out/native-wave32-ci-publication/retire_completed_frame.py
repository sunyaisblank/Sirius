import hashlib
import json
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
COMMIT = 'b8ecca0ff4d7fe399db5f37922416f589e822e98'
MANIFEST = 'attestations/qualification/2c8ef086/native-constant-platform-builds-and-frame/MANIFEST.json'


def run(*argv):
    return subprocess.check_output(argv, cwd=ROOT)


manifest = json.loads(run('git', 'show', COMMIT + ':' + MANIFEST))
receipt = json.loads((ROOT / 'out/native-constant-preservation/root-platform-frame-Git.json').read_text())
assert receipt['commit'] == COMMIT and receipt['pushed'] and receipt['fresh_remote_verified']
insight = json.loads((WORK / 'driver-insight-Git.json').read_text())
assert insight['pushed'] and insight['fresh_remote_verified']
expected = {name: dict(row['files']) for name, row in manifest['owned_source_roots'].items()}
for name, row in manifest['archive_backed_expanded_CI'].items():
    for member, item in row['members'].items():
        expected['out/native-constant-current-ci'][name + '/' + member] = item
stage = manifest['archive_backed_stage']
expected[stage['root']] = {stage['root'] + '/' + name: row for name, row in stage['members'].items()}
assert set(expected) == {'out/native-constant-current-ci', 'out/native-constant-original-frame',
                         'renders/constant-original-frame-2c8ef086',
                         'bin/windows-msvc/constant-original-frame-2c8ef086'}
tracked = run('git', 'ls-files', '-z').decode().split('\0')
inventories = {}
for name, rows in expected.items():
    path = ROOT / name
    assert path.is_dir() and not path.is_symlink() and path.resolve().is_relative_to(ROOT)
    assert not any(p.startswith(name + '/') for p in tracked)
    found = {}
    for item in path.rglob('*'):
        assert not item.is_symlink() and (item.is_dir() or item.is_file())
        if item.is_file():
            data = item.read_bytes()
            found[str(item.relative_to(ROOT))] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                                                  'mode': oct(stat.S_IMODE(item.stat().st_mode))}
    assert set(found) == set(rows), name
    for member, actual in found.items():
        assert all(actual[key] == rows[member][key] for key in ['bytes', 'sha256', 'mode']), member
    inventories[name] = found
assert not run('git', 'status', '--porcelain')
# Completed scientific child and both diagnostic query children have exited;
# saved owner/bridge records and independent shutdown/result reviews are conserved.
for name in expected:
    shutil.rmtree(ROOT / name)
    assert not (ROOT / name).exists()
record = {'schema': 'sirius-completed-constant-frame-retirement-v1',
          'preserved_commit': COMMIT, 'preserved_manifest': MANIFEST,
          'driver_insight_commit': insight['commit'], 'all_live_files_and_modes_matched_preserved_inventory': True,
          'retired_roots': inventories, 'retired_files': sum(len(rows) for rows in inventories.values()),
          'retired_bytes': sum(row['bytes'] for rows in inventories.values() for row in rows.values()),
          'scope': 'Only four completed owned CI/frame/render/stage roots; canonical builds, new candidate/producer/runtime sources and shared host providers remain.'}
(WORK / 'completed-frame-retirement.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({key: value for key, value in record.items() if key != 'retired_roots'}))
