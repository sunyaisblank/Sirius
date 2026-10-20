import hashlib
import io
import json
import stat
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
COMMIT = 'b9a5703a18f06070fe809a9bfd14a414c4f2f906'
PREFIX = 'attestations/qualification/2c8ef086/native-coefficient-sharing-controls-benefit-outcome'


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def signature(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def match(data, row):
    assert signature(data) == {k: row[k] for k in ('bytes', 'sha256')}


manifest_data = git('show', COMMIT + ':' + PREFIX + '/MANIFEST.json')
assert hashlib.sha256(manifest_data).hexdigest() == '0bdb55c81c690bcdcffeac1f3c05f781c8a0c0f921be7aad0419ffb222f283c4'
manifest = json.loads(manifest_data)
payloads = manifest['payloads']
assert len(payloads) == 242
receipt = json.loads((ROOT / 'out/native-coefficient-sharing-preservation/root-Git.json').read_text())
assert receipt['commit'] == COMMIT and receipt['pushed'] and receipt['fresh_remote_verified']
assert manifest['outcome']['checks_passed'] == 12
assert manifest['outcome']['checks_total'] == 22
assert manifest['outcome']['original_controls_passed'] == 1
assert manifest['outcome']['remaining_seven_state'] == 'not-run'
assert not manifest['outcome']['source_adopted']
assert not manifest['outcome']['full_scientific_or_release_qualification_claimed']
sources = {}
for name, row in payloads.items():
    source, copy = ROOT / row['source'], ROOT / PREFIX / name
    assert source.is_file() and copy.is_file()
    assert not source.is_symlink() and not copy.is_symlink()
    for path in (source, copy):
        match(path.read_bytes(), row)
        assert oct(stat.S_IMODE(path.stat().st_mode)) == row['mode']
    sources[row['source']] = row
for name, row in receipt['files'].items():
    match(git('show', COMMIT + ':' + name), row)
    assert git('ls-tree', COMMIT, '--', name).decode().split()[0] == row['Git_mode']
assert len(receipt['files']) == 246
stage = manifest['archive_backed_stage']
assert len(stage) == 122
archive_names = {row['archive_payload'] for row in stage.values()}
assert len(archive_names) == 1
archive_name = archive_names.pop()
archive_data = git('show', COMMIT + ':' + PREFIX + '/' + archive_name)
match(archive_data, payloads[archive_name])
with zipfile.ZipFile(io.BytesIO(archive_data)) as archive:
    infos = archive.infolist()
    assert len(infos) == 122 and len({info.filename for info in infos}) == 122
    assert {info.filename for info in infos} == {row['member'] for row in stage.values()}
    for name, row in stage.items():
        info = archive.getinfo(row['member'])
        assert info.external_attr == row['zip_external_attr']
        match(archive.read(info), row)
        path = ROOT / name
        assert path.is_file() and not path.is_symlink()
        match(path.read_bytes(), row)
        assert oct(stat.S_IMODE(path.stat().st_mode)) == row['mode'] == row['restore_mode']
for reference in manifest['external_Git_sources_not_copied']:
    assert git('cat-file', '-t', reference['commit']).strip() == b'commit'
    for prefix in ([reference['prefix']] if 'prefix' in reference else reference.get('prefixes', [])):
        assert git('ls-tree', reference['commit'], '--', prefix)
preservation = json.loads(git('show', COMMIT + ':' + PREFIX + '/PRESERVATION.json'))
assert preservation['outcome'] == manifest['outcome']
assert preservation['payload_count'] == len(payloads)
assert preservation['payload_bytes'] == sum(row['bytes'] for row in payloads.values())
assert preservation['archive_backed_stage_files'] == len(stage)
result = {'schema': 'sirius-native-coefficient-sharing-separate-root-conservation-v1',
          'passed': True, 'commit': COMMIT, 'manifest': signature(manifest_data),
          'payload_files': len(payloads), 'payload_bytes': preservation['payload_bytes'],
          'Git_files_bytes_and_exec_modes_verified': len(receipt['files']),
          'ZIP_members_bytes_attrs_and_live_restore_modes_verified': len(stage),
          'outcome': manifest['outcome'],
          'limits': 'Separate stdlib conservation implementation by ROOT, not an independent agent or new scientific execution. Host dependency hash seals are not backups. No deletion performed.'}
(WORK / 'conservation.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
