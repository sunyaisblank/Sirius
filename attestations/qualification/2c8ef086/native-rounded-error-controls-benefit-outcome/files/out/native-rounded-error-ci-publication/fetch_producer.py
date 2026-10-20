import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
assert len(sys.argv) == 5, 'ARTIFACTS_JSON ACTUAL_RUN_ID ACTUAL_PRIVATE_SHA NEW_STAGE'
artifacts_path = ROOT / sys.argv[1]
RUN = int(sys.argv[2])
HEAD = sys.argv[3]
assert RUN > 0 and re.fullmatch('[0-9a-f]{40}', HEAD)
stage = ROOT / sys.argv[4]
assert stage.is_relative_to(ROOT / 'bin/windows-msvc') and stage.resolve() == stage
artifacts = json.loads(artifacts_path.read_text())['artifacts']
assert len(artifacts) == 1
artifact = artifacts[0]
assert artifact['name'] == 'sirius-native-rounded-error-matched-producer'
assert artifact['workflow_run']['id'] == RUN and artifact['workflow_run']['head_sha'] == HEAD
assert artifact['expired'] is False and 0 < artifact['size_in_bytes'] <= 300 * 1024 ** 2
expected_hash = artifact['digest']
assert expected_hash.startswith('sha256:') and len(expected_hash) == 71
archive = WORK / ('producer-%d.zip' % artifact['id'])
assert not archive.exists() and not stage.exists()
argv = ['gh', 'api', 'repos/sunyaisblank/Sirius/actions/artifacts/%d/zip' % artifact['id']]
with archive.open('xb') as stream:
    process = subprocess.Popen(argv, cwd=ROOT, stdout=stream, stderr=subprocess.PIPE)
    try:
        _, errors = process.communicate(timeout=180)
    except subprocess.TimeoutExpired:
        process.kill()
        _, errors = process.communicate()
        (WORK / 'download.stderr').write_bytes(errors)
        raise
(WORK / 'download.stderr').write_bytes(errors)
assert process.returncode == 0
data_hash = hashlib.sha256()
with archive.open('rb') as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        data_hash.update(block)
assert archive.stat().st_size == artifact['size_in_bytes']
assert 'sha256:' + data_hash.hexdigest() == expected_hash
members = {}
with zipfile.ZipFile(archive) as bundle:
    rows = bundle.infolist()
    assert 0 < len(rows) <= 4096 and sum(row.file_size for row in rows) <= 1024 ** 3
    assert bundle.testzip() is None
    for row in rows:
        path = PurePosixPath(row.filename)
        assert not path.is_absolute() and '..' not in path.parts and '\\' not in row.filename
        assert path.parts and ':' not in row.filename
        assert not stat.S_ISLNK(row.external_attr >> 16)
        if row.is_dir():
            continue
        assert row.filename not in members and row.file_size <= 256 * 1024 ** 2
        data = bundle.read(row)
        target = stage / row.filename
        assert target.resolve().is_relative_to(stage)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(data)
        mode = 0o755 if target.suffix.lower() in ('.exe', '.dll') else 0o644
        os.chmod(target, mode)
        members[row.filename] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                                'mode': oct(mode), 'zip_external_attr': row.external_attr}
manifest = json.loads((stage / 'manifest.json').read_text())
producer = json.loads((stage / 'producer-result.json').read_text())
assert manifest['source_revision'] == producer['source_revision'] == HEAD
assert manifest['passed'] and producer['passed']
assert producer['producer_profile'] == 'native-rounded-error-five-shader-only'
assert type(producer['private_compile_object_count']) is int and type(producer['expected_private_compile_objects']) is int
assert producer['private_compile_object_count'] == producer['expected_private_compile_objects'] == 5
assert producer['canonical_CXX_source_revision'] == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert set(manifest['files']) | {'manifest.json'} == set(members)
for name, expected in manifest['files'].items():
    assert all(members[name][key] == expected[key] for key in ['bytes', 'sha256']), name
assert not (WORK / 'actual-download.json').exists()
record = {'schema': 'sirius-native-rounded-error-actual-download-v1', 'run_id': RUN, 'producer_source_revision': HEAD, 'artifact': artifact,
          'archive': str(archive.relative_to(ROOT)), 'archive_sha256': data_hash.hexdigest(),
          'stage': str(stage.relative_to(ROOT)), 'members': members,
          'API_digest_CRC_safe_members_and_manifest_exact': True,
          'scope': 'Transport/download and source/product manifest binding only; PE/MAP/runtime admission remains separate.'}
(WORK / 'actual-download.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({key: value for key, value in record.items() if key not in ['artifact', 'members']}))
