"""Preserve and verify one exact current official native build export."""
from pathlib import Path, PurePosixPath
import hashlib
import importlib.util
import json
import stat
import subprocess
import sys
import time
import zipfile

sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
RUN = 38077831907
HEAD = json.loads((WORK / 'source.json').read_text())['source_revision']

def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')

def api(path):
    return json.loads(subprocess.check_output(['gh', 'api', path], text=True, timeout=60))

assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'])
run = api(f'repos/sunyaisblank/Sirius/actions/runs/{RUN}')
write('export-ci-terminal.json', run)
assert run['head_sha'] == HEAD and run['status'] == 'completed' and run['conclusion'] == 'success'
jobs = api(f'repos/sunyaisblank/Sirius/actions/runs/{RUN}/jobs')
write('export-jobs-terminal.json', jobs)
expected = {'integration-no-render', 'integration-windows-no-render', 'integration-macos-no-render'}
assert {j['name'] for j in jobs['jobs'] if j['conclusion'] == 'success'} == expected
assert all(j['conclusion'] == ('success' if j['name'] in expected else 'skipped') for j in jobs['jobs'])
meta = api(f'repos/sunyaisblank/Sirius/actions/runs/{RUN}/artifacts')
write('export-artifacts.json', meta)
spec = importlib.util.spec_from_file_location('attestation', ROOT / 'scripts/verify-attestation.py')
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)
gate_verify = verify.load_build_gate_verifier()
for platform in ('windows', 'macos'):
    matches = [a for a in meta['artifacts'] if a['name'] == f'sirius-{platform}-native-attestation']
    assert len(matches) == 1
    a = matches[0]
    assert not a['expired'] and a['workflow_run']['id'] == RUN and a['workflow_run']['head_sha'] == HEAD
    archive = WORK / f'{platform}-export.zip'
    destination = ROOT / 'attestations/native-build' / HEAD[:7] / f'{platform}-build'
    assert not archive.exists() and not destination.exists()
    started = time.monotonic()
    argv = ['gh', 'api', f"repos/sunyaisblank/Sirius/actions/artifacts/{a['id']}/zip"]
    with archive.open('wb') as out, (WORK / f'{platform}-download.stderr').open('wb') as err:
        c = subprocess.run(argv, stdout=out, stderr=err, timeout=180)
    assert c.returncode == 0
    content = archive.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    assert a['digest'] == 'sha256:' + digest and len(content) == a['size_in_bytes']
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        names = set()
        for e in entries:
            p = PurePosixPath(e.filename)
            assert not p.is_absolute() and '..' not in p.parts and '\\' not in e.filename and ':' not in e.filename
            assert e.filename not in names and not stat.S_ISLNK(e.external_attr >> 16)
            names.add(e.filename)
        assert len([e for e in entries if not e.is_dir()]) == 40
        destination.mkdir(parents=True)
        z.extractall(destination)
    files = []
    for p in sorted(destination.rglob('*')):
        assert not p.is_symlink()
        if p.is_file():
            b = p.read_bytes()
            files.append({'path': str(p.relative_to(destination)), 'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()})
    assert len(files) == 40
    document = json.loads((destination / f'{platform}-build.json').read_text())
    verify.verify_path(destination / f'{platform}-build.json', ROOT)
    gate = gate_verify.validate_native_build_document(json.loads((destination / 'native_build_gate.json').read_text()))
    assert gate['source'] == {'revision': HEAD, 'clean': True}
    assert document['source_revision'] == HEAD
    result = {'source_revision': HEAD, 'platform': platform, 'artifact': a, 'argv': argv,
              'returncode': c.returncode, 'elapsed_seconds': time.monotonic() - started,
              'download_bytes': len(content), 'sha256': digest, 'extracted_files': files,
              'first_party_build_attestation_verified': True, 'runtime_qualification_claimed': False}
    write(f'{platform}-export-download.json', result)
    print(json.dumps({k: result[k] for k in ('platform', 'source_revision', 'download_bytes', 'sha256')}), flush=True)
