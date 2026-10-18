import datetime, hashlib, json, pathlib, stat, subprocess, sys, zipfile

WORK = pathlib.Path(__file__).resolve().parent
RUN = 39116392750
SOURCE = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
ALLOWED = {'sirius-windows-diagnostic-build-products', 'sirius-windows-native-attestation',
           'sirius-macos-native-attestation', 'sirius-gate-diagnostics-integration-windows-no-render',
           'sirius-gate-diagnostics-integration-macos-no-render'}
snapshot = pathlib.Path(sys.argv[1])
data = json.loads(snapshot.read_text())
assert {a['name'] for a in data['artifacts']} <= ALLOWED
for item in data['artifacts']:
    name, aid = item['name'], item['id']
    assert item['workflow_run']['head_sha'] == SOURCE and item['workflow_run']['id'] == RUN
    archive = WORK / 'archives' / f'{aid}-{name}.zip'
    receipt = WORK / f'artifact-{aid}.json'
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        assert hashlib.sha256(archive.read_bytes()).hexdigest() == saved['archive']['sha256']
        print(name + ': already verified', flush=True)
        continue
    archive.parent.mkdir(exist_ok=True)
    assert not archive.exists()
    with archive.open('xb') as stream:
        result = subprocess.run(['gh', 'api', f'repos/sunyaisblank/Sirius/actions/artifacts/{aid}/zip'],
                                stdout=stream, stderr=subprocess.PIPE, timeout=180, check=True)
    b = archive.read_bytes()
    pin = dict(path=str(archive), bytes=len(b), sha256=hashlib.sha256(b).hexdigest())
    assert len(b) == item['size_in_bytes']
    assert item['digest'] == 'sha256:' + pin['sha256']
    folder = WORK / 'files' / name
    assert not folder.exists()
    members, folded = {}, set()
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for info in z.infolist():
            p = pathlib.PurePosixPath(info.filename)
            assert not p.is_absolute() and '..' not in p.parts and '\\' not in info.filename
            assert ':' not in info.filename and info.filename.casefold() not in folded
            folded.add(info.filename.casefold())
            mode = info.external_attr >> 16
            assert not stat.S_ISLNK(mode)
            if info.is_dir():
                continue
            assert not stat.S_IFMT(mode) or stat.S_ISREG(mode)
            payload = z.read(info)
            dest = folder.joinpath(*p.parts)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(payload)
            restore = 0o755 if dest.suffix.lower() in ('.exe', '.dll') else (stat.S_IMODE(mode) or 0o644)
            dest.chmod(restore)
            assert dest.read_bytes() == payload
            members[info.filename] = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(),
                                          mode=oct(restore), zip_external_attr=info.external_attr)
    receipt.write_text(json.dumps(dict(schema='sirius-native-constant-canonical-CI-artifact-v1',
                                      created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                      source_revision=SOURCE, run_id=RUN, artifact=item, archive=pin,
                                      API_digest_equal=True, safe_CRC_regular_members=True,
                                      extracted_folder=str(folder), members=members), indent=2) + '\n')
    print(json.dumps(dict(name=name, members=len(members), archive_bytes=len(b), passed=True)), flush=True)
