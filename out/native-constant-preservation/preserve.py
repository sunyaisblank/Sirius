import datetime, hashlib, json, os, pathlib, shutil, stat, subprocess, zipfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT)
BASE = '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'
PRODUCER = 'd2e1a265106a8b27e976ccb47cef6303aaa51746'
BUNDLE = pathlib.Path('attestations/qualification/5b2df58b/native-constant-controls-benefit-pass')
RECEIPTS = pathlib.Path('out/native-constant-preservation')

def run(*args, data=None, env=None):
    return subprocess.run(args, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          check=True, env=env).stdout

def identity(p):
    p = pathlib.Path(p)
    assert p.is_file() and not p.is_symlink(), p
    return dict(bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                mode=oct(stat.S_IMODE(p.stat().st_mode)))

def inventory(root):
    return {str(p): identity(p) for p in sorted(pathlib.Path(root).rglob('*')) if p.is_file()}

def write(p, value):
    pathlib.Path(p).write_text(json.dumps(value, indent=2) + '\n')

assert run('git', 'rev-parse', 'HEAD').decode().strip() == BASE
assert not run('git', 'status', '--porcelain')
review = json.loads(pathlib.Path('out/native-constant-publication/independent-all8-completion.json').read_text())
assert review['passed'] and review['original_controls_passed'] == 8 and review['original22_benefit_passed']
for pin in review['critical_pins'].values():
    current = identity(pin['path'])
    assert all(current[k] == pin[k] for k in ('bytes', 'sha256')), pin
roots = sorted(str(p) for p in pathlib.Path('out').glob('native-constant-*')
               if p.is_dir() and p != RECEIPTS)
roots += sorted(str(p) for p in pathlib.Path('renders').glob('native-constant-runtime-*') if p.is_dir())
assert len(roots) == 20, roots
assert not BUNDLE.exists()
owned = {r: inventory(r) for r in roots}
payloads = {}
for entries in owned.values():
    for source, pin in entries.items():
        dest = BUNDLE / 'files' / source
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        assert identity(dest) == pin
        payloads[str(dest.relative_to(BUNDLE))] = dict(pin, source=source)

download = json.loads(pathlib.Path('out/native-constant-ci-publication/root-actual-download.json').read_text())
stage = pathlib.Path(download['stage'])
stage_inventory = {str(pathlib.Path(p).relative_to(stage)): pin for p, pin in inventory(stage).items()}
assert len(stage_inventory) == download['member_count'] == 122
assert {k: {q: v[q] for q in ('bytes', 'sha256')} for k, v in stage_inventory.items()} == download['members']
assert all(v['mode'] == ('0o755' if pathlib.Path(k).suffix.lower() in ('.exe', '.dll') else '0o644')
           for k, v in stage_inventory.items())
archive = download['archive']['path']
assert all(identity(archive)[k] == download['archive'][k] for k in ('bytes', 'sha256'))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    assert sorted(z.namelist()) == sorted(stage_inventory)
    for k, v in stage_inventory.items():
        b = z.read(k)
        assert len(b) == v['bytes'] and hashlib.sha256(b).hexdigest() == v['sha256']

manifest = dict(schema='sirius-native-constant-controls-benefit-pass-manifest-v1',
                created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                canonical_source_revision=BASE, private_producer_revision=PRODUCER,
                outcome='PASS_ORIGINAL_8_CONTROLS_AND_SOLE_22_BENEFIT_GATES_SOURCE_INTEGRATION_ELIGIBLE',
                payloads=payloads, owned_source_roots=owned,
                archive_backed_stage=dict(root=str(stage), archive_payload='files/' + archive,
                                          member_count=len(stage_inventory), members=stage_inventory))
write(BUNDLE / 'MANIFEST.json', manifest)
preservation = dict(schema='sirius-native-constant-controls-benefit-pass-preservation-v1',
                    payload_files=len(payloads), payload_bytes=sum(p['bytes'] for p in payloads.values()),
                    source_root_count=len(owned), archive_backed_stage_members=len(stage_inventory),
                    expanded_stage_bytes=sum(p['bytes'] for p in stage_inventory.values()),
                    decision=manifest['outcome'], audit_passed=True, benefit_passed=True,
                    original_controls_passed=8, original_benefit_gates_passed=22,
                    source_producer_and_raw_records_complete_within_declared_bounded_scope=True,
                    canonical_source_changed=False, full_qualification_or_release=False,
                    scope='Exact c035 bounded native RTConstant candidate; actual Windows producer, SDK/AMD, 8 original controls and sole original ABBA. No frame/UI/full science/Mandatory/all-provider/installed release claim. Archive stored once; member modes permit exact POSIX stage restoration.')
write(BUNDLE / 'PRESERVATION.json', preservation)
assert {r: inventory(r) for r in roots} == owned
assert {str(pathlib.Path(p).relative_to(stage)): pin for p, pin in inventory(stage).items()} == stage_inventory
print(json.dumps(preservation))
