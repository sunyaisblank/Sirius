"""Retire only complete, pushed, byte-verified native build downloads."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

sys.dont_write_bytecode = True
ROOT = Path.cwd()
TEMPLATE = ROOT / 'attestations/cleanup/2026-10-15/terminal-isolated-parity-preflight-hints/cleanup-terminal-parity-preflight-hints.py'
assert hashlib.sha256(TEMPLATE.read_bytes()).hexdigest() == 'b5124e2a63b6621206dd3889040f6e9a67a42740655b050ffe18f958dc548aa2'
spec = importlib.util.spec_from_file_location('preserved_cleanup_functions', TEMPLATE)
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
WORK = ROOT / 'out/native-provider-integration/build-domains'
common.ROOT = ROOT
common.TARGETS = (WORK,)
common.prefix_target = lambda path, targets=None: any(
    path == target or path.is_relative_to(target)
    for target in (common.TARGETS if targets is None else targets))
common.BUNDLE = ROOT / 'attestations/cleanup/2026-10-15/completed-native-build-downloads'
common.TAG = 'evidence/issue-48-completed-native-build-downloads-cleanup'
common.clean_source()
assert not common.BUNDLE.exists()
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', 'refs/tags/' + common.TAG], cwd=ROOT).returncode == 1
bundle = ROOT / 'attestations/native-vulkan/build-domains-5b5fad9'
preserved = common.preserved((
    'evidence/issue-48-native-build-domains-5b5fad9',
    '939916066066dfc16ff21171f791c48ae1792f72', str(bundle.relative_to(ROOT)),
    '4e7a0c6d2cab6981679300951a0e7d90d5ae10349aca62950dd9f52e18bc47d0', ()))
local = common.files_under(WORK)
expected = {}
for path in (bundle / 'verification').rglob('*'):
    if path.is_file():
        relative = str(path.relative_to(bundle / 'verification'))
        assert (WORK / relative).read_bytes() == path.read_bytes()
        expected[relative] = common.identity(path)
verification = json.loads((WORK / 'verification.json').read_text())
assert verification['passed'] and not verification['new_runtime_domain_admitted']
parts = json.loads((bundle / 'result.json').read_text())['archive_parts']
for role in ('macos', 'windows'):
    archive = WORK / (role + '.zip')
    raw = archive.read_bytes()
    assert raw == b''.join((bundle / row['path']).read_bytes() for row in parts[role])
    assert common.identity(archive) == verification['results'][role]['archive']
    expected[role + '.zip'] = common.identity(archive)
    with zipfile.ZipFile(archive) as zipped:
        members = verification['results'][role]['members']
        assert len(zipped.infolist()) == len(members) == 40
        assert set(zipped.namelist()) == set(members)
        for name, pin in members.items():
            relative = role + '/' + name
            assert relative in local
            assert common.identity(local[relative]) == pin
            assert zipped.read(name) == local[relative].read_bytes()
            expected[relative] = pin
assert set(local) == set(expected)
for name, pin in expected.items():
    assert common.identity(local[name]) == pin
initial = common.consumer_census()
common.clean_source()
final = common.consumer_census()
assert common.files_under(WORK).keys() == local.keys()
receipt = {'schema': 'sirius-preserved-native-build-download-cleanup-v1',
           'source': common.SOURCE, 'preservation': preserved,
           'review_inputs_released': True, 'initial_consumer_census': initial,
           'final_consumer_census': final, 'global_consumer_absence_claimed': False,
           'removed': {'path': str(WORK.relative_to(ROOT)), 'regular_files': len(local),
                       'regular_bytes': sum(pin['bytes'] for pin in expected.values())},
           'active_native_provider_source_and_ci_receipts_retained': True,
           'qualification_claimed': False, 'removal_completed': False}
common.BUNDLE.mkdir(parents=True)
shutil.copy2(Path(__file__), common.BUNDLE / Path(__file__).name)
common.write(common.BUNDLE / 'receipt.json', receipt)
assert WORK.resolve() == WORK and WORK.is_dir() and not WORK.is_symlink()
shutil.rmtree(WORK)
assert not WORK.exists()
common.clean_source()
receipt['removal_completed'] = True
common.write(common.BUNDLE / 'receipt.json', receipt)
common.write(common.BUNDLE / 'manifest.json', {
    'qualification_claimed': False,
    'files': {p.name: common.identity(p) for p in common.BUNDLE.iterdir() if p.is_file()}})
fd, index = tempfile.mkstemp(dir=ROOT / 'out', prefix='native-download-cleanup-index-')
os.close(fd)
os.unlink(index)
environment = dict(os.environ, GIT_INDEX_FILE=index)
try:
    subprocess.run(['git', 'read-tree', '--empty'], cwd=ROOT, env=environment, check=True, timeout=30)
    rows = []
    paths = sorted(common.BUNDLE.iterdir())
    for path in paths:
        blob = common.git('hash-object', '-w', '--', str(path)).decode().strip()
        rows.append(f'100644 {blob}\t{path.relative_to(ROOT)}\0')
    subprocess.run(['git', 'update-index', '-z', '--index-info'], cwd=ROOT, env=environment,
                   input=''.join(rows).encode(), check=True, timeout=30)
    tree = common.git('write-tree', environment=environment).decode().strip()
    commit = subprocess.check_output(['git', 'commit-tree', tree, '-p', common.SOURCE],
        cwd=ROOT, env=environment, input=b'Record cleanup of preserved native build downloads\n',
        timeout=30).decode().strip()
    subprocess.run(['git', 'tag', common.TAG, commit], cwd=ROOT, check=True, timeout=30)
    for path in paths:
        assert common.git('show', commit + ':' + str(path.relative_to(ROOT))) == path.read_bytes()
    print(json.dumps({'tag': common.TAG, 'commit': commit, 'all_git_bytes_verified': True,
                      'removed': receipt['removed']}))
finally:
    if Path(index).exists():
        Path(index).unlink()
