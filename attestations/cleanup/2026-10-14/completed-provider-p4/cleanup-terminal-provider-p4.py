from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import zipfile

root = Path.cwd()
source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert source == 'f5d649bf5067d1d68dc9215462f37830122f0540'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
tags = {'evidence/issue-53-p4-disk-volume-current-source': '3cd4d6cb65fe188a8cbae9ec75c596e551f63187',
        'evidence/issue-53-p4-cpu-interrupted-batch': '9722b47f957dee21afcbef87c6743122368495a2',
        'evidence/issue-48-provider-identity-fixed': '21f87e4806c36cd1ae9e45f912a5817d6878b774',
        'evidence/issue-48-provider-discovery-f5d649b': '626a5d80627c41704dbf1c2e0a8d2f6bf149547c',
        'evidence/issue-48-provider-ninja-escape-a7a0c9d': '30962d93c257c449dd1830d1e20506541ed1c786',
        'evidence/issue-48-provider-build-side-effects-bfefab7': 'e37b59b0d9c05b9d274e79f4fc52435aea62c554',
        'evidence/issue-48-provider-metadata-compile-a75d4b5': '4f8b9a00fa8db4e7e72a7ff0ad7cbba927e1894f'}
for tag, commit in tags.items():
    assert subprocess.check_output(['git', 'rev-parse', 'refs/tags/' + tag], text=True).strip() == commit
targets = [root / p for p in ('out/science-p4-disk-volume', 'out/science-p4-cpu-remainder',
                             'out/macos-provider-metadata-control', 'bin/linux-gcc/macos-provider-metadata-control',
                             'out/prepare-macos-metadata-control.py')]
for relative in ('out/science-p4-disk-volume/original-four-disk-volume-controls',
                 'out/science-p4-cpu-remainder/original-four-non-frame-cpu-controls',
                 'out/macos-provider-metadata-control/original-compile',
                 'out/macos-provider-metadata-control/corrected-compile',
                 'out/macos-provider-metadata-control/runtime-owner/corrected-run'):
    owner = json.loads((root / relative / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert all(row['absent'] for row in owner['observed_births'])
archives = {38556777741: '2b0dbcd768b646a9d258c1942c797e3a642a96ece71a4f2c2d51e2da00acf030',
            38557811215: '7cd11604400b28bbe3c39e55cf1f4d86ec8f1cd64dc2aa67bcb79c0f6ca74c6d',
            38559236199: '56d46acefdbffb7e84fdcf6835eb48608d286e04d1038f99c10a47d92d5fcb6d',
            38560779867: 'a86193c41a5f8fd5d5f431f31e8d4a9a3427169e25b92f2456386ef4e4e1f9cb'}
for run, digest in archives.items():
    path = root / f'out/macos-moltenvk-integration/native-{run}.zip'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    with zipfile.ZipFile(path) as z:
        prefix = 'out/macos-moltenvk-integration/observation/'
        phase = 'toolchain' if run == 38556777741 else 'provider-build'
        owner = json.loads(z.read(prefix + phase + '/owner.json'))
        assert owner['status'] == 'terminal' and not owner['cleanup_errors'] and not owner['remaining_observed']
        assert all(row['absent'] for row in owner['observed_births'])
    targets += [path, root / f'out/macos-moltenvk-integration/run-{run}-failure.log']

# Check actual future input lists; copied historical path strings are provenance.
for packet in ('out/forensic-alcubierre-completion', 'out/science-p4-cpu-remaining-separated', 'out/science-p4-doppler-raster'):
    work = root / packet
    actions = json.loads((work / 'actions.json').read_text())
    inputs = [Path(name) for action in actions.values() for name in action['frozen_inputs']]
    pins_path = work / 'pins.json'
    if pins_path.exists():
        inputs += [Path(name) for name in json.loads(pins_path.read_text())['files']]
    preparation = json.loads((work / 'preparation.json').read_text())
    inputs += [Path(name) for name in preparation.get('pins', {})]
    for item in inputs:
        full = item if item.is_absolute() else root / item
        assert all(full != target and not full.is_relative_to(target) for target in targets), (packet, str(full))
consumers, gaps = [], []
needles = [str(path).encode() for path in targets]
for entry in Path('/proc').glob('[0-9]*'):
    if int(entry.name) == os.getpid():
        continue
    try:
        data = (entry / 'cmdline').read_bytes() + (entry / 'maps').read_bytes()
        for name in ('cwd', 'exe'):
            try:
                data += os.readlink(entry / name).encode()
            except FileNotFoundError:
                pass
        for fd in (entry / 'fd').iterdir():
            try:
                data += os.readlink(fd).encode()
            except FileNotFoundError:
                pass
        if any(needle in data for needle in needles):
            consumers.append({'pid': int(entry.name)})
    except (FileNotFoundError, ProcessLookupError):
        pass
    except PermissionError:
        gaps.append({'pid': int(entry.name), 'error': 'PermissionError'})
assert not consumers, consumers
removed = []
for path in targets:
    assert path.exists() and not path.is_symlink() and path.resolve().is_relative_to(root)
    files = [p for p in path.rglob('*') if p.is_file() and not p.is_symlink()] if path.is_dir() else [path]
    removed.append({'path': str(path.relative_to(root)), 'regular_files': len(files),
                    'regular_bytes': sum(p.stat().st_size for p in files)})
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()
    assert not path.exists()
bundle = root / 'attestations/cleanup/2026-10-14/completed-provider-p4'
assert not bundle.exists()
bundle.mkdir(parents=True)
receipt = {'schema': 'sirius-preserved-terminal-scratch-cleanup-v1', 'runner_revision': source,
           'preserved_evidence_tags': tags, 'removed': removed, 'accessible_local_consumers': consumers,
           'local_census_permission_gaps': gaps, 'global_consumer_absence_claimed': False,
           'unfinished_alcubierre_and_p4_inputs_untouched': True, 'qualification_claimed': False}
(bundle / 'receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
shutil.copy2(Path(__file__), bundle / Path(__file__).name)
fd, index = tempfile.mkstemp(dir=root / 'out', prefix='cleanup-provider-p4-')
os.close(fd)
os.unlink(index)
env = dict(os.environ, GIT_INDEX_FILE=index)
try:
    subprocess.run(['git', 'read-tree', '--empty'], env=env, check=True)
    rows = []
    for path in sorted(bundle.iterdir()):
        blob = subprocess.check_output(['git', 'hash-object', '-w', '--', str(path)], text=True).strip()
        rows.append(f'100644 {blob}\t{path.relative_to(root)}\0')
    subprocess.run(['git', 'update-index', '-z', '--index-info'], input=''.join(rows).encode(), env=env, check=True)
    tree = subprocess.check_output(['git', 'write-tree'], env=env, text=True).strip()
    commit = subprocess.check_output(['git', 'commit-tree', tree, '-p', source], env=env,
               input=b'Remove preserved terminal provider and P4 scratch\n').decode().strip()
    tag = 'evidence/issue-48-53-terminal-scratch-cleanup'
    subprocess.run(['git', 'tag', tag, commit], check=True)
    for path in bundle.iterdir():
        assert subprocess.check_output(['git', 'show', commit + ':' + str(path.relative_to(root))]) == path.read_bytes()
    print(json.dumps({'tag': tag, 'commit': commit, 'removed': removed, 'permission_gaps': len(gaps)}))
finally:
    if Path(index).exists():
        Path(index).unlink()
