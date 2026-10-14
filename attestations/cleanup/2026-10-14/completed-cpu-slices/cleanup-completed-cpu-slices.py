"""Remove only preserved terminal scratch; retain a precise local census receipt."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile

root = Path.cwd()
source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert source == 'f5d649bf5067d1d68dc9215462f37830122f0540'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
targets = [root / p for p in ('out/science-p4-cpu-remaining-separated', 'out/science-p4-doppler-raster', 'out/science-e1-field-equations', 'out/science-p1-e4-oracle-current-source')]
tags = {'evidence/issue-53-p4-cpu-separated-current-source': '875f963cf532aa96cf3bb34ef832dc6fba57c0b6',
        'evidence/issue-53-p4-doppler-raster-current-source': '7ac9751de7603183bc66173672c5f219ce4d2e87',
        'evidence/issue-53-e1-field-equations-current-source': '28f91b8530a8d9de79402746b22b7f5f46bb3f1c',
        'evidence/issue-53-p1-e4-oracle-current-source': '0e7ac4b9667a61619e8d10cc4e3d5ce8c89c74fd'}
for tag, commit in tags.items():
    assert subprocess.check_output(['git', 'rev-parse', 'refs/tags/' + tag], text=True).strip() == commit
for packet in ('out/forensic-alcubierre-completion', 'out/science-p4-vulkan-session-scenes', 'out/science-p4-ordinary-fp32-parity'):
    work = root / packet
    actions = json.loads((work / 'actions.json').read_text())
    inputs = [Path(name) for action in actions.values() for name in action['frozen_inputs']]
    prep = json.loads((work / 'preparation.json').read_text())
    inputs += [Path(name) for name in prep.get('pins', {})]
    if (work / 'pins.json').exists():
        inputs += [Path(name) for name in json.loads((work / 'pins.json').read_text())['files']]
    for item in inputs:
        full = item if item.is_absolute() else root / item
        assert all(full != target and not full.is_relative_to(target) for target in targets), (packet, str(full))
for work in targets:
    for path in work.glob('*/owner.json'):
        owner = json.loads(path.read_text())
        assert owner['status'] == 'terminal' and owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert all(row['absent'] for row in owner['observed_births'])

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
        for descriptor in (entry / 'fd').iterdir():
            try:
                data += os.readlink(descriptor).encode()
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
    assert path.is_dir() and not path.is_symlink() and path.resolve().is_relative_to(root)
    files = [p for p in path.rglob('*') if p.is_file() and not p.is_symlink()]
    removed.append({'path': str(path.relative_to(root)), 'regular_files': len(files),
                    'regular_bytes': sum(p.stat().st_size for p in files)})
    shutil.rmtree(path)
    assert not path.exists()
bundle = root / 'attestations/cleanup/2026-10-14/completed-cpu-slices'
assert not bundle.exists()
bundle.mkdir(parents=True)
receipt = {'schema': 'sirius-preserved-terminal-scratch-cleanup-v1', 'runner_revision': source,
           'preserved_evidence_tags': tags, 'removed': removed, 'accessible_local_consumers': consumers,
           'local_census_permission_gaps': gaps, 'global_consumer_absence_claimed': False,
           'unfinished_alcubierre_inputs_untouched': True, 'qualification_claimed': False}
(bundle / 'receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
shutil.copy2(Path(__file__), bundle / Path(__file__).name)
fd, index = tempfile.mkstemp(dir=root / 'out', prefix='cleanup-receipt-')
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
               input=b'Record cleanup of preserved terminal CPU science slices\n').decode().strip()
    tag = 'evidence/issue-53-completed-cpu-slices-cleanup'
    subprocess.run(['git', 'tag', tag, commit], check=True)
    for path in bundle.iterdir():
        assert subprocess.check_output(['git', 'show', commit + ':' + str(path.relative_to(root))]) == path.read_bytes()
    print(json.dumps({'tag': tag, 'commit': commit, 'removed': removed, 'permission_gaps': len(gaps),
                      'receipt_sha256': hashlib.sha256((bundle / 'receipt.json').read_bytes()).hexdigest()}))
finally:
    if Path(index).exists():
        Path(index).unlink()
