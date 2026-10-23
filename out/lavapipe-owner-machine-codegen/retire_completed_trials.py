"""Retire archive-verified completed trial raw files, retaining current consumers."""
from pathlib import Path
import hashlib, json, os, subprocess, zipfile

R = Path.cwd()
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
tracked = {R / os.fsdecode(p) for p in subprocess.check_output(['git', 'ls-files', '-z']).split(b'\0') if p}
trials = [
    ('lavapipe-intergroup-row-owner', '263abe74/lavapipe-intergroup-row-owner-finite-outcome',
     'a2477e97', 'development/lavapipe-intergroup-row-owner-2026-10-23'),
    ('native-cooperative-transport-intergroup', 'b1048e57/native-cooperative-transport-intergroup-finite-outcome',
     'a72c6082', 'development/native-cooperative-transport-intergroup-2026-10-23')]
for name, directory, preserved, branch in trials:
    W = R / 'out' / name; D = R / 'attestations/qualification' / directory
    assert not (D / 'cleanup.json').exists()
    remote = subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/heads/' + branch], text=True).split()[0]
    assert subprocess.run(['git', 'merge-base', '--is-ancestor', preserved, remote]).returncode == 0
    manifest = doc(D / 'MANIFEST.json'); closure = doc(D / 'ownership-closure.json')
    assert closure['all_owned_runs_terminal_and_reaped']
    for owner in closure['owners'].values():
        for birth in owner['observed_births']:
            p = Path('/proc') / str(birth['pid']) / 'stat'
            if p.exists(): assert p.read_text().rsplit(')', 1)[1].split()[19] != birth['saved_start_ticks']
    for key, expected in manifest['archives'].items(): assert pin(D / key) == expected
    targets = {}
    with zipfile.ZipFile(D / 'diagnostics-and-trial.zip') as z:
        for key, expected in manifest['diagnostics'].items():
            p = R / key
            if p in tracked or not p.is_file(): continue
            if not (p.is_relative_to(W / 'runtime') or p.is_relative_to(W / 'readers')): continue
            # Keep the concise mixed failure owner used by its saved debugger
            # script/source contract; its complete raw run remains archived.
            if name == 'native-cooperative-transport-intergroup' and p.name == 'owner.json': continue
            assert pin(p) == expected
            data = z.read(key); assert len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256']
            targets[p] = expected['bytes']
    # Only the previous all-row host products are unused. The mixed ELF,
    # archive and selector object are the concrete next codegen consumer.
    if name == 'lavapipe-intergroup-row-owner':
        for leaf in ['audited_backend_tests', 'libsirius_backend.a']:
            key = 'bin/linux-gcc/native-cpu-row-ownership/matched/' + leaf
            p = R / key; expected = manifest['actual_product_paths'][key]
            assert pin(p) == {k: expected[k] for k in ['bytes', 'sha256']}
            blob = manifest['content_addressed_blobs'][expected['sha256']]
            archive = R / blob['archive'] if blob.get('reused') else D / blob['archive']
            with zipfile.ZipFile(archive) as z:
                data = z.read(blob['member'])
                assert len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256']
            targets[p] = expected['bytes']
    receipt = D / 'publication-receipts.zip'
    assert pin(receipt) == doc(D / 'publication.json')['receipt_archive']
    with zipfile.ZipFile(receipt) as z:
        for key in z.namelist():
            p = R / key
            if p in tracked or not p.is_file() or not p.is_relative_to(W): continue
            data = z.read(key); assert data == p.read_bytes()
            targets[p] = len(data)
    for p in targets: p.unlink()
    for folder in [W / 'runtime', W / 'readers', R / 'bin/linux-gcc/native-cpu-row-ownership/matched']:
        if folder.exists():
            for p in sorted(folder.rglob('*'), key=lambda p: len(p.parts), reverse=True):
                if p.is_dir() and not any(p.iterdir()): p.rmdir()
            if not any(folder.iterdir()): folder.rmdir()
    result = {'preserved_and_pushed_before_deletion': preserved,
        'removed_files': len(targets), 'removed_bytes': sum(targets.values()),
        'all_original_owned_births_absent': True,
        'retained': ['All tracked source and pushed bundles', 'Actual source/build context for the scoped machine-codegen implementation',
            'Current provider and complete three-code observations for actual IR comparison',
            'Original complete-input products and mixed ELF/archive/selector object for the new provider consumer',
            'Exact E/DP module products and partial diagnostic evidence'],
        'scope': 'Only archive-verified completed numerical/performance raw runs, unused previous all-row host products and duplicate publication copies. No active build/source, current mixed products or code-proof dependency removed.'}
    (D / 'cleanup.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'trial': name, 'removed_files': len(targets), 'removed_bytes': sum(targets.values())}))
