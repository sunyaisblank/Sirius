"""Retire archived rejected-trial raw files; retain the live compiler workspace."""
from pathlib import Path
import hashlib, json, os, subprocess, zipfile

R = Path.cwd(); W = R / 'out/lavapipe-acqrel-ordering'
D = R / 'attestations/qualification/1c95c98f/lavapipe-acqrel-ordering-finite-outcome'
PRESERVED = 'd324d909463615699517b9aa5034bb521ceddbd2'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
assert subprocess.run(['git', 'merge-base', '--is-ancestor', PRESERVED, 'HEAD']).returncode == 0
remote = subprocess.check_output(['git', 'ls-remote', 'origin',
    'refs/heads/development/lavapipe-acqrel-ordering-2026-10-23'], text=True).split()[0]
assert subprocess.run(['git', 'merge-base', '--is-ancestor', PRESERVED, remote]).returncode == 0
manifest = doc(D / 'MANIFEST.json'); closure = doc(D / 'ownership-closure.json')
assert closure['all_owned_runs_terminal_and_reaped'] and closure['owners_count'] == 9
for owner in closure['owners'].values():
    for birth in owner['observed_births']:
        p = Path('/proc') / str(birth['pid']) / 'stat'
        if p.exists(): assert p.read_text().rsplit(')', 1)[1].split()[19] != birth['saved_start_ticks']
for name, value in manifest['archives'].items(): assert pin(D / name) == value
tracked = {R / os.fsdecode(x) for x in subprocess.check_output(['git', 'ls-files', '-z']).split(b'\0') if x}
targets = {}
with zipfile.ZipFile(D / 'diagnostics-and-trial.zip') as z:
    for key, value in manifest['diagnostics'].items():
        p = R / key
        if p in tracked or not p.is_relative_to(W): continue
        assert not p.is_relative_to(W / 'upstream')
        assert pin(p) == value
        data = z.read(key)
        assert len(data) == value['bytes'] and hashlib.sha256(data).hexdigest() == value['sha256']
        targets[p] = value['bytes']
publication = W / 'publication-api.json'
if publication.exists():
    with zipfile.ZipFile(D / 'publication-receipts.zip') as z:
        assert z.read(str(publication.relative_to(R))) == publication.read_bytes()
    targets[publication] = publication.stat().st_size
for p in targets: p.unlink()
for p in sorted(W.rglob('*'), key=lambda p: len(p.parts), reverse=True):
    if p.is_dir() and not p.is_relative_to(W / 'upstream') and not any(p.iterdir()): p.rmdir()
result = {'preservation_commit': PRESERVED, 'pushed_before_deletion': True,
    'removed_files': len(targets), 'removed_bytes': sum(targets.values()),
    'all9_original_owned_runs_terminal_and_births_absent': True,
    'retained_active_source': str((W / 'upstream').relative_to(R)),
    'retained_active_build': 'bin/linux-gcc/lavapipe-acqrel-ordering',
    'actual_consumer': 'Distinct intergroup row-owner source implementation and subsequent incremental build; previous acquire-release modification is restored before composition.',
    'scope': 'Only archive-verified rejected-trial untracked diagnostics/publication copy removed. Active source/build context, tracked sources, original providers, restored exact row-owner products and all pushed bundles retained.'}
p = D / 'cleanup.json'; assert not p.exists(); p.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
