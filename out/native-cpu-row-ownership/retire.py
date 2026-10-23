"""Remove this rejected candidate's owned scratch only after pushed preservation."""
from pathlib import Path
import hashlib, json, os, subprocess, zipfile

R = Path.cwd(); W = Path(__file__).resolve().parent
D = R / 'attestations/qualification/bfa095ba/native-cpu-row-ownership-finite-outcome'
PRESERVED = '9e01dfdca7798c79927827334798c72c2774eebc'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == PRESERVED
assert subprocess.check_output(['git', 'rev-parse', '@{upstream}'], text=True).strip() == PRESERVED
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
manifest = doc(D / 'MANIFEST.json'); closure = doc(D / 'ownership-closure.json')
assert closure['all_owned_runs_terminal_and_reaped'] and closure['owners_count'] == 10
for owner in closure['owners'].values():
    for birth in owner['observed_births']:
        p = Path('/proc') / str(birth['pid']) / 'stat'
        if p.exists(): assert p.read_text().rsplit(')', 1)[1].split()[19] != birth['saved_start_ticks']
root = R / 'bin/linux-gcc/native-cpu-row-ownership'
needle = str(root).encode(); inaccessible = []; consumers = []
for proc in Path('/proc').glob('[0-9]*'):
    try:
        if needle in (proc / 'cmdline').read_bytes() or needle in (proc / 'maps').read_bytes():
            consumers.append(int(proc.name))
    except (FileNotFoundError, ProcessLookupError): pass
    except PermissionError: inaccessible.append(int(proc.name))
assert not consumers, consumers
for name, value in manifest['archives'].items(): assert pin(D / name) == value
for name, value in manifest['reused_archives'].items(): assert pin(R / name) == value['identity']
tracked = {R / os.fsdecode(x) for x in subprocess.check_output(['git', 'ls-files', '-z']).split(b'\0') if x}
targets = {}
with zipfile.ZipFile(D / 'diagnostics-and-trial.zip') as z:
    for key, value in manifest['diagnostics'].items():
        p = R / key
        if p in tracked: continue
        assert p.is_relative_to(W) and p.is_file()
        assert pin(p) == value
        data = z.read(key)
        assert len(data) == value['bytes'] and hashlib.sha256(data).hexdigest() == value['sha256']
        targets[p] = value['bytes']
for p in root.rglob('*'):
    if not p.is_file(): continue
    assert p not in tracked
    key = str(p.relative_to(R)); value = manifest['actual_product_paths'][key]
    assert pin(p) == {k: value[k] for k in ['bytes', 'sha256']}
    blob = manifest['content_addressed_blobs'][value['sha256']]
    archive = (R if blob['reused'] else D) / blob['archive']
    with zipfile.ZipFile(archive) as z:
        data = z.read(blob['member'])
        assert len(data) == value['bytes'] and hashlib.sha256(data).hexdigest() == value['sha256']
    targets[p] = value['bytes']
for p in sorted(targets, key=lambda p: len(p.parts), reverse=True): p.unlink()
for parent in [root, W]:
    for p in sorted(parent.rglob('*'), key=lambda p: len(p.parts), reverse=True):
        if p.is_dir() and not any(p.iterdir()): p.rmdir()
    if parent == root and parent.exists() and not any(parent.iterdir()): parent.rmdir()
result = {'preservation_commit': PRESERVED, 'pushed_before_deletion': True,
    'removed_files': len(targets), 'removed_bytes': sum(targets.values()),
    'retained_tracked_source_files': sum(p.is_relative_to(W) for p in tracked),
    'all10_original_owned_runs_terminal_and_births_absent': True,
    'readable_process_consumers': consumers, 'unreadable_other_process_ids': inaccessible,
    'scope': 'Only this rejected candidate archive-verified untracked diagnostics and actual product root removed. Existing private IEEE32 providers, previous admitted complete-input products, public sources, canonical products and all pushed bundles retained. Other inaccessible process state is not claimed universally visible.'}
p = D / 'cleanup.json'; assert not p.exists()
p.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
