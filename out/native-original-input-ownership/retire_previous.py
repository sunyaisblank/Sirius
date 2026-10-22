"""Retire only pushed, hash-verified previous scratch absent from current inputs."""
from pathlib import Path
import hashlib, json, os, subprocess
R = Path.cwd(); D = R / 'attestations/qualification/d9dbd657/native-phase-sharing-finite-pass-frame-timeout'
M = json.loads((D / 'MANIFEST.json').read_text())
subprocess.run(['git', 'merge-base', '--is-ancestor', '1bdfb8a8e8bb5a0865e1913c7febf483db48fe35',
                'origin/development/native-original-input-ownership-2026-10-23'], check=True)
tracked = set(subprocess.check_output(['git', 'ls-files', '-z']).decode().split('\0'))
active = set(tracked)
W = R / 'out/native-original-input-ownership'
for p in [W / 'actions.json', W / 'ordinary/actions.json', W / 'runtime/actions.json']:
    for a in json.loads(p.read_text()).values(): active.update(a['frozen_inputs'])
roots = [R / 'out/native-endpoint-phase-sharing', R / 'bin/linux-gcc/native-endpoint-phase-sharing',
         R / 'renders/native-endpoint-phase-d9dbd657']
records = dict(M['diagnostics']); records.update(M['actual_product_paths'])
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
candidates = {}
for root in roots:
    for p in root.rglob('*'):
        if not (p.is_file() or p.is_symlink()): continue
        key = str(p.relative_to(R))
        if key in tracked: continue
        assert key not in active, key
        assert key in records, ('unpreserved path', key)
        assert pin(p) == {k: records[key][k] for k in ['bytes', 'sha256']}, key
        if p.is_symlink(): assert os.readlink(p) == records[key]['symlink']
        candidates[key] = pin(p)
# Verify all original owned births are closed. Other readable processes supply
# an additional positive-consumer check, not universal process visibility.
closure = json.loads((D / 'ownership-closure.json').read_text())
for owner in closure['owners'].values():
    for birth in owner['observed_births']:
        p = Path('/proc') / str(birth['pid']) / 'stat'
        if p.exists(): assert p.read_text().rsplit(')', 1)[1].split()[19] != birth['saved_start_ticks']
observed = 0; inaccessible = []
for proc in Path('/proc').glob('[0-9]*'):
    try:
        status = (proc / 'status').read_text()
        uid = next(int(x.split()[1]) for x in status.splitlines() if x.startswith('Uid:'))
        if uid != os.getuid(): continue
        if proc.name == str(os.getpid()): continue
        maps = (proc / 'maps').read_text()
        command = (proc / 'cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace')
        cwd = os.readlink(proc / 'cwd')
        assert not any(str(root) in maps or str(root) in command or cwd.startswith(str(root)) for root in roots), ('active consumer', proc.name)
        observed += 1
    except (FileNotFoundError, ProcessLookupError): pass
    except PermissionError: inaccessible.append(int(proc.name))
for key in candidates: (R / key).unlink()
for root in roots:
    for p in sorted((p for p in root.rglob('*') if p.is_dir()), key=lambda p: len(p.parts), reverse=True):
        if not any(p.iterdir()): p.rmdir()
    if root.exists() and not any(root.iterdir()): root.rmdir()
assert subprocess.check_output(['git', 'status', '--porcelain']).strip() == b''
out = {'preservation': '1bdfb8a8e8bb5a0865e1913c7febf483db48fe35', 'archived_paths_verified_then_retired': candidates,
    'paths_retired': len(candidates), 'bytes_retired': sum(x['bytes'] for x in candidates.values()),
    'same_user_processes_checked': observed, 'current_explicit_and_tracked_inputs_preserved': True,
    'original17_owned_births_closed': True, 'other_processes_inaccessible': inaccessible,
    'process_coverage_scope': 'Original owned births plus readable same-user maps/commands/cwd; no universal process visibility claim.',
    'current_provider_Dopri_and_new_candidate_products_retained': True,
    'recovery': 'Pushed MANIFEST selects authentic diagnostic paths and content-addressed actual product blobs.'}
p = D / 'cleanup.json'; assert not p.exists(); p.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: out[k] for k in ['paths_retired', 'bytes_retired', 'current_explicit_and_tracked_inputs_preserved']}))
