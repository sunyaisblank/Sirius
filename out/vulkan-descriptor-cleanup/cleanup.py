from pathlib import Path
import hashlib, json, os, shutil, subprocess

R = Path.cwd()
C = Path(__file__).resolve().parent
W = R / 'out/vulkan-descriptor-reuse'
T = R / 'bin/linux-gcc/vulkan-descriptor-reuse'
A = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse'
D = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-cleanup-correction'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
load = lambda p: json.loads(p.read_text())

def repository_context():
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    assert revision == '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
    branch = subprocess.check_output(['git', 'branch', '--show-current'], text=True).strip()
    assert branch == 'development/committed-viewer-frame-2026-10-13'
    assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
    assert subprocess.check_output(['git', 'rev-list', '--left-right', '--count', 'HEAD...@{upstream}'], text=True).strip() == '0\t0'
    remote = subprocess.check_output(['git', 'ls-remote', '--heads', 'origin', branch], text=True).strip()
    assert remote == revision + '\trefs/heads/' + branch
    worktree = subprocess.check_output(['git', 'worktree', 'list', '--porcelain'], text=True)
    assert worktree.count('worktree ') == 1 and worktree.splitlines()[0] == 'worktree ' + str(R)
    restoration = load(A / 'files/out/vulkan-descriptor-reuse/restoration.json')
    for record in restoration['actual_primary_artifacts']:
        p = R / record['path']
        assert p.stat().st_size == record['bytes'] and sha(p) == record['sha256']
    pins = load(R / 'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/protected-before.json')
    for name, record in pins.items():
        p = R / name
        assert p.stat().st_size == record['bytes'] and sha(p) == record['sha256']
    renderer = restoration['renderer']
    p = R / renderer['path']
    assert p.stat().st_size == renderer['bytes'] and sha(p) == renderer['sha256']
    return {'revision': revision, 'branch': branch, 'remote': remote, 'one_canonical_worktree': True,
            'primary_artifact_bindings_exact': len(restoration['actual_primary_artifacts']),
            'protected_pins_exact': len(pins)}

def snapshot():
    result = {}
    for root in [W, T, C]:
        assert root.is_dir() and not root.is_symlink() and root.resolve() == root
        for p in [root, *root.rglob('*')]:
            assert not p.is_symlink(), str(p)
            stat = p.stat()
            assert p.is_file() or p.is_dir(), str(p)
            result[str(p)] = {'kind': 'file' if p.is_file() else 'directory',
                              'device': stat.st_dev, 'inode': stat.st_ino,
                              'mtime_ns': stat.st_mtime_ns}
    return result

assert C == R / 'out/vulkan-descriptor-cleanup'
assert C.is_dir() and not C.is_symlink() and C.resolve() == C
assert {p.name for p in C.iterdir()} == {'cleanup.py', 'source-preservation.json'}
assert D.is_dir() and not (D / 'manifest.json').exists()
preservation = load(C / 'source-preservation.json')
assert sha(C / 'cleanup.py') == preservation['source']['sha256']
assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + preservation['tag']], text=True).strip() == preservation['remote_readback']
assert (D / 'files/out/vulkan-descriptor-cleanup/cleanup.py').read_bytes() == (C / 'cleanup.py').read_bytes()
assert (D / 'files/out/vulkan-descriptor-cleanup/source-preservation.json').read_bytes() == (C / 'source-preservation.json').read_bytes()
context = repository_context()
initial = snapshot()
manifest = load(A / 'manifest.json')
expected = {r['path']: r for r in manifest['files']}
actual = {str(p.relative_to(A)): p for p in A.rglob('*') if p.is_file() and p != A / 'manifest.json'}
assert set(actual) == set(expected)
for key, p in actual.items():
    assert p.stat().st_size == expected[key]['bytes'] and sha(p) == expected[key]['sha256']
for p in W.rglob('*'):
    if p.is_file():
        target = A / 'files' / p.relative_to(R)
        assert target.read_bytes() == p.read_bytes(), str(p)
assert {str(p.relative_to(W)) for p in W.rglob('*') if p.is_file()} == {
    str(p.relative_to(A / 'files' / W.relative_to(R))) for p in
    (A / 'files' / W.relative_to(R)).rglob('*') if p.is_file()}
assert {str(p.relative_to(W)) for p in W.rglob('*') if p.is_dir()} == {
    str(p.relative_to(A / 'files' / W.relative_to(R))) for p in
    (A / 'files' / W.relative_to(R)).rglob('*') if p.is_dir()}
private = load(A / 'private-product-recovery.json')['records']
assert {r['origin'] for r in private} == {str(p.relative_to(R)) for p in T.rglob('*') if p.is_file()}
for record in private:
    live = R / record['origin']
    archived = R / record['recovery']['path']
    assert live.stat().st_size == record['bytes'] and sha(live) == record['sha256']
    assert archived.stat().st_size == record['bytes'] and sha(archived) == record['sha256']
    assert not archived.is_relative_to(W) and not archived.is_relative_to(T)
source = load(W / 'preserved-source.json')
assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + source['tag']], text=True).strip() == source['remote_readback']
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''

def proc(pid):
    words = Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()
    return {'pid': pid, 'ppid': int(words[1]), 'start_ticks': words[19]}

boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
births = set()
for p in W.rglob('owner.json'):
    owner = load(p)
    assert owner['status'] == 'terminal' and owner['remaining_owned_processes'] == []
    tokens = [(owner['boot_id'], int(birth['pid']), str(birth['saved_start_ticks']))
              for birth in owner['observed_births']]
    tokens.append((owner['boot_id'], int(owner['controller']['pid']),
                   str(owner['controller']['start_ticks'])))
    for token in tokens:
        births.add(token)
        if token[0] != boot:
            continue
        try:
            assert proc(token[1])['start_ticks'] != token[2], token
        except (FileNotFoundError, ProcessLookupError):
            pass
excluded = set()
pid = os.getpid()
while pid > 0 and pid not in excluded:
    excluded.add(pid)
    try:
        pid = proc(pid)['ppid']
    except (FileNotFoundError, ProcessLookupError):
        break
roots = [str(W), str(T), str(C)]
consumers = []
permission_limited = set()
for entry in Path('/proc').glob('[0-9]*'):
    pid = int(entry.name)
    if pid in excluded:
        continue
    try:
        identity = proc(pid)
        strings = [os.readlink(entry / 'cwd')]
        strings += [os.fsdecode(x) for x in (entry / 'cmdline').read_bytes().split(b'\0') if x]
        for fd in (entry / 'fd').iterdir():
            try:
                strings.append(os.readlink(fd))
            except (FileNotFoundError, ProcessLookupError):
                pass
        if any(root in value for root in roots for value in strings):
            assert proc(pid)['start_ticks'] == identity['start_ticks']
            consumers.append(identity)
    except PermissionError:
        permission_limited.add(pid)
    except (FileNotFoundError, ProcessLookupError):
        pass
assert not consumers, consumers
files = [p for root in [W, T, C] for p in root.rglob('*') if p.is_file()]
report = {'roots': roots, 'removed_files': len(files), 'removed_bytes': sum(p.stat().st_size for p in files),
          'sealed_manifest_sha256': sha(A / 'manifest.json'), 'saved_births_absent': len(births),
          'accessible_consumers': consumers, 'permission_limited_processes': len(permission_limited),
          'repository_context': context,
          'scope': 'Exact complete file/directory/root-identity checks, archive/private-root recovery, saved controller/child births and fresh scoped accessible-process check. Observer/ancestors excluded. No global hidden-process or scheduled-consumer absence claim.'}
# Recheck the complete mutable context and bytes immediately before unlinking.
assert snapshot() == initial
assert repository_context() == context
for p in W.rglob('*'):
    if p.is_file(): assert (A / 'files' / p.relative_to(R)).read_bytes() == p.read_bytes()
for record in private:
    p = R / record['origin']
    assert p.stat().st_size == record['bytes'] and sha(p) == record['sha256']
for key, p in actual.items():
    assert p.stat().st_size == expected[key]['bytes'] and sha(p) == expected[key]['sha256']
assert snapshot() == initial
assert {p.name for p in C.iterdir()} == {'cleanup.py', 'source-preservation.json'}
assert (D / 'files/out/vulkan-descriptor-cleanup/cleanup.py').read_bytes() == (C / 'cleanup.py').read_bytes()
assert (D / 'files/out/vulkan-descriptor-cleanup/source-preservation.json').read_bytes() == (C / 'source-preservation.json').read_bytes()
destination = D / 'cleanup-report.json'
assert not destination.exists()
shutil.rmtree(W)
shutil.rmtree(T)
destination.write_text(json.dumps(report, indent=2) + '\n')
supplement = [{'path': str(p.relative_to(D)), 'bytes': p.stat().st_size, 'sha256': sha(p)}
              for p in sorted(D.rglob('*')) if p.is_file()]
(D / 'manifest.json').write_text(json.dumps({'source_revision': context['revision'], 'files': supplement,
    'passed': True, 'scope': 'Separately pushed correction of the unexecuted cleanup helper; sealed original319-file evidence remains unchanged.'}, indent=2) + '\n')
(C / 'cleanup.py').unlink()
(C / 'source-preservation.json').unlink()
C.rmdir()
print(json.dumps(report))
