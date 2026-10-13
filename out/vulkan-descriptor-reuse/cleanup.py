from pathlib import Path
import hashlib, json, os, shutil, subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
T = R / 'bin/linux-gcc/vulkan-descriptor-reuse'
A = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
load = lambda p: json.loads(p.read_text())
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
    for birth in owner['observed_births']:
        token = (owner['boot_id'], int(birth['pid']), str(birth['start_ticks']))
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
roots = [str(W), str(T)]
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
files = [p for root in [W, T] for p in root.rglob('*') if p.is_file()]
report = {'roots': roots, 'removed_files': len(files), 'removed_bytes': sum(p.stat().st_size for p in files),
          'sealed_manifest_sha256': sha(A / 'manifest.json'), 'saved_births_absent': len(births),
          'accessible_consumers': consumers, 'permission_limited_processes': len(permission_limited),
          'scope': 'Exact archive/private-root recovery coverage and fresh scoped accessible-process check. Observer/ancestors excluded. No global hidden-process or scheduled-consumer absence claim.'}
destination = A.parent / 'vulkan-descriptor-reuse-cleanup.json'
assert not destination.exists()
shutil.rmtree(W)
shutil.rmtree(T)
destination.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
